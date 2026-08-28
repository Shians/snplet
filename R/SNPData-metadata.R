#' Add metadata to a SNPData object's info tables
#'
#' These functions provide a standardised interface for adding new columns (and, for
#' \code{add_donor_snp_metadata}, new rows) to the \code{barcode_info}, \code{snp_info},
#' \code{donor_info}, or \code{donor_snp_info} tables of a SNPData object. The functions
#' ensure data integrity by validating dimensions and preserving automatically
#' computed summary statistics.
#'
#' @param x A SNPData object, required.
#' @param metadata A data.frame, required, with new columns (and, for
#'   \code{add_donor_snp_metadata}, possibly new rows) to add.
#' @param join_by Character vector; default varies by function. Column(s) to
#'   join by: \code{"cell_id"} (default) or \code{"barcode"} for
#'   \code{add_barcode_metadata}; \code{"snp_id"} (default) for
#'   \code{add_snp_metadata}; \code{"donor"} (default) for
#'   \code{add_donor_metadata}; \code{c("snp_id", "donor", "zygosity_source")}
#'   (default) for \code{add_donor_snp_metadata}.
#' @param overwrite Logical (default \code{FALSE}, except \code{TRUE} for
#'   \code{add_donor_snp_metadata}, since \code{donor_snp_info}'s columns are
#'   fixed by its schema from construction, so treating them as pre-existing
#'   and requiring \code{overwrite = FALSE} would reject on every call).
#'   Whether to overwrite existing columns.
#'
#' @return A SNPData object with the updated table
#'
#' @examples
#' \dontrun{
#' snpdata <- get_example_snpdata()
#'
#' # Add new columns to barcode_info
#' new_barcode_info <- data.frame(
#'   cell_id = c("cell_1", "cell_2"),
#'   cell_type = c("T_cell", "B_cell"),
#'   treatment = c("control", "treated")
#' )
#' updated_snpdata <- add_barcode_metadata(snpdata, new_barcode_info)
#'
#' # Update existing columns in barcode_info by setting overwrite=TRUE
#' clonotype_info <- data.frame(
#'   cell_id = c("cell_1", "cell_2"),
#'   clonotype = c("clonotype_1", "clonotype_2")
#' )
#' updated_snpdata <- add_barcode_metadata(
#'   snpdata,
#'   clonotype_info,
#'   overwrite = TRUE
#' )
#'
#' # Add columns to barcode_info for subset of cells
#' # (cells not in metadata will have NA for new columns)
#' partial_barcode_info <- data.frame(
#'   cell_id = c("cell_1", "cell_2"),
#'   annotation = c("A", "B")
#' )
#' updated_snpdata <- add_barcode_metadata(snpdata, partial_barcode_info)
#'
#' # Add new columns to snp_info
#' new_snp_info <- data.frame(
#'   snp_id = c("snp_1", "snp_2"),
#'   gene_name = c("GENE1", "GENE2"),
#'   gene_biotype = c("protein_coding", "lncRNA")
#' )
#' updated_snpdata <- add_snp_metadata(snpdata, new_snp_info)
#' }
#'
#' @name add_metadata
NULL

.column_word <- function(join_by) if (length(join_by) == 1) "Column" else "Column(s)"

.update_metadata_info <- function(
    current_info,
    metadata,
    join_by,
    overwrite,
    preserved_cols,
    info_name,
    join_type = c("left", "full")
) {
    join_type <- match.arg(join_type)
    if (!all(join_by %in% colnames(current_info))) {
        stop(paste0(
            .column_word(join_by),
            " '",
            paste(join_by, collapse = ", "),
            "' not found in current ",
            info_name
        ))
    }

    if (any(duplicated(metadata[join_by]))) {
        stop(paste0(
            "Duplicate values found in join ",
            tolower(.column_word(join_by)),
            " '",
            paste(join_by, collapse = ", "),
            "' of metadata data.frame"
        ))
    }

    # The same check on the receiving side. Without it a key that repeats in
    # `current_info` silently copies one metadata row onto every row sharing that
    # key -- the row count is unchanged, so nothing downstream catches it. This
    # bites on `barcode`, which merge_snpdata() only guarantees unique within a
    # library, so a merged object can carry one barcode on two different cells.
    if (any(duplicated(current_info[join_by]))) {
        hint <- if ("barcode" %in% join_by) {
            paste0(
                " A barcode identifies a cell only within its own library, so a merged object may repeat one.",
                " Join by 'cell_id', or by c(\"library_id\", \"barcode\")."
            )
        } else {
            ""
        }
        stop(paste0(
            "Duplicate values found in join ",
            tolower(.column_word(join_by)),
            " '",
            paste(join_by, collapse = ", "),
            "' of ",
            info_name,
            "; each row must match at most one metadata row.",
            hint
        ))
    }

    if (!overwrite) {
        conflicting_cols <- intersect(
            setdiff(colnames(metadata), join_by),
            colnames(current_info)
        )
        if (length(conflicting_cols) > 0) {
            stop(paste0(
                "Column(s) already exist in ",
                info_name,
                ": ",
                paste(conflicting_cols, collapse = ", "),
                ". Set overwrite=TRUE to replace existing columns."
            ))
        }
    }

    join_fun <- if (join_type == "full") dplyr::full_join else dplyr::left_join
    updated_info <- current_info %>%
        join_fun(metadata, by = join_by, suffix = c("", ".new"))

    if (overwrite && any(grepl("\\.new$", colnames(updated_info)))) {
        new_cols <- grep("\\.new$", colnames(updated_info), value = TRUE)
        base_cols <- sub("\\.new$", "", new_cols)

        # A row absent from `metadata` produces NA in every ".new" column from
        # the join -- indistinguishable, by value alone, from metadata
        # legitimately setting that column to NA for a row it does cover. So
        # "was this row actually in metadata" has to be tracked directly
        # rather than inferred from NA-ness: only matched rows get
        # overwritten (NA included), unmatched rows are left untouched.
        key_of <- function(df, cols) do.call(paste, c(as.list(df[cols]), list(sep = "")))
        matched <- key_of(updated_info, join_by) %in% key_of(metadata, join_by)

        for (i in seq_along(new_cols)) {
            updated_info[[base_cols[i]]][matched] <- updated_info[[new_cols[i]]][matched]
            updated_info[[new_cols[i]]] <- NULL
        }
    }

    auto_computed <- current_info[preserved_cols[preserved_cols %in% colnames(current_info)]]

    for (col in preserved_cols) {
        if (col %in% colnames(updated_info)) {
            updated_info[[col]] <- NULL
        }
    }

    for (col in names(auto_computed)) {
        updated_info[[col]] <- auto_computed[[col]]
    }

    updated_info
}

#' Rebuild a SNPData object with updated metadata table(s) only
#'
#' The \code{add_*_metadata()} family only ever changes annotation columns
#' (and, for \code{donor_snp_info}, adds rows for new (snp_id, donor,
#' zygosity_source) triples); the count matrices and the identity/order of
#' SNPs and cells never change. Routing that through \code{new("SNPData",
#' ...)} would re-run \code{.recompute_metrics()}, which adds the full
#' \code{ref_count}/\code{alt_count} sparse matrices and takes
#' \code{rowSums}/\code{colSums} across the entire dataset, purely to
#' recompute stats that are unaffected by a metadata-only change. This
#' instead patches the requested slot(s) directly and revalidates, skipping
#' that recomputation entirely.
#'
#' @keywords internal
.clone_with_metadata <- function(x, snp_info = NULL, barcode_info = NULL, donor_info = NULL, donor_snp_info = NULL) {
    result <- x
    if (!is.null(snp_info)) {
        result@snp_info <- tibble::as_tibble(snp_info)
    }
    if (!is.null(barcode_info)) {
        result@barcode_info <- tibble::as_tibble(barcode_info)
        result <- .resync_library_info(result, result@barcode_info)
    }
    if (!is.null(donor_info)) {
        result@donor_info <- tibble::as_tibble(donor_info)
    }
    if (!is.null(donor_snp_info)) {
        result@donor_snp_info <- tibble::as_tibble(donor_snp_info)
    }
    methods::validObject(result)
    result
}

#' @rdname add_metadata
#' @export
add_barcode_metadata <- function(x, metadata, join_by = "cell_id", overwrite = FALSE) {
    # Validate input
    if (!methods::is(x, "SNPData")) {
        stop("Input must be a SNPData object")
    }

    if (!is.data.frame(metadata)) {
        stop("metadata must be a data.frame")
    }

    if (!all(join_by %in% colnames(metadata))) {
        stop(paste0(.column_word(join_by), " '", paste(join_by, collapse = ", "), "' not found in metadata data.frame"))
    }

    current_barcode_info <- barcode_info(x)
    updated_barcode_info <- .update_metadata_info(
        current_info = current_barcode_info,
        metadata = metadata,
        join_by = join_by,
        overwrite = overwrite,
        preserved_cols = c("library_size", "non_zero_snps"),
        info_name = "barcode_info"
    )

    result <- .clone_with_metadata(x, barcode_info = updated_barcode_info)
    .propagate_zygosity_source(result, x)
}

#' @rdname add_metadata
#' @export
add_snp_metadata <- function(x, metadata, join_by = "snp_id", overwrite = FALSE) {
    # Validate input
    if (!methods::is(x, "SNPData")) {
        stop("Input must be a SNPData object")
    }

    if (!is.data.frame(metadata)) {
        stop("metadata must be a data.frame")
    }

    if (!all(join_by %in% colnames(metadata))) {
        stop(paste0(.column_word(join_by), " '", paste(join_by, collapse = ", "), "' not found in metadata data.frame"))
    }

    current_snp_info <- snp_info(x)
    updated_snp_info <- .update_metadata_info(
        current_info = current_snp_info,
        metadata = metadata,
        join_by = join_by,
        overwrite = overwrite,
        preserved_cols = c("coverage", "non_zero_samples"),
        info_name = "snp_info"
    )

    result <- .clone_with_metadata(x, snp_info = updated_snp_info)
    .propagate_zygosity_source(result, x)
}

#' @rdname add_metadata
#' @export
add_donor_metadata <- function(x, metadata, join_by = "donor", overwrite = FALSE) {
    # Validate input
    if (!methods::is(x, "SNPData")) {
        stop("Input must be a SNPData object")
    }

    if (!is.data.frame(metadata)) {
        stop("metadata must be a data.frame")
    }

    if (!all(join_by %in% colnames(metadata))) {
        stop(paste0(.column_word(join_by), " '", paste(join_by, collapse = ", "), "' not found in metadata data.frame"))
    }

    current_donor_info <- donor_info(x)
    updated_donor_info <- .update_metadata_info(
        current_info = current_donor_info,
        metadata = metadata,
        join_by = join_by,
        overwrite = overwrite,
        preserved_cols = c("n_cells"),
        info_name = "donor_info"
    )

    result <- .clone_with_metadata(x, donor_info = updated_donor_info)
    .propagate_zygosity_source(result, x)
}

#' @rdname add_metadata
#' @export
add_donor_snp_metadata <- function(x, metadata, join_by = c("snp_id", "donor", "zygosity_source"), overwrite = TRUE) {
    # Validate input
    if (!methods::is(x, "SNPData")) {
        stop("Input must be a SNPData object")
    }

    if (!is.data.frame(metadata)) {
        stop("metadata must be a data.frame")
    }

    if (!all(join_by %in% colnames(metadata))) {
        stop(paste0(.column_word(join_by), " '", paste(join_by, collapse = ", "), "' not found in metadata data.frame"))
    }

    current_donor_snp_info <- donor_snp_info(x, source = "all")
    updated_donor_snp_info <- .update_metadata_info(
        current_info = current_donor_snp_info,
        metadata = metadata,
        join_by = join_by,
        overwrite = overwrite,
        preserved_cols = character(0),
        info_name = "donor_snp_info",
        join_type = "full"
    )

    result <- .clone_with_metadata(x, donor_snp_info = updated_donor_snp_info)
    .propagate_zygosity_source(result, x)
}

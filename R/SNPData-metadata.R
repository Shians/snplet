#' Add metadata to a SNPData object's info tables
#'
#' These functions provide a standardized interface for adding new columns (and, for
#' \code{add_donor_snp_metadata}, new rows) to the \code{barcode_info}, \code{snp_info},
#' \code{donor_info}, or \code{donor_snp_info} tables of a SNPData object. The functions
#' ensure data integrity by validating dimensions and preserving automatically
#' computed summary statistics.
#'
#' @param x A SNPData object
#' @param metadata A data.frame containing new columns (and, for \code{add_donor_snp_metadata},
#'   possibly new rows) to add
#' @param join_by Character vector specifying the column(s) to join by.
#'   For \code{barcode_info}: \code{"cell_id"} or \code{"barcode"}. For \code{snp_info}:
#'   \code{"snp_id"}. For \code{donor_info}: \code{"donor"}. For \code{donor_snp_info}:
#'   \code{c("snp_id", "donor")}.
#' @param overwrite Logical, whether to overwrite existing columns (default \code{FALSE},
#'   except for \code{add_donor_snp_metadata} where it defaults to \code{TRUE} since
#'   \code{donor_snp_info}'s columns are fixed by its schema from construction, so treating
#'   them as pre-existing and requiring \code{overwrite=FALSE} would reject on every call).
#'   Set to \code{TRUE} to update existing columns.
#'
#' @return A SNPData object with the updated table
#'
#' @examples
#' \dontrun{
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

        for (i in seq_along(new_cols)) {
            updated_info[[base_cols[i]]] <- updated_info[[new_cols[i]]]
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

    current_barcode_info <- get_barcode_info(x)
    updated_barcode_info <- .update_metadata_info(
        current_info = current_barcode_info,
        metadata = metadata,
        join_by = join_by,
        overwrite = overwrite,
        preserved_cols = c("library_size", "non_zero_snps"),
        info_name = "barcode_info"
    )

    # Create new SNPData object
    new(
        "SNPData",
        ref_count = x@ref_count,
        alt_count = x@alt_count,
        oth_count = x@oth_count,
        snp_info = x@snp_info,
        barcode_info = updated_barcode_info,
        donor_info = get_donor_info(x),
        donor_snp_info = get_donor_snp_info(x)
    )
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

    current_snp_info <- get_snp_info(x)
    updated_snp_info <- .update_metadata_info(
        current_info = current_snp_info,
        metadata = metadata,
        join_by = join_by,
        overwrite = overwrite,
        preserved_cols = c("coverage", "non_zero_samples"),
        info_name = "snp_info"
    )

    # Create new SNPData object
    new(
        "SNPData",
        ref_count = x@ref_count,
        alt_count = x@alt_count,
        oth_count = x@oth_count,
        snp_info = updated_snp_info,
        barcode_info = x@barcode_info,
        donor_info = get_donor_info(x),
        donor_snp_info = get_donor_snp_info(x)
    )
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

    current_donor_info <- get_donor_info(x)
    updated_donor_info <- .update_metadata_info(
        current_info = current_donor_info,
        metadata = metadata,
        join_by = join_by,
        overwrite = overwrite,
        preserved_cols = c("n_cells"),
        info_name = "donor_info"
    )

    # Create new SNPData object
    new(
        "SNPData",
        ref_count = x@ref_count,
        alt_count = x@alt_count,
        oth_count = x@oth_count,
        snp_info = x@snp_info,
        barcode_info = x@barcode_info,
        donor_info = updated_donor_info,
        donor_snp_info = get_donor_snp_info(x)
    )
}

#' @rdname add_metadata
#' @export
add_donor_snp_metadata <- function(x, metadata, join_by = c("snp_id", "donor"), overwrite = TRUE) {
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

    current_donor_snp_info <- get_donor_snp_info(x)
    updated_donor_snp_info <- .update_metadata_info(
        current_info = current_donor_snp_info,
        metadata = metadata,
        join_by = join_by,
        overwrite = overwrite,
        preserved_cols = character(0),
        info_name = "donor_snp_info",
        join_type = "full"
    )

    # Create new SNPData object
    new(
        "SNPData",
        ref_count = x@ref_count,
        alt_count = x@alt_count,
        oth_count = x@oth_count,
        snp_info = x@snp_info,
        barcode_info = x@barcode_info,
        donor_info = get_donor_info(x),
        donor_snp_info = updated_donor_snp_info
    )
}

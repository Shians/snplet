#' Add metadata to SNPData barcode_info or snp_info
#'
#' These functions provide a standardized interface for adding new columns
#' to the barcode_info or snp_info data frames of a SNPData object. The functions
#' ensure data integrity by validating dimensions and preserving automatically
#' computed summary statistics.
#'
#' @param x A SNPData object
#' @param metadata A data.frame containing new columns to add to barcode_info or snp_info
#' @param join_by Character string specifying the column name to join by.
#'   For barcode_info: "cell_id" or "barcode"
#'   For snp_info: "snp_id"
#' @param overwrite Logical, whether to overwrite existing columns (default FALSE).
#'   Set to TRUE to update existing columns in barcode_info or snp_info.
#'
#' @return A SNPData object with updated barcode_info or snp_info
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

.update_metadata_info <- function(
    current_info,
    metadata,
    join_by,
    overwrite,
    preserved_cols,
    info_name
) {
    if (!join_by %in% colnames(current_info)) {
        stop(paste0("Column '", join_by, "' not found in current ", info_name))
    }

    if (any(duplicated(metadata[[join_by]]))) {
        stop(paste0(
            "Duplicate values found in join column '",
            join_by,
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

    updated_info <- current_info %>%
        dplyr::left_join(metadata, by = join_by, suffix = c("", ".new"))

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

    if (!join_by %in% colnames(metadata)) {
        stop(paste0("Column '", join_by, "' not found in metadata data.frame"))
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
        barcode_info = updated_barcode_info
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

    if (!join_by %in% colnames(metadata)) {
        stop(paste0("Column '", join_by, "' not found in metadata data.frame"))
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
        barcode_info = x@barcode_info
    )
}

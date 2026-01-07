#' Merge two SNPData objects
#'
#' Combines two SNPData objects with flexible join strategies for SNPs and cells.
#' Overlapping entries are summed; unique entries are retained based on join type.
#'
#' @param x A SNPData object
#' @param y A SNPData object to merge with x
#' @param snp_join Join type for SNPs (rows). One of:
#'   \describe{
#'     \item{\code{"union"}}{Keep all SNPs from both objects (default)}
#'     \item{\code{"intersect"}}{Keep only SNPs present in both objects}
#'     \item{\code{"left"}}{Keep all SNPs from x, add y data where available}
#'     \item{\code{"right"}}{Keep all SNPs from y, add x data where available}
#'   }
#' @param cell_join Join type for cells (columns). One of:
#'   \describe{
#'     \item{\code{"union"}}{Keep all cells from both objects (default)}
#'     \item{\code{"intersect"}}{Keep only cells present in both objects}
#'     \item{\code{"left"}}{Keep all cells from x, add y data where available}
#'     \item{\code{"right"}}{Keep all cells from y, add x data where available}
#'   }
#'
#' @return A new SNPData object containing the merged data
#'
#' @details
#' This function provides independent control over which SNPs (rows) and cells
#' (columns) to retain when merging two SNPData objects. For positions present
#' in both objects, counts are summed element-wise. For positions present in
#' only one object (but retained by the join strategy), counts are kept with
#' zero-filling for the missing dimension.
#'
#' **Cell Merging Strategy:**
#' When both SNPData objects contain a \code{barcode} column in their
#' \code{barcode_info}, cells are matched and merged by their cell barcodes
#' rather than by \code{cell_id}. This allows proper merging of data from the
#' same physical cells across datasets, even if they have different internal
#' \code{cell_id} values. For cells with the same barcode, counts are summed.
#' If the \code{barcode} column is absent from either object, the function falls
#' back to merging by \code{cell_id}.
#'
#' Metadata is merged using the corresponding dplyr join function. For overlapping
#' SNPs or cells with conflicting metadata, values from x take priority. Auto-
#' computed columns (coverage, non_zero_samples, library_size, non_zero_snps)
#' are recalculated for the merged object.
#'
#' The independent join parameters enable fine-grained control:
#' \itemize{
#'   \item \code{snp_join="union", cell_join="union"}: Maximum data retention
#'   \item \code{snp_join="intersect", cell_join="intersect"}: Strict QC, only validated entries
#'   \item \code{snp_join="union", cell_join="intersect"}: Track more SNPs in same cells
#'   \item \code{snp_join="intersect", cell_join="union"}: Track consensus SNPs across more cells
#'   \item \code{snp_join="left", cell_join="left"}: Augment primary dataset
#' }
#'
#' @examples
#' \dontrun{
#' # Technical replicates - sum depth for validated entries only
#' rep1 <- import_cellsnp("replicate1/")
#' rep2 <- import_cellsnp("replicate2/")
#' combined <- merge_snpdata(rep1, rep2,
#'                          snp_join = "intersect",
#'                          cell_join = "intersect")
#'
#' # Batch integration - all cells, consensus SNPs only
#' batch1 <- import_cellsnp("donor1/")
#' batch2 <- import_cellsnp("donor2/")
#' integrated <- merge_snpdata(batch1, batch2,
#'                             snp_join = "intersect",
#'                             cell_join = "union")
#'
#' # SNP panel expansion - all SNPs, validated cells only
#' common <- import_cellsnp("common_vars/")
#' rare <- import_cellsnp("rare_vars/")
#' expanded <- merge_snpdata(common, rare,
#'                          snp_join = "union",
#'                          cell_join = "intersect")
#'
#' # Default: maximum data retention
#' dataset1 <- import_cellsnp("cohort_A/")
#' dataset2 <- import_cellsnp("cohort_B/")
#' combined <- merge_snpdata(dataset1, dataset2)
#' }
#'
#' @export
merge_snpdata <- function(
    x,
    y,
    snp_join = c("union", "intersect", "left", "right"),
    cell_join = c("union", "intersect", "left", "right")
) {
    # Match arguments
    snp_join <- match.arg(snp_join)
    cell_join <- match.arg(cell_join)

    # Validate both inputs are SNPData objects
    stopifnot(methods::is(x, "SNPData"))
    stopifnot(methods::is(y, "SNPData"))

    # Check for required identifier columns
    stopifnot("snp_id" %in% colnames(x@snp_info))
    stopifnot("snp_id" %in% colnames(y@snp_info))
    stopifnot("cell_id" %in% colnames(x@barcode_info))
    stopifnot("cell_id" %in% colnames(y@barcode_info))

    # Extract identifiers
    snp_ids_x <- rownames(x@ref_count)
    snp_ids_y <- rownames(y@ref_count)
    cell_ids_x <- colnames(x@ref_count)
    cell_ids_y <- colnames(y@ref_count)

    # Check if barcode column exists in both objects
    has_barcode_x <- "barcode" %in% colnames(x@barcode_info)
    has_barcode_y <- "barcode" %in% colnames(y@barcode_info)

    # Use barcodes for merging if available, otherwise fall back to cell_ids
    if (has_barcode_x && has_barcode_y) {
        # Extract barcodes from barcode_info
        barcodes_x <- x@barcode_info$barcode
        barcodes_y <- y@barcode_info$barcode

        # Create mapping from barcode to original cell_id (colnames)
        barcode_to_colname_x <- setNames(cell_ids_x, barcodes_x)
        barcode_to_colname_y <- setNames(cell_ids_y, barcodes_y)

        # Determine which barcodes to retain based on cell_join strategy
        barcodes_retained <- switch(
            cell_join,
            "union" = union(barcodes_x, barcodes_y),
            "intersect" = intersect(barcodes_x, barcodes_y),
            "left" = barcodes_x,
            "right" = barcodes_y
        )

        # Create mapping from original column names to retained barcode positions
        # This will be used to reorganize the matrices
        col_mapping_x <- match(barcodes_x, barcodes_retained)
        col_mapping_y <- match(barcodes_y, barcodes_retained)

        # Create new cell_ids for the retained barcodes
        cell_ids_retained <- paste0("cell_", seq_along(barcodes_retained))
    } else {
        # Fall back to using cell_ids from column names
        cell_ids_retained <- switch(
            cell_join,
            "union" = union(cell_ids_x, cell_ids_y),
            "intersect" = intersect(cell_ids_x, cell_ids_y),
            "left" = cell_ids_x,
            "right" = cell_ids_y
        )

        barcodes_retained <- cell_ids_retained
        col_mapping_x <- NULL
        col_mapping_y <- NULL
    }

    # Determine which SNPs to retain based on snp_join strategy
    snp_ids_retained <- switch(
        snp_join,
        "union" = union(snp_ids_x, snp_ids_y),
        "intersect" = intersect(snp_ids_x, snp_ids_y),
        "left" = snp_ids_x,
        "right" = snp_ids_y
    )

    # Handle edge case: empty result
    if (length(snp_ids_retained) == 0) {
        stop("No SNPs to retain after merge. Check your join strategy and input data.")
    }
    if (length(cell_ids_retained) == 0) {
        stop("No cells to retain after merge. Check your join strategy and input data.")
    }

    # Identify which entries need data from each object
    snp_from_x <- intersect(snp_ids_retained, snp_ids_x)
    snp_from_y <- intersect(snp_ids_retained, snp_ids_y)
    cell_from_x <- intersect(cell_ids_retained, cell_ids_x)
    cell_from_y <- intersect(cell_ids_retained, cell_ids_y)

    # Expand all matrices to retained dimensions
    ref_x_exp <- .expand_subset_matrix(x@ref_count, snp_ids_retained, cell_ids_retained, col_mapping_x)
    alt_x_exp <- .expand_subset_matrix(x@alt_count, snp_ids_retained, cell_ids_retained, col_mapping_x)
    oth_x_exp <- .expand_subset_matrix(x@oth_count, snp_ids_retained, cell_ids_retained, col_mapping_x)

    ref_y_exp <- .expand_subset_matrix(y@ref_count, snp_ids_retained, cell_ids_retained, col_mapping_y)
    alt_y_exp <- .expand_subset_matrix(y@alt_count, snp_ids_retained, cell_ids_retained, col_mapping_y)
    oth_y_exp <- .expand_subset_matrix(y@oth_count, snp_ids_retained, cell_ids_retained, col_mapping_y)

    # Merge by addition
    ref_merged <- ref_x_exp + ref_y_exp
    alt_merged <- alt_x_exp + alt_y_exp
    oth_merged <- oth_x_exp + oth_y_exp

    # Merge metadata using helper functions
    snp_info_merged <- .merge_snp_info(x, y, snp_ids_retained, snp_join)
    barcode_info_merged <- .merge_barcode_info(x, y, barcodes_retained, cell_join)

    # Create new SNPData object
    # The initialize method will automatically compute:
    # - coverage, non_zero_samples (in snp_info)
    # - library_size, non_zero_snps (in barcode_info)
    merged_obj <- SNPData(
        ref_count = ref_merged,
        alt_count = alt_merged,
        oth_count = oth_merged,
        snp_info = snp_info_merged,
        barcode_info = barcode_info_merged
    )

    logger::log_success(
        "Merged SNPData: {nrow(merged_obj)} SNPs x {ncol(merged_obj)} cells"
    )

    return(merged_obj)
}

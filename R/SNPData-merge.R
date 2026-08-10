#' Merge two SNPData objects
#'
#' Combines two SNPData objects with flexible join strategies for SNPs and cells.
#' Overlapping SNPs are summed, overlapping cells are resolved by
#' \code{cell_collision}, and unique entries are retained based on join type.
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
#'     \item{\code{"left"}}{Keep the cells of x, drawing on y where it covers
#'       the same cell: under the default \code{cell_collision}, a cell of x
#'       that y also carries can be replaced by y's copy rather than added to}
#'     \item{\code{"right"}}{Keep the cells of y, drawing on x on the same terms}
#'   }
#' @param cell_collision How to resolve a cell retained from both objects.
#'   Applied independently of \code{cell_join}: depth alone decides, so a
#'   \code{"left"} or \code{"right"} \code{cell_join} selects which cells appear
#'   but grants neither object priority over the other's counts. One of:
#'   \describe{
#'     \item{\code{"keep_best"}}{Keep whichever object's copy has the greater
#'       library size and discard the other, with an exact tie going to \code{x}
#'       (default). Correct when a shared identifier is a coincidence rather than
#'       a shared cell.}
#'     \item{\code{"sum"}}{Add the two copies' counts together, the unconditional
#'       behaviour of earlier versions. Correct only for technical replicates,
#'       the same physical cells sequenced twice.}
#'   }
#'
#' @return A new SNPData object containing the merged data
#'
#' @details
#' This function provides independent control over which SNPs (rows) and
#' cells (columns) to retain: SNPs present in both objects have their counts
#' summed element-wise, cells present in both are resolved according to
#' \code{cell_collision}, and positions present in only one object (but
#' retained by the join strategy) are zero-filled for the missing dimension.
#'
#' Metadata is merged using the corresponding dplyr join function: for
#' overlapping SNPs with conflicting metadata, values from x take priority.
#' A collided cell's metadata instead follows whichever copy
#' \code{cell_collision} retained, so its \code{donor} and other annotation
#' describe the cell whose counts survived (falling back to x priority under
#' \code{cell_collision = "sum"}). Auto-computed columns (\code{coverage},
#' \code{non_zero_samples}, \code{library_size}, \code{non_zero_snps}) are
#' recalculated for the merged object.
#'
#' \code{donor_info} is merged the same way (x priority, \code{n_cells}
#' recalculated). \code{donor_snp_info} instead \strong{errors} if x and y
#' disagree on \code{zygosity}, \code{zygosity_source}, \code{allele_on_x1},
#' or \code{xci_informative} for the same (SNP, donor) pair, since silently
#' picking one would hide a real reproducibility problem; numeric estimate
#' columns (p-values, GT posteriors, escape fractions) still take x priority
#' on conflict, since two independently fit or tested values can differ by
#' floating point alone. A donor with no cells left after \code{cell_join}
#' has its \code{donor_info}/\code{donor_snp_info} rows dropped, as when
#' subsetting.
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
#' @section Cell merging strategy:
#' Cells are matched and merged by barcode, not \code{cell_id}, when both
#' SNPData objects carry a \code{barcode} column in \code{barcode_info}; this
#' correctly merges the same physical cell across datasets even when its
#' internal \code{cell_id} differs between them. If \code{barcode} is absent
#' from either object, the function falls back to merging by \code{cell_id}.
#'
#' @section Cells that collide rather than correspond:
#' A shared barcode across independently barcoded libraries is often pure
#' chance, not a shared cell (10x's ~737,000-barcode whitelist gives two
#' 5,000-cell runs ~34 chance collisions). Summing such a pair fabricates a
#' cell carrying both donors' alleles, indistinguishable for X-linked genes
#' from XCI escape. The default \code{cell_collision = "keep_best"} avoids
#' this by discarding one member of the pair rather than summing;
#' \code{"sum"} is correct only for true technical replicates. To retain all
#' cells from independent libraries safely, suffix each object's barcodes
#' with a batch label before merging so they cannot collide.
#'
#' @examples
#' \dontrun{
#' # Technical replicates - sum depth for validated entries only
#' rep1 <- import_cellsnp("replicate1/")
#' rep2 <- import_cellsnp("replicate2/")
#' combined <- merge_snpdata(rep1, rep2,
#'                          snp_join = "intersect",
#'                          cell_join = "intersect",
#'                          cell_collision = "sum")
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
    cell_join = c("union", "intersect", "left", "right"),
    cell_collision = c("keep_best", "sum")
) {
    # Match arguments
    snp_join <- match.arg(snp_join)
    cell_join <- match.arg(cell_join)
    cell_collision <- match.arg(cell_collision)

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

    # Retained-position lookup for each object's rows/columns (NA when not
    # retained), computed once and shared across the ref/alt/oth channels
    # since all three carry identical dimnames.
    row_target_x <- .retained_position(snp_ids_x, snp_ids_retained)
    row_target_y <- .retained_position(snp_ids_y, snp_ids_retained)
    col_target_x <- col_mapping_x %||% .retained_position(cell_ids_x, cell_ids_retained)
    col_target_y <- col_mapping_y %||% .retained_position(cell_ids_y, cell_ids_retained)

    # A cell retained from both objects is summed only when the caller says the
    # two columns are the same physical cell; otherwise the deeper copy is kept
    # and the other discarded, so a barcode shared by chance between two
    # independently barcoded libraries cannot fuse two cells into one.
    winner_is_y <- NULL
    if (cell_collision == "keep_best") {
        resolution <- .resolve_cell_collisions(x, y, col_target_x, col_target_y, length(cell_ids_retained))
        col_target_x <- resolution$col_target_x
        col_target_y <- resolution$col_target_y
        winner_is_y <- resolution$winner_is_y
    }

    # Merge each allele channel directly from x's and y's triplets in one
    # pass, rather than expanding x and y to the merged dimensions separately
    # and adding -- sparseMatrix() sums duplicate (i, j) entries automatically.
    ref_merged <- .merge_count_matrix(
        x@ref_count,
        y@ref_count,
        row_target_x,
        col_target_x,
        row_target_y,
        col_target_y,
        snp_ids_retained,
        cell_ids_retained
    )
    alt_merged <- .merge_count_matrix(
        x@alt_count,
        y@alt_count,
        row_target_x,
        col_target_x,
        row_target_y,
        col_target_y,
        snp_ids_retained,
        cell_ids_retained
    )
    oth_merged <- .merge_count_matrix(
        x@oth_count,
        y@oth_count,
        row_target_x,
        col_target_x,
        row_target_y,
        col_target_y,
        snp_ids_retained,
        cell_ids_retained
    )

    # Merge metadata using helper functions
    snp_info_merged <- .merge_snp_info(x, y, snp_ids_retained, snp_join)
    barcode_info_merged <- .merge_barcode_info(x, y, barcodes_retained, cell_join, winner_is_y)

    # Donor genotypes are a property of the donor, not of the retained cells:
    # a donor with no surviving cells after the cell_join has its rows dropped
    # from both donor tables, same as `[` subsetting.
    donors_retained <- if ("donor" %in% colnames(barcode_info_merged)) {
        unique(stats::na.omit(barcode_info_merged$donor))
    } else {
        character(0)
    }
    donor_info_merged <- .merge_donor_info(x, y, donors_retained)
    donor_snp_info_merged <- .merge_donor_snp_info(x, y, snp_ids_retained, donors_retained)
    zygosity_source_merged <- .merge_zygosity_source(x, y)

    # Create new SNPData object
    # The initialize method will automatically compute:
    # - coverage, non_zero_samples (in snp_info)
    # - library_size, non_zero_snps (in barcode_info)
    # - n_cells (in donor_info)
    # zygosity_source is set explicitly below rather than left to the
    # constructor's auto-derivation, which collapses to NA once the merged
    # donor_snp_info carries more than one source -- see .merge_zygosity_source().
    merged_obj <- SNPData(
        ref_count = ref_merged,
        alt_count = alt_merged,
        oth_count = oth_merged,
        snp_info = snp_info_merged,
        barcode_info = barcode_info_merged,
        donor_info = donor_info_merged,
        donor_snp_info = donor_snp_info_merged
    )
    merged_obj@zygosity_source <- zygosity_source_merged

    logger::log_success(
        "Merged SNPData: {nrow(merged_obj)} SNPs x {ncol(merged_obj)} cells"
    )

    return(merged_obj)
}

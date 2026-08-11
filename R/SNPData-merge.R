#' Merge two SNPData objects
#'
#' Combines two SNPData objects with flexible join strategies for SNPs and cells.
#' Overlapping SNPs are summed, cells are matched on
#' \code{(library_id, barcode)}, and unique entries are retained based on join
#' type.
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
#' @param cell_join Join type for cells (columns). Two cells are "the same
#'   cell" only when they share both a \code{library_id} and a
#'   \code{barcode}; a cell in only one object is never merged into another,
#'   regardless of this setting. One of:
#'   \describe{
#'     \item{\code{"union"}}{Keep all cells from both objects (default)}
#'     \item{\code{"intersect"}}{Keep only cells present in both objects}
#'     \item{\code{"left"}}{Keep the cells of x, adding y's counts for
#'       whichever of them y also carries}
#'     \item{\code{"right"}}{Keep the cells of y, adding x's counts for
#'       whichever of them x also carries}
#'   }
#'
#' @return A new SNPData object containing the merged data
#'
#' @details
#' This function provides independent control over which SNPs (rows) and
#' cells (columns) to retain. A cell is identified by \code{(library_id,
#' barcode)}, not by \code{barcode} alone, because a 10x barcode is unique
#' only within the library it was drawn from (see \sQuote{Cells that collide
#' rather than correspond} below); every cell in both \code{x} and \code{y}
#' must therefore carry a non-\code{NA} \code{library_id}, or the function
#' errors. SNPs present in both objects have their counts summed
#' element-wise; cells present in both objects (identical \code{library_id}
#' and \code{barcode}) likewise have their counts summed, since a shared key
#' means the same physical cell. Positions present in only one object, but
#' retained by the join strategy, are zero-filled for the missing dimension.
#'
#' Metadata is merged using the corresponding dplyr join function: for
#' overlapping SNPs or cells with conflicting metadata, values from x take
#' priority. Auto-computed columns (\code{coverage}, \code{non_zero_samples},
#' \code{library_size}, \code{non_zero_snps}) are recalculated for the merged
#' object.
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
#' Cells are matched and merged by \code{(library_id, barcode)}, not
#' \code{cell_id}, when both SNPData objects carry a \code{barcode} column in
#' \code{barcode_info}; this correctly merges the same physical cell across
#' datasets even when its internal \code{cell_id} differs between them. If
#' \code{barcode} is absent from either object, the function falls back to
#' matching on \code{(library_id, cell_id)}. A cell's \code{cell_id} is kept
#' from the input where possible, and only renumbered (\code{"cell_1"},
#' \code{"cell_2"}, ...) when the retained identifiers would otherwise
#' collide, e.g. when two libraries happen to reuse the same \code{cell_id}.
#'
#' @section Cells that collide rather than correspond:
#' A shared barcode across independently barcoded libraries is often pure
#' chance, not a shared cell (10x's ~737,000-barcode whitelist gives two
#' 5,000-cell runs ~34 chance collisions). Summing such a pair fabricates a
#' cell carrying both donors' alleles, indistinguishable for X-linked genes
#' from XCI escape. \code{library_id} is what prevents this: cells are
#' matched on \code{(library_id, barcode)}, so a barcode shared across two
#' different library labels is kept as two separate cells rather than
#' summed. \code{merge_snpdata()} therefore requires every cell in both
#' objects to carry a non-\code{NA} \code{library_id}; set it at import with
#' \code{import_cellsnp(..., library_id = )} or afterwards via
#' \code{barcode_info(obj)$library_id <- }.
#'
#' @examples
#' \dontrun{
#' # Two cellSNP runs over the SAME library (e.g. rerun against a wider SNP
#' # panel): the shared library_id means a repeated barcode is the same cell,
#' # so depth is pooled.
#' rep1 <- import_cellsnp("replicate1/", gene_anno, library_id = "lib1")
#' rep2 <- import_cellsnp("replicate2/", gene_anno, library_id = "lib1")
#' combined <- merge_snpdata(rep1, rep2,
#'                          snp_join = "intersect",
#'                          cell_join = "intersect")
#'
#' # Batch integration - two distinct libraries, so no cell is shared and
#' # chance barcode collisions stay as separate cells
#' batch1 <- import_cellsnp("donor1/", gene_anno, library_id = "batch1")
#' batch2 <- import_cellsnp("donor2/", gene_anno, library_id = "batch2")
#' integrated <- merge_snpdata(batch1, batch2,
#'                             snp_join = "intersect",
#'                             cell_join = "union")
#'
#' # SNP panel expansion over one library - all SNPs, validated cells only
#' common <- import_cellsnp("common_vars/", gene_anno, library_id = "lib1")
#' rare <- import_cellsnp("rare_vars/", gene_anno, library_id = "lib1")
#' expanded <- merge_snpdata(common, rare,
#'                          snp_join = "union",
#'                          cell_join = "intersect")
#'
#' # Default: maximum data retention
#' dataset1 <- import_cellsnp("cohort_A/", gene_anno, library_id = "cohort_A")
#' dataset2 <- import_cellsnp("cohort_B/", gene_anno, library_id = "cohort_B")
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

    # A cell is identified by (library_id, barcode), not by barcode alone: a 10x
    # barcode is unique only within the library it was drawn from, so the pair is
    # what distinguishes the same cell sequenced twice from two different cells
    # that happened to draw the same barcode.
    .check_library_ids(x, y)

    # Barcodes are the real cell identity; cell_id is a positional label the
    # constructor generates, so it is only a fallback for objects built without
    # barcodes (see the class docs).
    id_col <- if ("barcode" %in% colnames(x@barcode_info) && "barcode" %in% colnames(y@barcode_info)) {
        "barcode"
    } else {
        "cell_id"
    }

    keys_x <- .cell_keys(x@barcode_info, id_col)
    keys_y <- .cell_keys(y@barcode_info, id_col)

    keys_retained <- switch(
        cell_join,
        "union" = union(keys_x, keys_y),
        "intersect" = intersect(keys_x, keys_y),
        "left" = keys_x,
        "right" = keys_y
    )

    col_mapping_x <- match(keys_x, keys_retained)
    col_mapping_y <- match(keys_y, keys_retained)

    # Determine which SNPs to retain based on snp_join strategy
    snp_ids_retained <- switch(
        snp_join,
        "union" = union(snp_ids_x, snp_ids_y),
        "intersect" = intersect(snp_ids_x, snp_ids_y),
        "left" = snp_ids_x,
        "right" = snp_ids_y
    )

    # Handle edge case: empty result. Tested on the retained keys rather than on
    # the cell_ids derived from them below, since paste0() over a zero-length
    # index still yields one element and would mask an empty cell set.
    if (length(snp_ids_retained) == 0) {
        stop("No SNPs to retain after merge. Check your join strategy and input data.")
    }
    if (length(keys_retained) == 0) {
        stop("No cells to retain after merge. Check your join strategy and input data.")
    }

    # Column names come from the retained keys' own identifiers where those stay
    # unique, so a merge that does not actually mix libraries keeps its cell_ids
    # recognisable; two libraries sharing an identifier force a renumber, since
    # column names and cell_id must stay unique.
    ids_retained <- .key_identifiers(
        keys_retained,
        c(keys_x, keys_y),
        c(x@barcode_info[[id_col]], y@barcode_info[[id_col]])
    )
    cell_ids_retained <- if (id_col == "cell_id" && !anyDuplicated(ids_retained)) {
        ids_retained
    } else {
        paste0("cell_", seq_along(keys_retained))
    }

    # Retained-position lookup for each object's rows (NA when not retained),
    # computed once and shared across the ref/alt/oth channels since all three
    # carry identical dimnames. The column equivalents are the key mappings
    # built above, so a cell present in both objects lands on one column and
    # .merge_count_matrix() sums its two contributions.
    row_target_x <- .retained_position(snp_ids_x, snp_ids_retained)
    row_target_y <- .retained_position(snp_ids_y, snp_ids_retained)
    col_target_x <- col_mapping_x
    col_target_y <- col_mapping_y

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
    barcode_info_merged <- .merge_barcode_info(x, y, keys_retained, cell_join, id_col, cell_ids_retained)

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
    merged_obj <- .merge_library_bams(merged_obj, x, y)

    logger::log_success(
        "Merged SNPData: {nrow(merged_obj)} SNPs x {ncol(merged_obj)} cells"
    )

    return(merged_obj)
}

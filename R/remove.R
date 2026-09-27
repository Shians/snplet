#' Remove doublet cells from a SNPData object
#'
#' Removes cells whose \code{donor} in \code{barcode_info} is \code{"doublet"}:
#' barcodes Vireo judged to hold two cells from different donors. Their reads
#' mix two genotypes, so they are unreliable at every level, including
#' per-barcode.
#'
#' @param x A SNPData object, required.
#' @param drop_na Logical (default \code{TRUE}). Whether to also remove
#'   cells with NA donor assignments.
#'
#' @return A filtered SNPData object with doublets removed
#' @seealso \code{\link{remove_unassigned}} for cells Vireo could not assign,
#'   and \code{\link{keep_singlets}} to remove both.
#' @family data cleaning functions
#' @export
#'
#' @examples
#' \dontrun{
#' snp_data <- get_example_snpdata()
#' # Remove doublets from SNPData object
#' filtered_data <- remove_doublets(snp_data)
#' }
remove_doublets <- function(x, drop_na = TRUE) {
    .remove_donor_labels(x, "doublet", drop_na)
}

#' Remove unassigned cells from a SNPData object
#'
#' Removes cells whose \code{donor} in \code{barcode_info} is
#' \code{"unassigned"}: single cells Vireo could not confidently assign to a
#' donor. Unlike doublets their per-barcode counts are valid, so removing them
#' is optional; donor-level functions exclude them regardless.
#'
#' @param x A SNPData object, required.
#' @param drop_na Logical (default \code{TRUE}). Whether to also remove
#'   cells with NA donor assignments.
#'
#' @return A filtered SNPData object with unassigned cells removed
#' @seealso \code{\link{remove_doublets}} for doublets, and
#'   \code{\link{keep_singlets}} to remove both.
#' @family data cleaning functions
#' @export
#'
#' @examples
#' \dontrun{
#' snp_data <- get_example_snpdata()
#' # Remove unassigned cells from SNPData object
#' filtered_data <- remove_unassigned(snp_data)
#' }
remove_unassigned <- function(x, drop_na = TRUE) {
    .remove_donor_labels(x, "unassigned", drop_na)
}

#' Keep only singlet cells in a SNPData object
#'
#' Keeps cells that Vireo assigned to exactly one donor, removing those whose
#' \code{donor} in \code{barcode_info} is \code{"doublet"}, \code{"unassigned"},
#' or \code{NA}. These are the same cells that donor-level functions such as
#' \code{\link{donor_count_df}} and \code{\link{assign_xci}} use.
#'
#' @param x A SNPData object, required. Its \code{barcode_info} must have a
#'   \code{donor} column.
#'
#' @return A filtered SNPData object holding only singlet cells.
#' @seealso \code{\link{remove_doublets}} and \code{\link{remove_unassigned}}
#'   to remove one label only.
#' @family data cleaning functions
#' @export
#'
#' @examples
#' \dontrun{
#' snp_data <- get_example_snpdata()
#' # Keep only cells assigned to a single donor
#' singlets <- keep_singlets(snp_data)
#' }
keep_singlets <- function(x) {
    if (!methods::is(x, "SNPData")) {
        stop("Input must be a SNPData object")
    }
    barcode_info <- barcode_info(x)
    if (!"donor" %in% colnames(barcode_info)) {
        stop(
            "Donor information not available. Add donor data using add_barcode_metadata() or import_cellsnp() with vireo_folder parameter."
        )
    }

    .drop_barcodes(x, !.is_real_donor(barcode_info$donor), "doublet, unassigned or NA-donor")
}

# Removes cells whose donor label is one of `labels` (and, if drop_na, cells
# with no donor).
.remove_donor_labels <- function(x, labels, drop_na) {
    if (!methods::is(x, "SNPData")) {
        stop("Input must be a SNPData object")
    }

    barcode_info <- barcode_info(x)
    if (!"donor" %in% colnames(barcode_info)) {
        warning("No 'donor' column found in barcode_info, returning original object")
        return(x)
    }

    cells_to_remove <- barcode_info$donor %in% labels | (drop_na & is.na(barcode_info$donor))
    .drop_barcodes(x, cells_to_remove, paste(labels, collapse = "/"))
}

# Subsets out the cells flagged in `cells_to_remove`, logging how many were
# removed under the description `what`.
.drop_barcodes <- function(x, cells_to_remove, what) {
    barcodes_total <- ncol(x)
    barcodes_removed <- sum(cells_to_remove)
    barcodes_remaining <- barcodes_total - barcodes_removed
    removed_perc <- scales::percent(barcodes_removed / barcodes_total, accuracy = 0.01)
    logger::log_info(
        "{barcodes_removed} ({removed_perc}) {what} barcodes removed. {barcodes_remaining} barcodes remaining."
    )

    x[, !cells_to_remove]
}

#' Remove SNPs with NA gene names
#'
#' This function removes SNPs that have NA values in the gene_name column
#' of the snp_info slot. This is useful for focusing analysis on SNPs
#' that are associated with known genes.
#'
#' @param x A SNPData object, required.
#' @param gene_col Character scalar (default \code{"gene_name"}). Column
#'   name in \code{snp_info} containing gene names.
#'
#' @return A filtered SNPData object with NA gene SNPs removed
#' @family data cleaning functions
#' @export
#'
#' @examples
#' snp_data <- get_example_snpdata()
#' # Filter out SNPs with NA gene names
#' filtered_data <- remove_na_genes(snp_data)
remove_na_genes <- function(x, gene_col = "gene_name") {
    # Check that x is a SNPData object
    if (!methods::is(x, "SNPData")) {
        stop("Input must be a SNPData object")
    }

    # Get SNP info
    snp_info <- snp_info(x)

    # Check if gene column exists
    if (!gene_col %in% colnames(snp_info)) {
        warning(paste0("No '", gene_col, "' column found in snp_info, returning original object"))
        return(x)
    }

    # Identify SNPs with non-NA gene names
    snps_to_keep <- !is.na(snp_info[[gene_col]])

    snps_total <- nrow(snp_info)
    snps_removed <- sum(!snps_to_keep)
    snps_remaining <- snps_total - snps_removed
    removed_perc <- scales::percent(snps_removed / snps_total, accuracy = 0.01)
    logger::log_info(
        "{snps_removed} ({removed_perc}) SNPs with NA gene_name values removed. {snps_remaining} SNPs remaining."
    )

    # Return filtered data
    return(x[snps_to_keep, ])
}

#' Remove barcodes with NA clonotype values
#'
#' This function filters out barcodes that have NA values in the clonotype column
#' of the barcode_info slot. This is useful for analyses that require valid
#' clonotype assignments for all barcodes.
#'
#' @param x A SNPData object, required.
#' @param clonotype_col Character scalar (default \code{"clonotype"}).
#'   Column name in \code{barcode_info} containing clonotype values.
#'
#' @return A filtered SNPData object with NA clonotype barcodes removed
#' @family data cleaning functions
#' @export
#'
#' @examples
#' snp_data <- get_example_snpdata()
#' # Remove barcodes with NA clonotype values
#' filtered_data <- remove_na_clonotypes(snp_data)
remove_na_clonotypes <- function(x, clonotype_col = "clonotype") {
    # Check that x is a SNPData object
    if (!methods::is(x, "SNPData")) {
        stop("Input must be a SNPData object")
    }

    # Get sample info
    barcode_info <- barcode_info(x)

    # Check if clonotype column exists
    if (!clonotype_col %in% colnames(barcode_info)) {
        warning(paste0(
            "No '",
            clonotype_col,
            "' column found in barcode_info. ",
            "Add clonotype information using add_barcode_metadata() or import_cellsnp() with vdj_file. ",
            "Returning original object."
        ))
        return(x)
    }

    # Identify barcodes with NA clonotype values
    barcodes_to_remove <- is.na(barcode_info[[clonotype_col]])

    clonotypes_total <- nrow(barcode_info)
    clonotypes_removed <- sum(barcodes_to_remove)
    clonotypes_remaining <- clonotypes_total - clonotypes_removed
    removed_perc <- scales::percent(clonotypes_removed / clonotypes_total, accuracy = 0.01)
    logger::log_info(
        "{clonotypes_removed} ({removed_perc}) barcodes with NA clonotype values removed. {clonotypes_remaining} barcodes remaining."
    )

    # Return filtered data
    return(x[, !barcodes_to_remove])
}

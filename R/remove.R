#' Remove doublet cells from a SNPData object
#'
#' This function filters out cells that are identified as doublets from a SNPData object.
#' Doublets are identified based on the 'donor' column in the barcode_info slot, where
#' cells labeled as 'doublet' are removed.
#'
#' @param x A SNPData object
#' @param drop_na Logical, whether to also remove cells with NA donor assignments (default TRUE)
#'
#' @return A filtered SNPData object with doublets removed
#' @export
#'
#' @examples
#' \dontrun{
#' snp_data <- get_example_snpdata()
#' # Remove doublets from SNPData object
#' filtered_data <- remove_doublets(snp_data)
#' }
remove_doublets <- function(x, drop_na = TRUE) {
    # Check that x is a SNPData object
    if (!methods::is(x, "SNPData")) {
        stop("Input must be a SNPData object")
    }

    # Get sample info
    barcode_info <- get_barcode_info(x)

    # Check if donor column exists
    if (!"donor" %in% colnames(barcode_info)) {
        warning("No 'donor' column found in barcode_info, returning original object")
        return(x)
    }

    # Identify doublets
    cells_to_remove <- barcode_info$donor == "doublet"

    # Handle NA values
    if (drop_na) {
        cells_to_remove[is.na(cells_to_remove)] <- TRUE
    } else {
        cells_to_remove[is.na(cells_to_remove)] <- FALSE
    }

    barcodes_total <- ncol(x)
    barcodes_removed <- sum(cells_to_remove)
    barcodes_remaining <- barcodes_total - barcodes_removed
    removed_perc <- scales::percent(barcodes_removed / barcodes_total, accuracy = 0.01)
    logger::log_info(
        "{barcodes_removed} ({removed_perc}) doublet barcodes removed. {barcodes_remaining} barcodes remaining."
    )

    # Return filtered data
    return(x[, !cells_to_remove])
}

#' Remove SNPs with NA gene names
#'
#' This function removes SNPs that have NA values in the gene_name column
#' of the snp_info slot. This is useful for focusing analysis on SNPs
#' that are associated with known genes.
#'
#' @param x A SNPData object
#' @param gene_col The column name in snp_info containing gene names (default: "gene_name")
#'
#' @return A filtered SNPData object with NA gene SNPs removed
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
    snp_info <- get_snp_info(x)

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
#' @param x A SNPData object
#' @param clonotype_col The column name in barcode_info containing clonotype values (default: "clonotype")
#'
#' @return A filtered SNPData object with NA clonotype barcodes removed
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
    barcode_info <- get_barcode_info(x)

    # Check if clonotype column exists
    if (!clonotype_col %in% colnames(barcode_info)) {
        warning(paste0(
            "No '", clonotype_col, "' column found in barcode_info. ",
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

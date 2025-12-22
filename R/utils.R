# Suppress R CMD check notes about global variables used in data.frame operations
utils::globalVariables(c(
    ".snp_row",
    "adj_p_val",
    "alt",
    "barcode",
    "binom_p_value",
    "cell",
    "cell_id",
    "chrom",
    "clonotype",
    "donor",
    "donor_id",
    "gene_name",
    "library_size",
    "maf",
    "minor_count",
    "n",
    "n_het_donors",
    "pos",
    "ref",
    "seqnames",
    "signif_snps_clonotype",
    "snp_id",
    "start",
    "total_count"
))

percentile_summary <- function(x, percentiles = c(0.1, 0.25, 0.75, 0.9, 0.95, 0.99)) {
    out <- c(min = min(x))
    past_median <- FALSE
    for (p in percentiles) {
        if (p > 0.5 && !past_median) {
            out <- c(out, median(x))
            names(out)[length(out)] <- "median"
            past_median <- TRUE
        }
        out <- c(out, as.numeric(quantile(x, p)))
        names(out)[length(out)] <- paste0("p", as.character(round(p * 100)))
    }

    out <- c(out, max = max(x))
    names(out)[length(out)] <- "max"

    out
}

groupedRowMeans <- function(x, groups) {
    # calculate the mean of each group
    out <- matrix(NA, ncol = length(unique(groups)), nrow = nrow(x))
    colnames(out) <- sort(unique(groups))
    rownames(out) <- rownames(x)

    for (i in unique(groups)) {
        out[, i] <- Matrix::rowMeans(
            x[, groups == i, drop = FALSE],
            na.rm = TRUE
        )
    }

    out
}

groupedRowSums <- function(x, groups) {
    # calculate the sum of each group
    out <- matrix(NA, ncol = length(unique(groups)), nrow = nrow(x))
    colnames(out) <- sort(unique(groups))
    rownames(out) <- rownames(x)

    for (i in unique(groups)) {
        out[, i] <- Matrix::rowSums(
            x[, groups == i, drop = FALSE],
            na.rm = TRUE
        )
    }

    out
}

#' Generate standardized SNP IDs from genomic coordinates
#'
#' Creates SNP identifiers in the format "chr:pos:ref:alt" for consistent
#' identification across datasets. This format is deterministic and contains
#' all information needed for variant identification.
#'
#' @param chrom Character vector of chromosome names
#' @param pos Integer vector of genomic positions
#' @param ref Character vector of reference alleles
#' @param alt Character vector of alternate alleles
#'
#' @return Character vector of SNP IDs in format "chr:pos:ref:alt"
#'
#' @examples
#' \dontrun{
#' snp_ids <- make_snp_id(
#'     chrom = c("chr1", "chr2"),
#'     pos = c(12345, 67890),
#'     ref = c("A", "G"),
#'     alt = c("G", "T")
#' )
#' # Returns: c("chr1:12345:A:G", "chr2:67890:G:T")
#' }
#'
#' @keywords internal
make_snp_id <- function(chrom, pos, ref, alt) {
    paste(chrom, pos, ref, alt, sep = ":")
}

#' Check if a file exists
#' @param path Path to the file to check
check_file <- function(path) {
    if (file.exists(path)) {
        logger::log_info("File found: {path}")
    } else {
        stop(glue::glue("Required file not found: {path}"))
    }
}

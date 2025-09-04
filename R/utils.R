# Suppress R CMD check notes about global variables used in data.frame operations
utils::globalVariables(c(
    "adj_p_val",
    "alt",
    "barcode",
    "cell",
    "cell_id",
    "chrom",
    "clonotype",
    "donor",
    "donor_id",
    "gene_name",
    "library_size",
    "maf",
    "n",
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

#' Check if a file exists
#' @param path Path to the file to check
check_file <- function(path) {
    if (file.exists(path)) {
        logger::log_info("File found: {path}")
    } else {
        stop(glue::glue("Required file not found: {path}"))
    }
}

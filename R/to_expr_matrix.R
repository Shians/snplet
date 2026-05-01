#' Convert SNPData to expression-like matrix
#'
#' Transforms SNP allele counts into an expression-like matrix that captures both
#' allelic imbalance and sequencing depth. This matrix can be used for downstream
#' analyses such as dimensionality reduction, clustering, or differential expression.
#'
#' @details
#' The expression matrix is calculated using the following formula:
#'
#' \deqn{E = \frac{REF - ALT}{depth + 1} \times \log(depth + 1)}
#'
#' Where:
#' \itemize{
#'   \item \code{depth = REF + ALT} is the total read depth at each SNP position
#'   \item \code{(REF - ALT) / (depth + 1)} represents the allelic proportion, ranging
#'         from -1 (all ALT) to +1 (all REF), with 0 indicating balanced expression
#'   \item \code{log(depth + 1)} weights the signal by sequencing depth, similar to
#'         log-normalization in RNA-seq analysis
#'   \item The pseudocount of 1 prevents division by zero and log(0)
#' }
#'
#' This transformation has several useful properties:
#' \itemize{
#'   \item SNPs with higher depth contribute more to the signal
#'   \item Allelic imbalance direction is preserved (positive = REF bias, negative = ALT bias)
#'   \item Low-depth SNPs are naturally down-weighted
#'   \item The log transformation stabilizes variance across depth ranges
#' }
#'
#' @param x SNPData object
#' @param level Aggregation level: "barcode", "clonotype", or "donor"
#' @return Matrix with snp_id as rows and cell/clonotype/donor as columns
#' @export
setGeneric("to_expr_matrix", function(x, level = c("barcode", "clonotype", "donor")) standardGeneric("to_expr_matrix"))

#' @rdname to_expr_matrix
#' @include SNPData-class.R
setMethod("to_expr_matrix", signature(x = "SNPData"), function(x, level = c("barcode", "clonotype", "donor")) {
    level <- match.arg(level)
    barcode_info <- get_barcode_info(x)

    if (level != "barcode" && !level %in% colnames(barcode_info)) {
        if (level == "clonotype") {
            stop(
                "Clonotype information not available. Add clonotype data using add_barcode_metadata() or import_cellsnp() with vdj_file parameter."
            )
        } else {
            stop("No ", level, " column in barcode_info.")
        }
    }

    if (level == "barcode") {
        ref <- ref_count(x)
        alt <- alt_count(x)
        depth <- ref + alt
        prop <- (ref - alt) / (depth + 1)
        mat <- prop * log1p(depth)
        rownames(mat) <- rownames(ref)
        colnames(mat) <- colnames(ref)
        return(mat)
    } else if (level == "clonotype") {
        grouped <- .prepare_grouped_counts(x, barcode_info$clonotype, "clonotype")
        ref <- groupedRowSums(ref_count(grouped$x), grouped$groups)
        alt <- groupedRowSums(alt_count(grouped$x), grouped$groups)
        depth <- ref + alt
        prop <- (ref - alt) / (depth + 1)
        mat <- prop * log1p(depth)
        rownames(mat) <- rownames(ref)
        colnames(mat) <- colnames(ref)
        return(mat)
    } else if (level == "donor") {
        grouped <- .prepare_grouped_counts(x, barcode_info$donor, "donor")
        ref <- groupedRowSums(ref_count(grouped$x), grouped$groups)
        alt <- groupedRowSums(alt_count(grouped$x), grouped$groups)
        depth <- ref + alt
        prop <- (ref - alt) / (depth + 1)
        mat <- prop * log1p(depth)
        rownames(mat) <- rownames(ref)
        colnames(mat) <- colnames(ref)
        return(mat)
    }
})

#' @keywords internal
.prepare_grouped_counts <- function(x, groups, level) {
    if (all(is.na(groups))) {
        if (level == "clonotype") {
            stop(
                "All clonotype values are NA. Cannot perform clonotype-level aggregation. Add clonotype data using add_barcode_metadata() or import_cellsnp() with vdj_file parameter."
            )
        }
        stop(
            "All donor values are NA. Cannot perform donor-level aggregation. Add donor data using add_barcode_metadata() or import_cellsnp() with vireo_file parameter."
        )
    }

    keep <- !is.na(groups)
    if (any(!keep)) {
        x <- x[, keep, drop = FALSE]
        groups <- groups[keep]
    }

    list(x = x, groups = groups)
}

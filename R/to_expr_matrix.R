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
#'         log-normalisation in RNA-seq analysis
#'   \item The pseudocount of 1 prevents division by zero and log(0)
#' }
#'
#' This transformation has several useful properties:
#' \itemize{
#'   \item SNPs with higher depth contribute more to the signal
#'   \item Allelic imbalance direction is preserved (positive = REF bias, negative = ALT bias)
#'   \item Low-depth SNPs are naturally down-weighted
#'   \item The log transformation stabilises variance across depth ranges
#' }
#'
#' @param x A SNPData object, required.
#' @param level Character scalar, one of \code{"barcode"}, \code{"clonotype"},
#'   or \code{"donor"} (default \code{"barcode"}). Aggregation level.
#' @return Matrix with snp_id as rows and cell/clonotype/donor as columns
#' @export
#'
#' @examples
#' \dontrun{
#' snp_data <- get_example_snpdata()
#' expr_matrix <- to_expr_matrix(snp_data)
#' }
setGeneric("to_expr_matrix", function(x, level = c("barcode", "clonotype", "donor")) standardGeneric("to_expr_matrix"))

#' @rdname to_expr_matrix
#' @include SNPData-class.R
setMethod("to_expr_matrix", signature(x = "SNPData"), function(x, level = c("barcode", "clonotype", "donor")) {
    level <- match.arg(level)
    barcode_info <- barcode_info(x)

    if (level != "barcode" && !level %in% colnames(barcode_info)) {
        .stop_level_unavailable(level)
    }

    if (level == "barcode") {
        mat <- .compute_expr_matrix(ref_count(x), alt_count(x))
    } else {
        grouped <- .prepare_grouped_counts(x, barcode_info[[level]], level)
        ref <- groupedRowSums(ref_count(grouped$x), grouped$groups)
        alt <- groupedRowSums(alt_count(grouped$x), grouped$groups)
        mat <- .compute_expr_matrix(ref, alt)
    }

    mat
})

# Implements the allelic-imbalance formula documented in @details: an
# allelic proportion in [-1, 1] weighted by log-depth, with a pseudocount of
# 1 guarding the division and the log against zero depth.
.compute_expr_matrix <- function(ref, alt) {
    depth <- ref + alt
    prop <- (ref - alt) / (depth + 1)
    mat <- prop * log1p(depth)
    rownames(mat) <- rownames(ref)
    colnames(mat) <- colnames(ref)
    mat
}

.stop_level_unavailable <- function(level) {
    if (level == "clonotype") {
        stop(
            "Clonotype information not available. Add clonotype data using add_barcode_metadata() or import_cellsnp() with vdj_file parameter."
        )
    } else {
        stop(
            "Donor information not available. Add donor data using add_barcode_metadata() or import_cellsnp() with vireo_folder parameter."
        )
    }
}

#' @keywords internal
.prepare_grouped_counts <- function(x, groups, level) {
    if (all(is.na(groups))) {
        if (level == "clonotype") {
            stop(
                "All clonotype values are NA. Cannot perform clonotype-level expression matrix conversion. Add clonotype data using add_barcode_metadata() or import_cellsnp() with vdj_file parameter."
            )
        }
        stop(
            "All donor values are NA. Cannot perform donor-level expression matrix conversion. Add donor data using add_barcode_metadata() or import_cellsnp() with vireo_folder parameter."
        )
    }

    keep <- !is.na(groups)
    if (any(!keep)) {
        x <- x[, keep, drop = FALSE]
        groups <- groups[keep]
    }

    list(x = x, groups = groups)
}

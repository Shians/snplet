#' Test for minor allele frequency significance
#'
#' Performs a one-sided exact binomial test for each row to assess whether the
#' observed minor allele count is significantly greater than expected under a null
#' minor allele frequency (MAF). Assumes each row's counts are independent; the
#' returned p-values are adjusted for multiple testing across all rows.
#'
#' @param x A data.frame (or tibble) with one row per observation, containing the
#'   integer columns \code{ref_count}, \code{alt_count} and \code{total_count}.
#'   All counts must be non-negative and satisfy
#'   \code{ref_count + alt_count <= total_count}.
#' @param p Numeric scalar giving the null hypothesis minor allele frequency, in
#'   \code{[0, 1]} (default: 0.10).
#'
#' @return The input \code{x} with three columns appended:
#'   \describe{
#'     \item{\code{minor_allele_count}}{\code{pmin(ref_count, alt_count)} for each row.}
#'     \item{\code{p_val}}{One-sided exact binomial p-value for observing at least
#'       \code{minor_allele_count} minor alleles out of \code{ceiling(total_count)}
#'       trials under \code{Binomial(n, p)}.}
#'     \item{\code{adj_p_val}}{Benjamini-Hochberg (BH) adjusted p-value.}
#'   }
#'
#' @details
#' The test is an exact binomial test, not a beta-binomial test: it does not model
#' overdispersion in the allele counts. The one-sided p-value \eqn{P(X \ge x)} for
#' \eqn{X \sim \mathrm{Binomial}(n, p)} is evaluated exactly via the regularised
#' incomplete beta function rather than by simulation. The number of trials is taken
#' as \code{ceiling(total_count)} so that fractional counts (e.g. from aggregation)
#' yield a valid integer size, and the major allele count is floored at zero.
#'
#' @keywords internal
#' @examples
#' df <- tibble::tibble(
#'     ref_count = c(10, 5),
#'     alt_count = c(2, 8),
#'     total_count = c(12, 13)
#' )
#' test_maf(df)
test_maf <- function(x, p = 0.10) {
    # validate input
    stopifnot(is(x, "data.frame"))

    req_cols <- c("ref_count", "alt_count", "total_count")
    missing_cols <- setdiff(req_cols, colnames(x))
    if (length(missing_cols) > 0) {
        missing_list <- paste0(missing_cols, collapse = ", ")
        stop(glue::glue("Missing required columns: {missing_list}"))
    }

    # validate that all counts are non-negative
    negative_ref <- which(x$ref_count < 0)
    if (length(negative_ref) > 0) {
        if (length(negative_ref) <= 5) {
            row_list <- paste0(negative_ref, collapse = ", ")
        } else {
            row_list <- paste0(paste0(head(negative_ref, 5), collapse = ", "), ", ...")
        }
        stop(glue::glue("Invalid data: ref_count < 0 in row(s): {row_list}"))
    }

    negative_alt <- which(x$alt_count < 0)
    if (length(negative_alt) > 0) {
        if (length(negative_alt) <= 5) {
            row_list <- paste0(negative_alt, collapse = ", ")
        } else {
            row_list <- paste0(paste0(head(negative_alt, 5), collapse = ", "), ", ...")
        }
        stop(glue::glue("Invalid data: alt_count < 0 in row(s): {row_list}"))
    }

    negative_total <- which(x$total_count < 0)
    if (length(negative_total) > 0) {
        if (length(negative_total) <= 5) {
            row_list <- paste0(negative_total, collapse = ", ")
        } else {
            row_list <- paste0(paste0(head(negative_total, 5), collapse = ", "), ", ...")
        }
        stop(glue::glue("Invalid data: total_count < 0 in row(s): {row_list}"))
    }

    # validate that ref_count + alt_count <= total_count
    invalid_rows <- which(x$ref_count + x$alt_count > x$total_count)
    if (length(invalid_rows) > 0) {
        if (length(invalid_rows) <= 5) {
            row_list <- paste0(invalid_rows, collapse = ", ")
        } else {
            row_list <- paste0(paste0(head(invalid_rows, 5), collapse = ", "), ", ...")
        }
        stop(glue::glue("Invalid data: ref_count + alt_count > total_count in row(s): {row_list}"))
    }

    minor_allele_count <- pmin(x$ref_count, x$alt_count)
    total_count <- ceiling(x$total_count)
    major_allele_count <- pmax(total_count - minor_allele_count, 0)

    p_val <- binom_test(minor_allele_count, total_count, p, alternative = "greater")

    result <- x %>%
        dplyr::mutate(
            minor_allele_count = minor_allele_count,
            p_val = p_val,
            adj_p_val = p.adjust(p_val, method = "BH")
        )

    return(result)
}

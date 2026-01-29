#' Test for minor allele frequency significance
#'
#' Performs a binomial test for each row to assess whether the observed minor allele count
#' is significantly greater than expected under a null minor allele frequency (MAF).
#'
#' @param x A data.frame containing columns 'ref_count', 'alt_count', and 'total_count'
#' @param p Null hypothesis minor allele frequency (default: 0.10)
#' @return A data.frame with columns for minor_allele_count, p_val, and adj_p_val (BH adjusted)
#' @examples
#' \dontrun{
#' df <- tibble::tibble(ref_count = c(10, 5), alt_count = c(2, 8), total_count = c(12, 13))
#' test_maf(df)
#' }
test_maf <- function(x, p = 0.10) {
    # validate input
    stopifnot(is(x, "data.frame"))

    req_cols <- c("ref_count", "alt_count", "total_count")
    missing_cols <- setdiff(req_cols, colnames(x))
    if (length(missing_cols) > 0) {
        missing_list <- paste0(missing_cols, collapse = ", ")
        stop(glue::glue("Missing required columns: {missing_list}"))
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

    # equivalent to binom.test(minor_allele_count, total_count, p, alternative = "greater")$p.value
    # pbeta is vectorised and much faster for this test
    p_val <- stats::pbeta(p, minor_allele_count, major_allele_count + 1, lower.tail = TRUE)

    result <- x %>%
        dplyr::mutate(
            minor_allele_count = minor_allele_count,
            p_val = p_val,
            adj_p_val = p.adjust(p_val, method = "BH")
        )

    return(result)
}

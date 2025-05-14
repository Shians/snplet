test_maf <- function(x, p = 0.10) {
    # validate input
    stopifnot(is(x, "data.frame"))

    req_cols <- c("ref_count", "alt_count", "total_count")
    missing_cols <- setdiff(req_cols, colnames(x))
    if (length(missing_cols) > 0) {
        stop(glue::glue("Missing required columns: {paste(missing_cols, collapse = ', ')}"))
    }

    minor_allele_count <- pmin(x$ref_count, x$alt_count)

    p_val <- furrr::future_map2_dbl(
        minor_allele_count,
        ceiling(x$total_count),
        ~binom.test(.x, .y, p, alternative = "greater")$p.value
    )

    result <- x %>%
        dplyr::mutate(
            minor_allele_count = minor_allele_count,
            p_val = p_val,
            adj_p_val = p.adjust(p_val, method = "BH")
        )

    return(result)
}

test_maf <- function(x, p = 0.10) {
    # validate input
    stopifnot(is(x, "data.frame"))
    stopifnot(all(c("ref_count", "alt_count", "total_count") %in% colnames(x)))

    minor_allele_count <- pmin(x$ref_count, x$alt_count)

    p_val <- purrr::map2_dbl(
        minor_allele_count,
        ceiling(x$total_count/2),
        ~binom.test(.x, .y, p, alternative = "greater")$p.value
    )

    x %>%
        dplyr::mutate(
            minor_allele_count = minor_allele_count,
            p_val = p_val,
            adj_p_val = p.adjust(p_val, method = "BH")
        )
}

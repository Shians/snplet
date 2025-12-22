#' Get donor-level SNP heterozygosity status
#'
#' Returns a long-format data frame with heterozygosity status per SNP and donor.
#'
#' @param x A SNPData object
#' @param min_total_count Minimum total read depth (ref + alt) required per donor to test for heterozygosity (default: 10)
#' @param p_value_threshold P-value threshold for binomial test (default: 0.05). P-values are multiple-testing corrected and SNPs with p < threshold reject monoallelic expression.
#' @param minor_allele_prop Minor allele proportion used as the null threshold for monoallelic expression testing (default: 0.1).
#' @return A tibble with columns: snp_id, donor, ref_count, alt_count, total_count, minor_allele_count, p_val, adj_p_val, tested, status
#' @export
#'
#' @examples
#' \dontrun{
#' snp_data <- get_example_snpdata()
#' donor_het_status_df(snp_data)
#' }
#' @rdname donor_het_status_df
setGeneric("donor_het_status_df", function(
    x,
    min_total_count = 10,
    p_value_threshold = 0.05,
    minor_allele_prop = 0.1
) standardGeneric("donor_het_status_df"))

donor_het_status_df_impl <- function(
    x,
    min_total_count = 10,
    p_value_threshold = 0.05,
    minor_allele_prop = 0.1
) {
    stopifnot(min_total_count >= 1)
    stopifnot(p_value_threshold > 0 && p_value_threshold <= 1)
    stopifnot(minor_allele_prop > 0 && minor_allele_prop < 0.5)

    donor_counts <- donor_count_df(x, test_maf = FALSE) %>%
        dplyr::mutate(tested = total_count >= min_total_count)

    tested_counts <- donor_counts %>%
        dplyr::filter(tested)

    if (nrow(tested_counts) > 0) {
        tested_counts <- tested_counts %>%
            test_maf(p = minor_allele_prop) %>%
            dplyr::select(
                snp_id,
                donor,
                minor_allele_count,
                p_val,
                adj_p_val
            )
    } else {
        tested_counts <- donor_counts %>%
            dplyr::select(snp_id, donor) %>%
            dplyr::mutate(
                minor_allele_count = NA_real_,
                p_val = NA_real_,
                adj_p_val = NA_real_
            )
    }

    donor_counts %>%
        dplyr::left_join(tested_counts, by = c("snp_id", "donor")) %>%
        dplyr::mutate(
            status = dplyr::if_else(
                tested & !is.na(adj_p_val) & adj_p_val < p_value_threshold,
                "het",
                "hom"
            )
        )
}

#' @rdname donor_het_status_df
setMethod(
    "donor_het_status_df",
    signature(x = "SNPData"),
    donor_het_status_df_impl
)

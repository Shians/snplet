#' Get donor-level SNP heterozygosity status
#'
#' Returns a long-format data frame with heterozygosity status per SNP and donor.
#'
#' @param x A SNPData object
#' @param min_total_count Minimum total read depth (ref + alt) required per donor to test for heterozygosity (default: 10)
#' @param p_value_threshold P-value threshold for binomial test (default: 0.05). P-values are multiple-testing corrected and SNPs with p < threshold reject monoallelic expression.
#' @param minor_allele_prop Minor allele proportion used as the null threshold for monoallelic expression testing (default: 0.1).
#' @return A tibble with columns: snp_id, gene_name, chrom, pos, strand (if available in snp_info), donor, ref_count, alt_count, total_count, ref_ratio, maf, minor_allele_count, p_val, adj_p_val, tested, zygosity
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
    old_threshold <- logger::log_threshold()
    logger::log_threshold(logger::WARN)
    on.exit(logger::log_threshold(old_threshold), add = TRUE)
    
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
            zygosity = dplyr::case_when(
                tested & !is.na(adj_p_val) & adj_p_val < p_value_threshold ~ "het",
                tested ~ "hom",
                TRUE ~ "unknown"
            )
        )
}

#' @rdname donor_het_status_df
#' @include SNPData-class.R
setMethod(
    "donor_het_status_df",
    signature(x = "SNPData"),
    donor_het_status_df_impl
)

#' Get heterozygous SNP data for a specific donor
#'
#' Filters a SNPData object to a specific donor and only includes SNPs that are 
#' heterozygous for that donor.
#'
#' @param snp_data A SNPData object
#' @param donor Character string specifying the donor to filter for
#' @param ... Additional arguments passed to `donor_het_status_df`
#' @return A filtered SNPData object containing only the specified donor and 
#'   their heterozygous SNPs
#' @export
#'
#' @examples
#' \dontrun{
#' snp_data <- get_example_snpdata()
#' donor_het_data <- get_donor_het_snpdata(snp_data, "donor0")
#' }
get_donor_het_snpdata <- function(snp_data, donor, ...) {
    # Filter to the specified donor first
    donor_data <- filter_barcodes(snp_data, donor == !!donor)
    
    # Get heterozygosity status for all donors (but we only care about our donor)
    het_status <- donor_het_status_df(donor_data, ...)
    
    # Get SNPs that are heterozygous for this donor
    het_snps <- het_status %>%
        dplyr::filter(donor == !!donor, zygosity == "het") %>%
        dplyr::pull(snp_id)
    
    # Filter to heterozygous SNPs
    if (length(het_snps) > 0) {
        filter_snps(donor_data, snp_id %in% het_snps)
    } else {
        # Return empty SNPData object if no heterozygous SNPs
        filter_snps(donor_data, FALSE)
    }
}


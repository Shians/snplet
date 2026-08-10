# Runs the binomial monoallelic-expression test and classifies the result as
# "het"/"hom". `donor_counts` must already be filtered to exactly the
# (snp_id, donor) population meant to be BH-corrected together -- the
# correction scope is a caller decision (donor_het_status_df() and
# infer_zygosity() intentionally test different populations).
.binomial_zygosity_calls <- function(donor_counts, p_value_threshold, minor_allele_prop) {
    if (nrow(donor_counts) == 0) {
        return(
            donor_counts %>%
                dplyr::mutate(
                    minor_allele_count = NA_real_,
                    p_val = NA_real_,
                    adj_p_val = NA_real_,
                    zygosity = NA_character_
                )
        )
    }
    donor_counts %>%
        test_maf(p = minor_allele_prop) %>%
        dplyr::mutate(zygosity = dplyr::if_else(adj_p_val < p_value_threshold, "het", "hom"))
}

#' Get donor-level SNP heterozygosity status
#'
#' Returns a long-format data frame with heterozygosity status per SNP and donor.
#' A \code{(snp_id, donor)} pair with a zygosity call already stored in
#' \code{donor_snp_info(x)} (e.g. from Vireo genotypes read at import, or from
#' \code{\link{infer_zygosity}}) is returned as-is rather than recomputed, so the same
#' pair's zygosity doesn't flip without warning depending on which chromosomes happen to be
#' loaded alongside it. Pairs with no stored call fall back to the binomial test as
#' before. Errors if \code{x} has no active zygosity source at all (see
#' \code{\link{zygosity_source}}); call \code{\link{infer_zygosity}} first.
#'
#' @param x A SNPData object, required.
#' @param min_total_count Numeric (default 10). Minimum total read depth
#'   (ref + alt) required per donor to test for heterozygosity.
#' @param p_value_threshold Numeric, in \code{[0, 1]} (default 0.05).
#'   P-value threshold for the binomial test; p-values are multiple-testing
#'   corrected, and SNPs with p < threshold reject monoallelic expression.
#' @param minor_allele_prop Numeric, in \code{[0, 1]} (default 0.1). Minor
#'   allele proportion used as the null threshold for monoallelic expression
#'   testing.
#' @return A tibble with columns: snp_id, gene_name, chrom, pos, strand (if available in
#'   snp_info), donor, ref_count, alt_count, total_count, ref_ratio, maf,
#'   minor_allele_count, p_val, adj_p_val, tested, zygosity, zygosity_source. For a pair
#'   with a stored call, \code{p_val}/\code{adj_p_val} are \code{NA} (no binomial test was
#'   run) and \code{zygosity_source} names the origin (e.g. \code{"vireo_gt"}); for a
#'   binomial-tested pair \code{zygosity_source} is \code{"binomial"}.
#' @export
#'
#' @examples
#' \dontrun{
#' snp_data <- get_example_snpdata()
#' donor_het_status_df(snp_data)
#' }
#' @rdname donor_het_status_df
setGeneric(
    "donor_het_status_df",
    function(
        x,
        min_total_count = 10,
        p_value_threshold = 0.05,
        minor_allele_prop = 0.1
    ) {
        standardGeneric("donor_het_status_df")
    }
)

donor_het_status_df_impl <- function(
    x,
    min_total_count = 10,
    p_value_threshold = 0.05,
    minor_allele_prop = 0.1
) {
    old_threshold <- logger::log_threshold()
    logger::log_threshold(logger::WARN)
    on.exit(logger::log_threshold(old_threshold), add = TRUE)

    if (is.na(zygosity_source(x))) {
        stop(
            "No zygosity source established for this SNPData object. Call infer_zygosity(x) to ",
            "compute binomial-derived heterozygosity calls, or import with Vireo genotypes ",
            "(import_cellsnp(vireo_folder = ...))."
        )
    }

    stopifnot(min_total_count >= 1)
    stopifnot(p_value_threshold > 0 && p_value_threshold <= 1)
    stopifnot(minor_allele_prop > 0 && minor_allele_prop < 0.5)

    donor_counts <- donor_count_df(x, test_maf = FALSE) %>%
        dplyr::mutate(tested = total_count >= min_total_count)

    stored_calls <- donor_snp_info(x) %>%
        dplyr::select(snp_id, donor, zygosity, zygosity_source) %>%
        dplyr::filter(!is.na(zygosity))

    # Only binomial-test (snp_id, donor) pairs with no trusted stored call.
    to_test <- donor_counts %>%
        dplyr::anti_join(stored_calls, by = c("snp_id", "donor")) %>%
        dplyr::filter(tested)

    tested_counts <- .binomial_zygosity_calls(to_test, p_value_threshold, minor_allele_prop) %>%
        dplyr::select(snp_id, donor, minor_allele_count, p_val, adj_p_val)

    donor_counts %>%
        dplyr::left_join(tested_counts, by = c("snp_id", "donor")) %>%
        dplyr::left_join(stored_calls, by = c("snp_id", "donor")) %>%
        dplyr::mutate(
            zygosity = dplyr::case_when(
                !is.na(zygosity) ~ zygosity,
                tested & !is.na(adj_p_val) & adj_p_val < p_value_threshold ~ "het",
                tested ~ "hom",
                TRUE ~ "unknown"
            ),
            zygosity_source = dplyr::if_else(
                is.na(zygosity_source) & zygosity %in% c("het", "hom"),
                "binomial",
                zygosity_source
            ),
            tested = tested | zygosity %in% c("het", "hom")
        )
}

#' @rdname donor_het_status_df
#' @include SNPData-class.R
setMethod(
    "donor_het_status_df",
    signature(x = "SNPData"),
    donor_het_status_df_impl
)

#' Compute and persist binomial-derived zygosity calls
#'
#' The one entry point for binomial-derived heterozygosity calls: computes the
#' monoallelic-expression binomial test (see \code{\link{test_maf}}) for every
#' \code{(snp_id, donor)} pair meeting \code{min_total_count}, and writes the results into
#' \code{donor_snp_info(x, source = "all")} tagged \code{zygosity_source = "binomial"}.
#' Because \code{donor_snp_info} keys on \code{(snp_id, donor, zygosity_source)}, this can
#' never touch a call from any other source (e.g. \code{"vireo_gt"}) for the same pair;
#' it is purely additive. If \code{x} has no active \code{\link{zygosity_source}} yet, it
#' becomes \code{"binomial"}; otherwise the active source is left untouched, so calling
#' this on an object that already has Vireo genotypes adds a second, independently
#' selectable source without changing what \code{donor_het_status_df()}/\code{assign_xci()}
#' use by default. Switch to it with \code{zygosity_source(x) <- "binomial"}.
#'
#' @param x A SNPData object, required.
#' @param min_total_count Numeric (default 10). Minimum total read depth
#'   (ref + alt) required per donor to test for heterozygosity.
#' @param p_value_threshold Numeric, in \code{[0, 1]} (default 0.05).
#'   P-value threshold for the binomial test; p-values are multiple-testing
#'   corrected, and SNPs with p < threshold reject monoallelic expression.
#' @param minor_allele_prop Numeric, in \code{[0, 1]} (default 0.1). Minor
#'   allele proportion used as the null threshold for monoallelic expression
#'   testing.
#' @return A SNPData object with binomial-derived rows added to \code{donor_snp_info}.
#' @export
#'
#' @examples
#' \dontrun{
#' snp_data <- get_example_snpdata()
#' snp_data <- infer_zygosity(snp_data)
#'
#' # Compare against Vireo genotype calls already present on the same object
#' snp_data_vireo <- import_cellsnp(cellsnp_dir, gene_annotation, vireo_folder = vireo_dir)
#' snp_data_vireo <- infer_zygosity(snp_data_vireo)
#' vireo_calls <- donor_snp_info(snp_data_vireo, source = "vireo_gt")
#' binomial_calls <- donor_snp_info(snp_data_vireo, source = "binomial")
#' dplyr::inner_join(
#'     vireo_calls,
#'     binomial_calls,
#'     by = c("snp_id", "donor"),
#'     suffix = c("_vireo", "_binomial")
#' )
#' }
infer_zygosity <- function(x, min_total_count = 10, p_value_threshold = 0.05, minor_allele_prop = 0.1) {
    old_threshold <- logger::log_threshold()
    logger::log_threshold(logger::WARN)
    on.exit(logger::log_threshold(old_threshold), add = TRUE)

    stopifnot(min_total_count >= 1)
    stopifnot(p_value_threshold > 0 && p_value_threshold <= 1)
    stopifnot(minor_allele_prop > 0 && minor_allele_prop < 0.5)

    testable <- donor_count_df(x, test_maf = FALSE) %>%
        dplyr::filter(total_count >= min_total_count)

    if (nrow(testable) == 0) {
        logger::log_warn(
            "infer_zygosity(): no (snp_id, donor) pair meets min_total_count = {min_total_count}; ",
            "SNPData unchanged."
        )
        return(x)
    }

    calls <- .binomial_zygosity_calls(testable, p_value_threshold, minor_allele_prop) %>%
        dplyr::transmute(
            snp_id,
            donor,
            zygosity_source = "binomial",
            zygosity,
            zygosity_p_val = p_val,
            zygosity_adj_p_val = adj_p_val
        )

    # add_donor_snp_metadata() sets the active zygosity_source slot to
    # "binomial" automatically if x had no active source yet (see
    # .propagate_zygosity_source()); an already-active source (e.g.
    # "vireo_gt") is left untouched.
    add_donor_snp_metadata(x, calls, join_by = c("snp_id", "donor", "zygosity_source"), overwrite = TRUE)
}

#' Get heterozygous SNP data for a specific donor
#'
#' Filters a SNPData object to a specific donor and only includes SNPs that are
#' heterozygous for that donor.
#'
#' @param snp_data A SNPData object, required.
#' @param donor Character scalar, required. Donor to filter for.
#' @param ... Additional arguments passed to \code{donor_het_status_df}.
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

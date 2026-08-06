#' Test for XCI escape significance
#'
#' Performs a one-sided exact beta-binomial test for each row to assess whether the
#' phase-corrected inactive-X allele count is significantly greater than expected
#' under a null escape fraction, modelling the overdispersion in single-cell allelic
#' counts via a beta-binomial rather than binomial null.
#'
#' Unlike a MAF test on raw \code{ref_count}/\code{alt_count}, this operates on counts
#' already reoriented to the inactive allele per cell (see
#' \code{\link{haplotype_expression}}'s \code{active_count}/\code{inactive_count}).
#' Pooling raw ref/alt across cells without that reorientation confounds escape with
#' between-cell heterogeneity in which X is active: a gene under strict XCI in a
#' population split between X1-active and X2-active cells pools to a MAF near 0.5,
#' indistinguishable from a gene that escapes in every cell. Reorienting to
#' active/inactive before pooling removes that confound.
#'
#' @details
#' \code{\link{haplotype_expression}} returns one row per (donor, phased SNP,
#' active-X group) -- multiple rows per gene whenever it has more than one
#' informative SNP, and always two rows per SNP (the X1-active and X2-active
#' cell groups, both already reoriented to active/inactive terms). Testing that
#' table directly over-counts: it multiplies the number of tests (inflating the
#' BH correction) and splits a single gene's evidence across redundant rows
#' instead of pooling it. Collapse to one row per (donor, gene) -- summing
#' \code{active_count}/\code{inactive_count} -- before calling \code{test_escape()};
#' see the second example.
#'
#' The BH correction in \code{adj_p_val} is scoped to whatever rows are passed
#' in a single call. Passing every donor's genes in one call corrects across the
#' whole multi-donor table; calling once per donor (e.g. via
#' \code{dplyr::group_by(donor) \%>\% dplyr::group_modify()}) corrects within each
#' donor independently. Neither is more "correct" in general -- it depends
#' whether donors should be treated as independent experiments or pooled -- but
#' the two give different \code{adj_p_val} for the same gene, so pick
#' deliberately rather than defaulting to whichever grouping the input happens
#' to arrive in.
#'
#' @param x A data.frame (or tibble) with one row per observation, containing the
#'   integer columns \code{active_count} and \code{inactive_count}. Typically
#'   \code{\link{haplotype_expression}} output collapsed to one row per (donor,
#'   gene) -- see Details.
#' @param p Numeric scalar or vector in \code{[0, 1]} giving the null hypothesis
#'   escape fraction (default: 0.10). Recycled against the rows of \code{x} if a
#'   vector, e.g. to supply a different, empirically estimated null per donor
#'   (such as \code{donor_info(x)$xci_median_pi_g}, the median escape fraction
#'   among that donor's informative genes -- a data-driven background/noise
#'   level rather than an assumed constant).
#' @param rho Numeric scalar or vector in \code{[0, 1)} giving the beta-binomial
#'   overdispersion, typically \code{donor_info(x)$xci_rho} as fitted by
#'   \code{\link{assign_xci}} -- a donor-pooled estimate, deliberately not the
#'   EM's per-cell overdispersion (see \code{\link{.fit_pooled_rho_by_donor}}
#'   for why the two differ). Recycled against the rows of \code{x} if a
#'   vector, e.g. to supply a different \code{rho} per donor.
#'
#' @return The input \code{x} with three columns appended:
#'   \describe{
#'     \item{\code{coverage}}{\code{active_count + inactive_count} for each row.}
#'     \item{\code{p_val}}{One-sided exact beta-binomial p-value for observing at
#'       least \code{inactive_count} inactive-allele reads out of
#'       \code{ceiling(coverage)} trials under \code{BetaBinomial(n, p, rho)}.}
#'     \item{\code{adj_p_val}}{Benjamini-Hochberg (BH) adjusted p-value.}
#'   }
#'
#' @export
#' @examples
#' df <- tibble::tibble(
#'     donor = c("donor0", "donor0", "donor1"),
#'     active_count = c(10, 5, 8),
#'     inactive_count = c(1, 6, 2)
#' )
#'
#' # Fixed null escape fraction and overdispersion, shared across all rows
#' test_escape(df, p = 0.10, rho = 0.3)
#'
#' # Donor-specific null and overdispersion, taken from assign_xci()'s fit
#' # (here fabricated; normally donor_info(snp_data) after assign_xci(snp_data))
#' donor_null <- tibble::tibble(
#'     donor = c("donor0", "donor1"),
#'     xci_median_pi_g = c(0.02, 0.05),
#'     xci_rho = c(0.3, 0.4)
#' )
#' df_with_null <- dplyr::left_join(df, donor_null, by = "donor")
#' test_escape(df_with_null, p = df_with_null$xci_median_pi_g, rho = df_with_null$xci_rho)
#'
#' \dontrun{
#' # Full recommended pipeline on real haplotype_expression() output: collapse
#' # the SNP x active-X-group rows to one row per (donor, gene), attach each
#' # donor's empirical null, then BH-correct within each donor separately.
#' snp_data <- assign_xci(snp_data)
#' hap_by_gene <- haplotype_expression(snp_data) %>%
#'     dplyr::summarise(
#'         active_count = sum(active_count),
#'         inactive_count = sum(inactive_count),
#'         .by = c(donor, gene_name)
#'     ) %>%
#'     dplyr::left_join(
#'         dplyr::select(donor_info(snp_data), donor, xci_median_pi_g, xci_rho),
#'         by = "donor"
#'     )
#'
#' escape_result <- hap_by_gene %>%
#'     dplyr::group_by(donor) %>%
#'     dplyr::group_modify(~ test_escape(.x, p = .x$xci_median_pi_g, rho = .x$xci_rho)) %>%
#'     dplyr::ungroup()
#' }
test_escape <- function(x, p = 0.10, rho = 0.05) {
    # validate input
    stopifnot(is(x, "data.frame"))

    req_cols <- c("active_count", "inactive_count")
    missing_cols <- setdiff(req_cols, colnames(x))
    if (length(missing_cols) > 0) {
        missing_list <- paste0(missing_cols, collapse = ", ")
        stop(glue::glue("Missing required columns: {missing_list}"))
    }

    # validate that all counts are non-negative
    negative_active <- which(x$active_count < 0)
    if (length(negative_active) > 0) {
        if (length(negative_active) <= 5) {
            row_list <- paste0(negative_active, collapse = ", ")
        } else {
            row_list <- paste0(paste0(head(negative_active, 5), collapse = ", "), ", ...")
        }
        stop(glue::glue("Invalid data: active_count < 0 in row(s): {row_list}"))
    }

    negative_inactive <- which(x$inactive_count < 0)
    if (length(negative_inactive) > 0) {
        if (length(negative_inactive) <= 5) {
            row_list <- paste0(negative_inactive, collapse = ", ")
        } else {
            row_list <- paste0(paste0(head(negative_inactive, 5), collapse = ", "), ", ...")
        }
        stop(glue::glue("Invalid data: inactive_count < 0 in row(s): {row_list}"))
    }

    inactive_count <- x$inactive_count
    coverage <- ceiling(x$active_count + x$inactive_count)

    p_val <- betabinom_test(inactive_count, coverage, p, rho, alternative = "greater")

    result <- x %>%
        dplyr::mutate(
            coverage = coverage,
            p_val = p_val,
            adj_p_val = p.adjust(p_val, method = "BH")
        )

    return(result)
}

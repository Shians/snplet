#' Test for XCI escape significance
#'
#' Performs a one-sided exact beta-binomial test for each row to assess whether the
#' phase-corrected inactive-X allele count is significantly greater than expected
#' under a null escape fraction, modelling the overdispersion in single-cell allelic
#' counts via a beta-binomial rather than binomial null.
#'
#' This operates on counts already reoriented to the inactive allele per cell
#' (see \code{\link{haplotype_expression}}'s \code{active_count}/\code{inactive_count}),
#' unlike a MAF test on raw \code{ref_count}/\code{alt_count}. Pooling raw ref/alt
#' across cells without that reorientation confounds escape with between-cell
#' heterogeneity in which X is active: a gene under strict XCI in a population
#' split between X1-active and X2-active cells pools to a MAF near 0.5,
#' indistinguishable from a gene that escapes in every cell. Reorienting to
#' active/inactive before pooling removes that confound.
#'
#' @details
#' Passing the SNPData object itself is the recommended path: it takes the
#' counts at the right grain, the null from the donor's own fit, and the
#' multiplicity correction within that donor, none of which the data.frame
#' method can infer from a bare table. The data.frame method exists for
#' hand-built or externally derived counts.
#'
#' Every count must be pooled to one row per gene before testing. Two failure
#' modes are worth naming, both of which the SNPData method avoids by
#' construction:
#' \itemize{
#'   \item Testing the two active-X groups as separate rows doubles the number
#'     of tests BH corrects over and splits a gene's evidence between two
#'     half-coverage rows, so it can reach significance in one group and not
#'     the other purely by which cells fell where.
#'     \code{\link{haplotype_expression}} pools them by default; only
#'     \code{by_active_x = TRUE} splits them.
#'   \item Summing a gene's SNPs counts a read spanning several of them once
#'     per SNP, inflating well-covered multi-SNP genes. Never sum the per-SNP
#'     rows from \code{by_snp = TRUE}; the default output already elects one
#'     representative SNP per gene, and
#'     \code{\link{haplotype_expression_by_molecule}} is how to use a gene's
#'     other SNPs without double-counting.
#' }
#'
#' Rows with zero coverage cannot be tested and are given \code{NA} rather
#' than a p-value of 1, which keeps them out of the BH denominator instead of
#' inflating it.
#'
#' @inheritSection assign_xci Phase is inferred from expression, not genotyped
#'
#' @param x A SNPData object that \code{\link{assign_xci}} has annotated, or a
#'   data.frame with integer columns \code{active_count} and
#'   \code{inactive_count} and one row per gene (required).
#' @param p Numeric, in \code{[0, 1]}, scalar or one value per row of \code{x}
#'   (default \code{NULL}). Null hypothesis escape fraction. \code{NULL} means
#'   each donor's own \code{xci_median_pi_g} for a SNPData, and 0.10 for a
#'   data.frame, which has no donor fit to draw on.
#' @param rho Numeric, in \code{[0, 1)}, scalar or one value per row of
#'   \code{x} (default \code{NULL}). Beta-binomial overdispersion. \code{NULL}
#'   means each donor's own \code{xci_rho} for a SNPData, and 0.05 for a
#'   data.frame. \code{xci_rho} is a donor-pooled estimate fitted at this
#'   function's own grain, not the EM's per-cell overdispersion; see
#'   \code{\link{.fit_pooled_rho_by_donor}}.
#' @param by_donor Logical (default \code{TRUE}). Correct within each donor
#'   separately, treating donors as independent experiments, which matches the
#'   per-donor \code{p} and \code{rho}. \code{FALSE} corrects across the whole
#'   table at once, appropriate when donors are replicates of one experiment.
#'   Ignored when \code{x} carries no \code{donor} column.
#'
#' @return For a data.frame, the input with three columns appended; for a
#'   SNPData, the gene-level count table those same columns are appended to
#'   (see \code{\link{haplotype_expression}}), plus \code{count_source}.
#'   \describe{
#'     \item{\code{coverage}}{\code{active_count + inactive_count} for each row.}
#'     \item{\code{p_val}}{One-sided exact beta-binomial p-value for observing at
#'       least \code{inactive_count} inactive-allele reads out of
#'       \code{ceiling(coverage)} trials under \code{BetaBinomial(n, p, rho)};
#'       \code{NA} for an untestable row (zero coverage, or a donor whose fit
#'       yielded no \code{p}/\code{rho}).}
#'     \item{\code{adj_p_val}}{Benjamini-Hochberg (BH) adjusted p-value, over
#'       the testable rows of the donor (or of the whole table when
#'       \code{by_donor = FALSE}).}
#'     \item{\code{count_source}}{SNPData only: \code{"molecule"} or
#'       \code{"snp"}, naming where the counts came from.}
#'   }
#'
#' @family X-chromosome inactivation functions
#' @export
#' @examples
#' df <- tibble::tibble(
#'     donor = c("donor0", "donor0", "donor1"),
#'     active_count = c(10, 5, 8),
#'     inactive_count = c(1, 6, 2)
#' )
#'
#' # Fixed null escape fraction and overdispersion, shared across all rows.
#' # donor0's two genes are corrected together, donor1's on its own.
#' test_escape(df, p = 0.10, rho = 0.3)
#'
#' # One correction across every donor's genes instead
#' test_escape(df, p = 0.10, rho = 0.3, by_donor = FALSE)
#'
#' # Donor-specific null and overdispersion, taken from assign_xci()'s fit
#' # (here fabricated; the SNPData method reads these off the object itself)
#' donor_null <- tibble::tibble(
#'     donor = c("donor0", "donor1"),
#'     xci_median_pi_g = c(0.02, 0.05),
#'     xci_rho = c(0.3, 0.4)
#' )
#' df_with_null <- dplyr::left_join(df, donor_null, by = "donor")
#' test_escape(df_with_null, p = df_with_null$xci_median_pi_g, rho = df_with_null$xci_rho)
#'
#' \dontrun{
#' # The recommended path: counts, null and correction all taken from the object
#' snp_data <- assign_xci(snp_data)
#' escape_result <- test_escape(snp_data)
#'
#' # Uses read-backed molecule counts automatically once they are available
#' snp_data <- add_molecule_phase(snp_data, bam_files = c(lib1 = "lib1.bam"))
#' escape_result <- test_escape(snp_data)
#' }
setGeneric(
    "test_escape",
    function(x, p = NULL, rho = NULL, by_donor = TRUE) {
        standardGeneric("test_escape")
    }
)

#' @rdname test_escape
setMethod("test_escape", signature(x = "data.frame"), function(x, p = NULL, rho = NULL, by_donor = TRUE) {
    # The data.frame method has no donor fit to read a null off, so it falls
    # back to fixed values rather than erroring: it is the method for counts
    # that came from somewhere other than assign_xci().
    p <- p %||% 0.10
    rho <- rho %||% 0.05

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

    # Recycle up front so that subsetting to the testable rows below keeps p and
    # rho aligned with the counts, whether they arrived as a scalar or per row.
    p <- rep_len(p, nrow(x))
    rho <- rep_len(rho, nrow(x))

    # A row with no reads carries no evidence either way, and a donor whose fit
    # produced no p or rho has no null to test against. Both get NA rather than
    # a p-value of 1: p.adjust() drops NAs from the denominator, so an
    # untestable row no longer costs the testable ones power.
    testable <- coverage > 0 & !is.na(p) & !is.na(rho)
    p_val <- rep(NA_real_, nrow(x))
    if (any(testable)) {
        p_val[testable] <- betabinom_test(
            inactive_count[testable],
            coverage[testable],
            p[testable],
            rho[testable],
            alternative = "greater"
        )
    }

    result <- dplyr::mutate(x, coverage = coverage, p_val = p_val)

    # Correcting within a donor keeps the correction on the same footing as the
    # per-donor p and rho above: one experiment per donor.
    if (by_donor && "donor" %in% colnames(result)) {
        result <- dplyr::mutate(result, adj_p_val = p.adjust(p_val, method = "BH"), .by = donor)
    } else {
        result <- dplyr::mutate(result, adj_p_val = p.adjust(p_val, method = "BH"))
    }

    result
})

#' @rdname test_escape
#' @include SNPData-class.R
setMethod("test_escape", signature(x = "SNPData"), function(x, p = NULL, rho = NULL, by_donor = TRUE) {
    if (!.has_xci_diagnostics(x)) {
        stop("No stored XCI diagnostics found. Run assign_xci(x) first.")
    }

    counts <- .escape_counts(x)

    donor_fit <- donor_info(x)
    missing_fit <- setdiff(c("xci_median_pi_g", "xci_rho"), colnames(donor_fit))
    if (length(missing_fit) > 0 && (is.null(p) || is.null(rho))) {
        stop(
            "donor_info(x) has no ",
            paste(missing_fit, collapse = " or "),
            " to take the null from; re-run assign_xci(x), or pass p and rho explicitly."
        )
    }

    # Left join rather than filter: a donor whose fit yielded no null keeps its
    # genes in the output with NA p_val, so it reads as untested rather than
    # silently vanishing.
    if ("donor" %in% colnames(counts)) {
        counts <- dplyr::left_join(
            counts,
            dplyr::select(donor_fit, dplyr::any_of(c("donor", "xci_median_pi_g", "xci_rho"))),
            by = "donor"
        )
    } else {
        # No donor column means a single-donor object; its sole fit applies to
        # every row.
        counts$xci_median_pi_g <- dplyr::first(donor_fit$xci_median_pi_g)
        counts$xci_rho <- dplyr::first(donor_fit$xci_rho)
    }

    unfitted <- unique(counts$donor[is.na(counts$xci_median_pi_g) | is.na(counts$xci_rho)])
    if (length(unfitted) > 0) {
        logger::log_warn(
            "No stored null for donor(s) {paste(unfitted, collapse = ', ')}; ",
            "their genes are returned untested (NA p_val)"
        )
    }

    result <- test_escape(
        dplyr::select(counts, -xci_median_pi_g, -xci_rho),
        p = p %||% counts$xci_median_pi_g,
        rho = rho %||% counts$xci_rho,
        by_donor = by_donor
    )
    result
})

#' Take gene-level escape counts from whichever source the object carries
#'
#' Prefers read-backed molecule counts, which count a molecule spanning several
#' of a gene's het SNPs once instead of once per SNP and so use every SNP's
#' evidence, and falls back to the per-SNP counts of
#' \code{\link{haplotype_expression}}, which elect a single representative SNP
#' per gene. Both are returned at one row per (donor, gene).
#'
#' Which source was used is logged and returned as \code{count_source} rather
#' than left implicit, since it depends on the object's state: the same call on
#' the same data before and after \code{\link{add_molecule_phase}} gives
#' different counts, and the column is what makes that visible in the result.
#'
#' @param x A SNPData object with stored XCI diagnostics, required.
#'
#' @return The chosen count table with a \code{count_source} column of
#'   \code{"molecule"} or \code{"snp"}.
#'
#' @keywords internal
.escape_counts <- function(x) {
    donor_snp_info <- donor_snp_info(x)
    has_molecules <- !is.null(attr(x, "molecule_calls")) &&
        all(c("phase_block", "allele_on_x1") %in% colnames(donor_snp_info)) &&
        nrow(snp_gene_map(x)) > 0

    if (has_molecules) {
        logger::log_info("Counting escape from read-backed molecules")
        return(dplyr::mutate(haplotype_expression_by_molecule(x), count_source = "molecule"))
    }

    logger::log_info("No molecule phase stored; counting escape from per-SNP reads")
    dplyr::mutate(haplotype_expression(x), count_source = "snp")
}

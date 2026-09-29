#' Extract per-cell active-X assignments
#'
#' Pulls the stored X-chromosome inactivation call and posterior for every cell
#' out of a SNPData object's barcode metadata. Requires that
#' \code{\link{assign_xci}} or \code{\link{assign_xci_by_clonotype}}
#' has been run first; errors otherwise.
#'
#' @inheritSection assign_xci Phase is inferred from expression, not genotyped
#'
#' @param x A SNPData object, required, that had XCI diagnostics stored by
#'   \code{\link{assign_xci}} or \code{\link{assign_xci_by_clonotype}}.
#'
#' @return A tibble with one row per cell and columns \code{cell_id}
#'   (character), \code{donor} (character; present only when the object carries
#'   donor information), \code{active_x} (character "X1"/"X2", or \code{NA}
#'   for cells below the confidence threshold), and \code{xci_post_X1_active}
#'   (numeric posterior probability that X1 is the active chromosome).
#'
#' @family X-chromosome inactivation functions
#' @export
setGeneric("xci_assignments", function(x) standardGeneric("xci_assignments"))

#' @rdname xci_assignments
#' @include SNPData-class.R
setMethod("xci_assignments", signature(x = "SNPData"), function(x) {
    barcode_info <- barcode_info(x)
    if (!"active_x" %in% colnames(barcode_info)) {
        stop("No stored XCI diagnostics found. Run assign_xci(x) first.")
    }
    barcode_info %>%
        dplyr::select(cell_id, dplyr::any_of("donor"), active_x, xci_post_X1_active)
})

#' Extract SNP haplotypes from a SNPData object with stored XCI diagnostics
#'
#' Returns the inferred phase and escape fraction for each SNP in the final
#' gene set used by the EM model. Phase is donor-specific (a SNP retained in
#' several donors carries a distinct X1-allele in each), so the result has one
#' row per SNP and donor rather than collapsing them to a single value.
#'
#' Despite the name, these are not genotyped haplotypes: they are the model's
#' expression-derived assignment of alleles to two clusters. See the section
#' below before treating \code{allele_on_x1} as chromosomal phase.
#'
#' @inheritSection assign_xci Phase is inferred from expression, not genotyped
#'
#' @param x A SNPData object, required, that had XCI diagnostics stored by
#'   \code{\link{assign_xci}} or \code{\link{assign_xci_by_clonotype}}.
#'
#' @return A tibble with one row per informative SNP and donor, with columns
#'   \code{snp_id} (character), \code{gene_name} (character; present only when
#'   the object carries gene annotation), \code{donor} (character),
#'   \code{allele_on_x1} (character, "REF" or "ALT": the allele the model
#'   assigned to the X1 cluster in that donor, i.e. the one expressed by
#'   X1-active cells, not a genotyped haplotype), and \code{escape_fraction}
#'   (numeric estimated fraction of reads from the minority allele in that
#'   donor, bounded below 0.5; see the section below).
#'
#' @family X-chromosome inactivation functions
#' @importFrom tidyr separate_longer_delim
#' @export
setGeneric("xci_haplotypes", function(x) standardGeneric("xci_haplotypes"))

#' @rdname xci_haplotypes
#' @include SNPData-class.R
setMethod("xci_haplotypes", signature(x = "SNPData"), function(x) {
    if (!.has_xci_diagnostics(x)) {
        stop("No stored XCI diagnostics found. Run assign_xci(x) first.")
    }
    snp_info <- snp_info(x)
    donor_snp_info(x) %>%
        dplyr::filter(xci_informative) %>%
        dplyr::left_join(dplyr::select(snp_info, snp_id, dplyr::any_of("gene_name")), by = "snp_id") %>%
        dplyr::select(
            snp_id,
            dplyr::any_of("gene_name"),
            donor,
            allele_on_x1,
            escape_fraction = xci_escape_fraction
        )
})

#' Write XCI diagnostics from a fit into a SNPData object
#'
#' Promotes the fit diagnostics into the object's indexable metadata slots so
#' they survive subsetting: barcode metadata gains \code{active_x} and
#' \code{xci_post_X1_active}; \code{donor_snp_info} gains one row per phased
#' SNP x donor with \code{xci_informative} (whether the gene drove active-X
#' calling), \code{allele_on_x1} and \code{xci_escape_fraction}; \code{donor_info}
#' gains one row per donor with \code{xci_median_pi_g} (median escape fraction
#' among informative genes, a per-donor empirical background/noise level),
#' \code{xci_rho}, and \code{xci_skew} (the EM's fitted X1-active prior: the
#' model's own estimate of what fraction of this donor's cells have X1
#' active, i.e. the quantified XCI skew, fit jointly with phase and escape
#' rather than read off the confidence-thresholded hard calls afterwards).
#' \code{donor_info} also gains fit-quality diagnostics: \code{xci_n_het_genes}
#' (genes entering the EM), \code{xci_n_informative} (genes that drove active-X
#' calling), \code{xci_frac_assigned} (fraction of modelling units, cells or
#' clonotypes, meeting the confidence threshold) and \code{xci_minority_frac}
#' (fraction of assigned units on the minority X). Last, \code{xci_reliable} and
#' \code{xci_flags} label donors failing any check in
#' \code{\link{.flag_unreliable_xci_fits}}, which are also logged as warnings
#' and excluded from escape results by default.
#' For a clonotype-level fit the per-cell projection is used for barcode
#' annotation.
#'
#' \code{xci_rho} is \emph{not} the EM's per-cell beta-binomial overdispersion
#' (\code{xci_result$rho} inside \code{.infer_xci()}); that value is fit
#' from individual cells' small read counts and, applied directly to
#' donor-pooled counts (as \code{\link{test_escape}} needs), implicitly assumes
#' complete correlation across all of a donor's cells, which is far too
#' conservative once thousands of cells are pooled (see
#' \code{\link{.fit_pooled_rho_by_donor}}). \code{xci_rho} is instead refit
#' directly at the donor-pooled level, against the population of informative
#' (non-escaping) genes' pooled \code{active_count}/\code{inactive_count} from
#' \code{\link{haplotype_expression}}, whose default grain is the aggregation
#' level \code{test_escape} actually operates at.
#'
#' @keywords internal
.store_xci_fit <- function(x, fit) {
    donor_fits <- purrr::keep(fit, ~ !is.null(.x))
    logger::log_info("Storing XCI diagnostics for {length(donor_fits)} donor(s)")

    barcode_diag <- purrr::map(donor_fits, function(f) {
        assignments <- if (!is.null(f$cell_assignments)) f$cell_assignments else f$assignments
        assignments %>%
            dplyr::transmute(
                cell_id,
                # assignment already names the active X; "unassigned" means no call.
                active_x = dplyr::na_if(assignment, "unassigned"),
                xci_post_X1_active = post_X1_active,
                # Record the unit the model was fit on so downstream plotting can
                # aggregate cells back to clonotypes when appropriate.
                xci_fit_unit = if (is.null(f$unit)) "cell" else f$unit
            )
    }) %>%
        dplyr::bind_rows()

    snp_diag <- purrr::map(donor_fits, function(f) {
        f$haplotypes %>%
            dplyr::transmute(
                snp_id,
                donor = f$donor,
                zygosity_source = zygosity_source(x),
                xci_informative,
                allele_on_x1,
                xci_escape_fraction = escape_fraction
            )
    }) %>%
        dplyr::bind_rows()

    donor_diag <- purrr::map(donor_fits, function(f) {
        tibble::tibble(
            donor = f$donor,
            xci_median_pi_g = f$median_pi_g,
            xci_skew = f$skew,
            xci_n_het_genes = f$n_het_genes,
            xci_n_informative = f$n_informative,
            xci_frac_assigned = f$frac_assigned,
            xci_minority_frac = f$minority_frac
        )
    }) %>%
        dplyr::bind_rows()

    if (nrow(barcode_diag) > 0) {
        x <- add_barcode_metadata(x, barcode_diag, join_by = "cell_id", overwrite = TRUE)
    }
    if (nrow(snp_diag) > 0) {
        x <- add_donor_snp_metadata(x, snp_diag, join_by = c("snp_id", "donor", "zygosity_source"), overwrite = TRUE)
    }
    if (nrow(donor_diag) > 0) {
        x <- add_donor_metadata(x, donor_diag, join_by = "donor", overwrite = TRUE)
    }

    # xci_rho depends on xci_median_pi_g and haplotype_expression(), both of
    # which require the diagnostics just written above, so this must run last.
    logger::log_info("Fitting donor-pooled overdispersion (xci_rho)")
    pooled_rho <- .fit_pooled_rho_by_donor(x)
    if (nrow(pooled_rho) > 0) {
        x <- add_donor_metadata(x, pooled_rho, join_by = "donor", overwrite = TRUE)
    }
    if (nrow(donor_diag) > 0) {
        reliability <- donor_diag %>%
            dplyr::left_join(pooled_rho, by = "donor") %>%
            .flag_unreliable_xci_fits()
        x <- add_donor_metadata(x, reliability, join_by = "donor", overwrite = TRUE)
    }
    logger::log_info("XCI diagnostics stored")

    x
}

# Majority-X fraction above which a donor's XCI results are labelled unreliable.
# Beyond it the fitted prior shifts a minority cell's log-odds by more than
# about three quarters of a read's evidence, so single-read minority cells go
# uncalled, and an RNA-derived het call can no longer separate an inactivated
# het from background noise (pooled minor-allele fraction about m + 0.04).
.xci_max_reliable_skew <- 0.9

#' Label donors whose XCI fit is unlikely to be valid
#'
#' Two kinds of donor are labelled. A heavily skewed donor (majority X above
#' \code{.xci_max_reliable_skew}) undercalls its minority-X cells, and its
#' zygosity calls, all derived from RNA, miss het SNPs in X-inactivated genes;
#' its results are unreliable without DNA genotypes. A failed fit is otherwise
#' stored as if it were valid, the typical cause being a skewed donor left with
#' only escape genes, which express both alleles in every cell and so cannot
#' separate cells. Logs one warning per labelled donor naming every check it
#' failed.
#'
#' @param donor_diag A data frame with one row per donor and columns
#'   \code{donor}, \code{xci_n_informative}, \code{xci_skew},
#'   \code{xci_minority_frac}, \code{xci_median_pi_g} and \code{xci_rho}.
#'
#' @return A tibble with columns \code{donor}, \code{xci_reliable} (logical)
#'   and \code{xci_flags} (the failed checks joined by \code{"; "}, \code{NA}
#'   when none failed).
#'
#' @keywords internal
.flag_unreliable_xci_fits <- function(donor_diag) {
    reasons <- purrr::pmap(dplyr::select(donor_diag, -donor), .xci_fit_failure_reasons)
    flags <- tibble::tibble(
        donor = donor_diag$donor,
        xci_reliable = lengths(reasons) == 0,
        xci_flags = purrr::map_chr(reasons, ~ if (length(.x) == 0) NA_character_ else paste(.x, collapse = "; "))
    )
    purrr::pwalk(dplyr::filter(flags, !xci_reliable), function(donor, xci_flags, ...) {
        logger::log_warn(
            "[{donor}] XCI fit is unreliable: {xci_flags}. Its escape results are highly unreliable without DNA ",
            "genotypes, and haplotype_expression(), test_escape() and as_escape_experiment() exclude it unless ",
            "include_unreliable = TRUE. See donor_info(x)$xci_flags."
        )
    })
    flags
}

#' Names the fit-quality checks one donor's XCI fit fails
#'
#' @keywords internal
.xci_fit_failure_reasons <- function(
    xci_n_informative,
    xci_skew,
    xci_minority_frac,
    xci_median_pi_g,
    xci_rho,
    ...
) {
    majority_frac <- max(xci_skew, 1 - xci_skew)
    too_few_genes <- xci_n_informative < 10
    one_sided <- is.na(xci_minority_frac) || xci_minority_frac == 0
    extreme_skew <- majority_frac > .xci_max_reliable_skew
    escape_dominated <- !is.na(xci_median_pi_g) && xci_median_pi_g > 0.2
    no_rho <- is.na(xci_rho)

    failed <- c(too_few_genes, one_sided, extreme_skew, escape_dominated, no_rho)
    messages <- c(
        "only {xci_n_informative} informative genes (< 10)",
        "no assigned cells on the minority X",
        paste0(
            "skew {round(majority_frac, 2)} (> {.xci_max_reliable_skew}): minority-X cells are undercalled ",
            "and RNA-derived het calls miss X-inactivated genes"
        ),
        "median escape fraction {round(xci_median_pi_g, 2)} (> 0.2), so the het SNPs are mostly escape genes",
        "pooled overdispersion (xci_rho) could not be fitted"
    )
    purrr::map_chr(messages[failed], ~ glue::glue(.x))
}

#' Drop donors labelled unreliable by the XCI fit
#'
#' @param x A SNPData object.
#' @param include_unreliable Logical. \code{TRUE} returns \code{x} unchanged.
#'
#' @return \code{x} without the cells of donors whose \code{xci_reliable} is
#'   \code{FALSE}. Unchanged when the label has not been stored.
#'
#' @keywords internal
.drop_unreliable_xci_donors <- function(x, include_unreliable) {
    donor_info <- donor_info(x)
    if (include_unreliable || !"xci_reliable" %in% colnames(donor_info)) {
        return(x)
    }
    unreliable <- donor_info$donor[donor_info$xci_reliable %in% FALSE]
    if (length(unreliable) == 0) {
        return(x)
    }
    cell_donor <- barcode_info(x)$donor
    if (all(cell_donor %in% unreliable)) {
        stop(
            "Every donor's XCI fit is labelled unreliable (see donor_info(x)$xci_flags). ",
            "Pass include_unreliable = TRUE to use them anyway."
        )
    }
    logger::log_warn(
        "Excluding {length(unreliable)} donor(s) with unreliable XCI fits: {paste(unreliable, collapse = ', ')}. ",
        "See donor_info(x)$xci_flags, or pass include_unreliable = TRUE."
    )
    x[, !cell_donor %in% unreliable]
}

#' Fit donor-pooled beta-binomial overdispersion for escape testing
#'
#' Estimates \code{rho} at the same aggregation level \code{\link{test_escape}}
#' operates at (donor-pooled gene counts), rather than reusing the EM's
#' per-cell \code{rho}. Restricts to genes that passed
#' \code{.filter_uninformative_genes()} (the trusted, genuinely
#' non-escaping population, same population \code{xci_median_pi_g} summarises)
#' and fits a single scalar \code{rho} per donor via exact beta-binomial MLE
#' against that donor's \code{xci_median_pi_g} as a fixed null probability:
#' how much do these genes' pooled inactive-read fractions actually vary
#' around the shared background rate, beyond binomial sampling.
#'
#' @param x A SNPData object with \code{active_x}, \code{allele_on_x1},
#'   \code{xci_informative} and \code{xci_median_pi_g} already written (i.e.
#'   called from \code{\link{.store_xci_fit}} after those are stored).
#'
#' @return A tibble with columns \code{donor} and \code{xci_rho}. A donor with
#'   fewer than two informative genes (not enough points to estimate spread)
#'   gets \code{NA}.
#'
#' @keywords internal
.fit_pooled_rho_by_donor <- function(x) {
    if (!.has_xci_diagnostics(x)) {
        return(tibble::tibble(donor = character(0), xci_rho = double(0)))
    }

    logger::log_info("Aggregating haplotype expression for pooled rho fitting")
    # The default grain is already one row per (donor, gene): it elects one
    # representative SNP per gene rather than summing across them (a read
    # spanning several of a gene's het SNPs would otherwise be counted once per
    # SNP, inflating well-covered multi-SNP genes and skewing the spread rho is
    # fitted to), and pools the two active-X groups, which are disjoint cell
    # sets. That is the grain test_escape() tests at, which is the whole point
    # of fitting rho here rather than reusing the EM's per-cell value.
    # Reliability is labelled after, and partly from, this fit.
    hap_by_gene <- haplotype_expression(x, xci_informative_only = TRUE, include_unreliable = TRUE) %>%
        dplyr::left_join(dplyr::select(donor_info(x), donor, xci_median_pi_g), by = "donor") %>%
        dplyr::filter(!is.na(xci_median_pi_g))

    if (nrow(hap_by_gene) == 0) {
        return(tibble::tibble(donor = character(0), xci_rho = double(0)))
    }

    logger::log_info("Fitting pooled rho for {dplyr::n_distinct(hap_by_gene$donor)} donor(s)")
    hap_by_gene %>%
        dplyr::summarise(
            xci_rho = if (dplyr::n() >= 2) {
                .fit_pooled_rho(inactive_count, active_count + inactive_count, dplyr::first(xci_median_pi_g))
            } else {
                NA_real_
            },
            .by = donor
        )
}

#' Exact beta-binomial MLE of rho against a fixed null probability
#'
#' @param x Integer vector of pooled inactive-allele counts, one per gene.
#' @param n Integer vector of pooled coverage, one per gene.
#' @param p Numeric scalar, the fixed null probability (e.g. \code{xci_median_pi_g}).
#' @param rho_bounds Numeric length-2 vector, the search interval for \code{rho}.
#'
#' @keywords internal
.fit_pooled_rho <- function(x, n, p, rho_bounds = c(1e-4, 0.999)) {
    neg_ll <- function(rho) {
        -sum(VGAM::dbetabinom(x, size = n, prob = p, rho = rho, log = TRUE))
    }
    stats::optimize(neg_ll, interval = rho_bounds)$minimum
}

#' Whether a SNPData object carries stored XCI fit diagnostics
#'
#' \code{\link{.store_xci_fit}} is the only writer of \code{xci_informative},
#' so at least one \code{TRUE} in \code{donor_snp_info} is a reliable signal
#' that \code{assign_xci}/\code{assign_xci_by_clonotype} has run, even though
#' the column itself may already exist (as \code{NA}) from zygosity calling.
#'
#' @keywords internal
.has_xci_diagnostics <- function(x) {
    donor_snp_info <- donor_snp_info(x)
    "xci_informative" %in% colnames(donor_snp_info) && any(donor_snp_info$xci_informative, na.rm = TRUE)
}

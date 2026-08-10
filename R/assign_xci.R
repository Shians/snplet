#' Assign the active X chromosome to cells
#'
#' Identifies which X chromosome is active in female cells based on allelic
#' imbalance at heterozygous SNPs using an Expectation-Maximisation (EM)
#' algorithm with a beta-binomial likelihood.
#'
#' @details
#' This function infers which X chromosome remains active in each cell by
#' fitting an EM model to allelic read counts at heterozygous SNPs on the X
#' chromosome. X-chromosome inactivation (XCI) is the dosage-compensation
#' mechanism this reconstructs: in female mammals, one of the two X
#' chromosomes is randomly silenced in each cell.
#'
#' The algorithm works as follows:
#' \enumerate{
#'   \item Filters to heterozygous SNPs on the X chromosome for each donor
#'   \item Selects the SNP with highest coverage per gene to avoid redundancy
#'   \item Removes outlier genes with atypical allelic skew
#'   \item Runs a beta-binomial EM algorithm with multiple random initialisations
#'   \item Assigns cells to X1 or X2 based on posterior probability
#' }
#'
#' Cells that do not meet \code{confidence_threshold} receive \code{NA} in
#' the \code{active_x} column.
#'
#' @section Phase is inferred from expression, not genotyped:
#' snplet infers phase from scRNA-seq allelic counts alone: no DNA-based
#' haplotype is ever observed, and that limitation propagates into every
#' downstream result.
#'
#' \code{"X1"} and \code{"X2"} name two expression-defined haplotype
#' clusters, not identified physical chromosomes: exchangeable labels, fixed
#' only by the convention that X1 is the larger active-X group, and not
#' comparable across donors. \code{allele_on_x1} records which allele the
#' model assigned to the X1 cluster in one donor; it is not a genotyped
#' haplotype.
#'
#' The active allele is \emph{defined} as the majority-expressed one: the
#' EM's phase step orients each gene so the silenced allele is the minority,
#' which bounds the per-gene escape fraction below 0.5. Three consequences
#' follow:
#' \itemize{
#'   \item Escape above 0.5 cannot be represented; the measurable range runs
#'     from 0 (strict inactivation) to 0.5 (complete escape), and 0.5 is the
#'     maximum, not the midpoint.
#'   \item A gene transcribed mainly from the inactive X is stored with
#'     inverted phase and surfaces as unusually \emph{clean} inactivation
#'     rather than as an escapee. \code{XIST} is the textbook case.
#'   \item "active"/"inactive" are majority/minority labels, coinciding with
#'     the true active and inactive X only where canonical inactivation
#'     holds (precisely what an escape analysis is questioning).
#' }
#'
#' Separating these cases needs phase independent of expression (DNA-based
#' genotyping, or trio/population phasing), which this package does not
#' have. \code{\link{add_molecule_phase}} supplies true physical linkage
#' between SNPs co-observed on one molecule, but a molecule is a single
#' transcript: it cannot link across genes, and its blocks are still
#' oriented to X1/X2 using the expression-derived fit, so it refines phase
#' within a gene without escaping this constraint.
#'
#' @param x A SNPData object, required, with donor assignments and
#'   heterozygosity information for the X chromosome.
#' @param n_inits Integer (default 10). Number of random initialisations for
#'   the EM algorithm; the run with the highest log-likelihood is returned.
#' @param confidence_threshold Numeric, in \code{[0, 1]} (default 0.95).
#'   Posterior probability threshold for hard assignment: cells whose
#'   posterior probability that a given X is the active one reaches
#'   \code{confidence_threshold} are assigned that X; cells where neither X
#'   reaches the threshold receive \code{NA}.
#'
#' @return SNPData object with an additional \code{active_x} column in
#'   barcode metadata, with values "X1" or "X2" indicating the inferred
#'   active X chromosome state. Cells that do not meet the confidence
#'   threshold receive \code{NA}. The full XCI diagnostics are also written
#'   into the object's metadata slots, so the result can be
#'   passed directly to \code{\link{plot_xci_heatmap}},
#'   \code{\link{xci_assignments}}, and \code{\link{xci_haplotypes}}.
#'
#' @family X-chromosome inactivation functions
#' @export
#'
#' @examples
#' \dontrun{
#' # Assign the active X chromosome to cells
#' snp_data <- assign_xci(snp_data)
#'
#' # View results
#' barcode_info(snp_data) %>%
#'   count(donor, active_x)
#'
#' # Diagnostics are stored, so plotting works directly
#' plot_xci_heatmap(snp_data, donor = "donor1")
#' }
setGeneric(
    "assign_xci",
    function(
        x,
        n_inits = 10,
        confidence_threshold = 0.95
    ) {
        standardGeneric("assign_xci")
    }
)

#' @rdname assign_xci
#' @include SNPData-class.R
setMethod(
    "assign_xci",
    signature(x = "SNPData"),
    function(x, n_inits = 10, confidence_threshold = 0.95) {
        .fit_xci(
            x,
            n_inits = n_inits,
            confidence_threshold = confidence_threshold,
            by = "cell"
        )
    }
)

#' Assign the active X chromosome to cells by clonotype
#'
#' Identifies which X chromosome is active in female cells based on allelic
#' imbalance at heterozygous SNPs, aggregating counts to the clonotype level
#' before running the EM model and projecting assignments back to individual cells.
#'
#' @details
#' This function is similar to \code{\link{assign_xci}} but runs the
#' beta-binomial EM on clonotype-aggregated counts rather than per-cell counts.
#' This approach:
#' \itemize{
#'   \item Increases statistical power by aggregating counts across clonally related cells
#'   \item Reduces noise from sparse cell-level data
#'   \item Is biologically motivated since clonally related cells should share X-inactivation state
#' }
#'
#' The algorithm works as follows:
#' \enumerate{
#'   \item Filters to heterozygous SNPs on the X chromosome for each donor
#'   \item Selects the SNP with highest coverage per gene to avoid redundancy
#'   \item Aggregates ALT and REF counts by clonotype
#'   \item Removes outlier genes with atypical allelic skew
#'   \item Runs a beta-binomial EM algorithm with multiple random initialisations
#'   \item Assigns clonotypes to X1 or X2 based on posterior probability
#'   \item Projects clonotype assignments back to individual cells
#' }
#'
#' @inheritSection assign_xci Phase is inferred from expression, not genotyped
#'
#' @inheritParams assign_xci
#' @param x A SNPData object, required, with donor assignments, clonotype
#'   information, and heterozygosity information for the X chromosome.
#' @param confidence_threshold Numeric, in \code{[0, 1]} (default 0.95).
#'   Posterior probability threshold for hard assignment: clonotypes whose
#'   posterior probability that a given X is the active one reaches
#'   \code{confidence_threshold} are assigned that X; clonotypes where
#'   neither X reaches the threshold receive \code{NA}.
#'
#' @return SNPData object with an additional \code{active_x} column in
#'   barcode metadata, with values "X1" or "X2" indicating the inferred
#'   active X chromosome state. Cells from clonotypes that do not meet the
#'   confidence threshold receive \code{NA}. The full XCI diagnostics are also
#'   written into the object's metadata slots,
#'   so the result can be passed directly to
#'   \code{\link{plot_xci_heatmap}},
#'   \code{\link{xci_assignments}}, and \code{\link{xci_haplotypes}}.
#'
#' @family X-chromosome inactivation functions
#' @export
#'
#' @examples
#' \dontrun{
#' # Assign the active X chromosome to cells by clonotype
#' snp_data <- assign_xci_by_clonotype(snp_data)
#'
#' # View results
#' barcode_info(snp_data) %>%
#'   count(donor, clonotype, active_x)
#'
#' # Diagnostics are stored, so plotting works directly
#' plot_xci_heatmap(snp_data, donor = "donor1")
#' }
setGeneric(
    "assign_xci_by_clonotype",
    function(
        x,
        n_inits = 10,
        confidence_threshold = 0.95
    ) {
        standardGeneric("assign_xci_by_clonotype")
    }
)

#' @rdname assign_xci_by_clonotype
#' @include SNPData-class.R
setMethod(
    "assign_xci_by_clonotype",
    signature(x = "SNPData"),
    function(x, n_inits = 10, confidence_threshold = 0.95) {
        .fit_xci(
            x,
            n_inits = n_inits,
            confidence_threshold = confidence_threshold,
            by = "clonotype"
        )
    }
)

#' Fit X-chromosome inactivation model (internal engine)
#'
#' Runs the beta-binomial EM algorithm for each donor and returns a full fit
#' object containing per-unit assignments, per-gene haplotypes, and the final
#' allele count matrices used in the model. The modelling unit is either the
#' cell (\code{by = "cell"}) or the clonotype (\code{by = "clonotype"}).
#'
#' This is the shared engine behind \code{\link{assign_xci}} and
#' \code{\link{assign_xci_by_clonotype}}. It is not exported; the public
#' entry points return a SNPData object carrying the diagnostics in its metadata
#' slots.
#'
#' @param x SNPData object containing X chromosome SNP data with donor
#'   assignments and heterozygosity information. For \code{by = "clonotype"},
#'   clonotype information must also be present.
#' @param n_inits Number of random initialisations for the EM algorithm.
#' @param confidence_threshold Posterior probability threshold for hard
#'   assignment.
#' @param by Modelling unit, \code{"cell"} or \code{"clonotype"}.
#'
#' @return The input SNPData object with diagnostics written into its metadata
#'   slots.
#'
#' @keywords internal
#' @include SNPData-class.R
.fit_xci <- function(
    x,
    n_inits = 10,
    confidence_threshold = 0.95,
    by = c("cell", "clonotype")
) {
    by <- match.arg(by)
    if (by == "clonotype") {
        .check_clonotype_available(x)
    }

    unique_donors <- barcode_info(x)$donor %>%
        unique() %>%
        setdiff(c("doublet", "unassigned")) %>%
        sort()

    # Split the data per donor up front so each worker only receives its own
    # subset rather than the full SNPData object, keeping serialised globals small.
    donor_data <- purrr::map(unique_donors, function(d) filter_samples(x, donor == d))

    result <- furrr::future_map2(
        donor_data,
        unique_donors,
        function(dd, d) {
            tryCatch(
                .fit_xci_donor(
                    dd,
                    n_inits,
                    confidence_threshold,
                    by = by
                ),
                error = function(e) {
                    logger::log_warn("Failed to fit XCI for donor {d}: {conditionMessage(e)}")
                    NULL
                }
            )
        },
        .options = furrr::furrr_options(packages = "snplet", seed = TRUE)
    ) %>%
        magrittr::set_names(unique_donors)

    fit <- structure(result, class = "xci_fit")
    .store_xci_fit(x, fit)
}

#' @keywords internal
.check_clonotype_available <- function(x) {
    barcode_info <- barcode_info(x)
    if (!"clonotype" %in% colnames(barcode_info)) {
        stop(
            "Clonotype information not available. Add clonotype data using add_barcode_metadata() or import_cellsnp() with vdj_file parameter."
        )
    }
    if (all(is.na(barcode_info$clonotype))) {
        stop(
            "All clonotype values are NA. Cannot perform clonotype-level X-inactivation assignment. Add clonotype data using add_barcode_metadata() or import_cellsnp() with vdj_file parameter."
        )
    }
    invisible(TRUE)
}

#' Fit the XCI model for one donor, at cell or clonotype level
#'
#' Shared per-donor fitter for both modelling units. For \code{by = "cell"} the
#' count matrices are used as-is; for \code{by = "clonotype"} they are first
#' aggregated by clonotype and the resulting clonotype-level posteriors are
#' projected back down to individual cells.
#'
#' @keywords internal
.fit_xci_donor <- function(
    snp_data,
    n_inits = 10,
    confidence_threshold = 0.95,
    by = c("cell", "clonotype")
) {
    by <- match.arg(by)
    donor <- unique(barcode_info(snp_data)$donor)

    if (by == "clonotype") {
        logger::log_info("[{donor}] Fitting XCI model at clonotype level")
    } else {
        logger::log_info("[{donor}] Fitting XCI model")
    }

    snp_data <- .filter_to_informative_het_snps(snp_data, donor)

    unit_data <- .xci_unit_matrices(snp_data, by, donor)
    ref_mat <- unit_data$ref_mat
    alt_mat <- unit_data$alt_mat
    unit_ids <- colnames(ref_mat)
    snp_info <- snp_info(snp_data)

    if (nrow(ref_mat) == 0 || ncol(ref_mat) == 0) {
        logger::log_warn("Insufficient data for donor {donor}")
        return(NULL)
    }

    xci_result <- .infer_xci(
        ref_mat,
        alt_mat,
        n_inits = n_inits,
        confidence_threshold = confidence_threshold,
        donor = donor
    )

    fit <- .assemble_xci_fit(
        xci_result = xci_result,
        donor = donor,
        ref_mat = ref_mat,
        alt_mat = alt_mat,
        snp_info = snp_info,
        unit_ids = unit_ids,
        unit = by
    )

    .project_to_cells(fit, unit_data$cell_to_clonotype, by)
}

#' Build the count matrices for the requested modelling unit
#'
#' For \code{by = "cell"} the per-cell count matrices are used as-is. For
#' \code{by = "clonotype"} cells with no clonotype are dropped and the
#' remaining cells' counts are summed within each clonotype, so the EM sees
#' one column per clonotype instead of per cell.
#'
#' @keywords internal
.xci_unit_matrices <- function(snp_data, by, donor = NULL) {
    if (by != "clonotype") {
        return(list(ref_mat = ref_count(snp_data), alt_mat = alt_count(snp_data), cell_to_clonotype = NULL))
    }

    has_clonotype <- !is.na(barcode_info(snp_data)$clonotype)
    if (!any(has_clonotype)) {
        stop(glue::glue("No cells with non-NA clonotype values for donor {donor}"))
    }
    snp_data <- snp_data[, has_clonotype]

    cell_to_clonotype <- barcode_info(snp_data) %>%
        dplyr::select(cell_id, clonotype)
    list(
        ref_mat = groupedRowSums(ref_count(snp_data), cell_to_clonotype$clonotype),
        alt_mat = groupedRowSums(alt_count(snp_data), cell_to_clonotype$clonotype),
        cell_to_clonotype = cell_to_clonotype
    )
}

#' Project a clonotype-level fit's posteriors down to individual cells
#'
#' A cell-level fit (\code{by = "cell"}) already keys \code{fit$assignments} on
#' \code{cell_id}, so it is returned unchanged. A clonotype-level fit keys on
#' the clonotype id instead; this adds \code{fit$cell_assignments}, one row per
#' cell, by having each cell inherit its clonotype's posterior.
#'
#' @keywords internal
.project_to_cells <- function(fit, cell_to_clonotype, by) {
    if (by != "clonotype") {
        return(fit)
    }
    fit$cell_assignments <- fit$assignments %>%
        dplyr::rename(clonotype = cell_id) %>%
        dplyr::inner_join(cell_to_clonotype, by = "clonotype") %>%
        dplyr::select(cell_id, post_X1_active, post_X2_active, assignment)
    fit
}

#' Assemble a per-donor XCI fit from an EM result
#'
#' Shared post-EM assembly used by both the cell- and clonotype-level donor
#' fitters. \code{unit_ids} label the columns of the count matrices (cell IDs
#' or clonotype IDs); the resulting \code{assignments} tibble keys on
#' \code{cell_id} regardless so accessors and the heatmap treat the modelling
#' unit uniformly.
#'
#' @keywords internal
.assemble_xci_fit <- function(xci_result, donor, ref_mat, alt_mat, snp_info, unit_ids, unit = c("cell", "clonotype")) {
    unit <- match.arg(unit)
    unit_label <- if (unit == "clonotype") "clonotypes" else "cells"

    assignments <- xci_result$post %>%
        dplyr::mutate(cell_id = unit_ids[cell], post_X2_active = 1 - post_X1_active) %>%
        dplyr::select(cell_id, post_X1_active, post_X2_active, assignment)

    # Phase and escape fraction for every gene, not just the informative subset that drove
    # calling: a gene with LLR <= 0 (inconsistent with the fitted phase) or an outlier escape
    # fraction carries no information about which X is active (correctly excluded from that
    # job) — genes with real but moderate escape are not excluded — but once active-X is
    # called, any gene can still be phased against it.
    all_genes <- .rephase_all_genes(ref_mat, alt_mat, post = xci_result$post, rho = xci_result$rho)
    haplotypes <- tibble::tibble(
        snp_id = snp_info$snp_id,
        gene_name = snp_info$gene_name,
        xci_informative = xci_result$gene_keep,
        allele_on_x1 = dplyr::case_when(
            all_genes$h_g == 0 ~ "REF",
            all_genes$h_g == 1 ~ "ALT",
            TRUE ~ NA_character_
        ),
        escape_fraction = all_genes$pi_g
    ) %>%
        dplyr::filter(!is.na(allele_on_x1))

    # Empirical background escape/noise level for this donor: the typical residual
    # escape fraction among genes that actually drove active-X calling (informative,
    # non-outlier). Recorded alongside rho so callers can use a donor-specific null
    # (e.g. in test_escape()) instead of an arbitrary fixed value.
    median_pi_g <- if (any(haplotypes$xci_informative, na.rm = TRUE)) {
        stats::median(haplotypes$escape_fraction[haplotypes$xci_informative], na.rm = TRUE)
    } else {
        NA_real_
    }

    # assignment already names the active X; "unassigned" carries no active call.
    active <- dplyr::na_if(assignments$assignment, "unassigned")
    counts <- table(factor(active, c("X1", "X2")))
    logger::log_info(
        "[{donor}] XCI fit complete: {sum(xci_result$gene_keep)} informative / {nrow(haplotypes)} phased genes, ",
        "{nrow(assignments)} {unit_label} (active X1={counts[['X1']]}, active X2={counts[['X2']]}, ",
        "unassigned={sum(is.na(active))})"
    )

    snp_info_filtered <- snp_info[xci_result$gene_keep, ]
    ref_mat_filtered <- ref_mat[xci_result$gene_keep, , drop = FALSE]
    alt_mat_filtered <- alt_mat[xci_result$gene_keep, , drop = FALSE]
    rownames(ref_mat_filtered) <- snp_info_filtered$snp_id
    rownames(alt_mat_filtered) <- snp_info_filtered$snp_id
    colnames(ref_mat_filtered) <- unit_ids
    colnames(alt_mat_filtered) <- unit_ids

    list(
        donor = donor,
        unit = unit,
        assignments = assignments,
        haplotypes = haplotypes,
        ref_mat = ref_mat_filtered,
        alt_mat = alt_mat_filtered,
        rho = xci_result$rho,
        median_pi_g = median_pi_g
    )
}

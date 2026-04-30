#' Assign inactive X chromosome to cells
#'
#' Identifies which X chromosome is inactive in female cells based on allelic
#' imbalance at heterozygous SNPs using an Expectation-Maximisation (EM)
#' algorithm with a beta-binomial likelihood.
#'
#' @details
#' X-chromosome inactivation (XCI) is a dosage compensation mechanism in female
#' mammals where one of the two X chromosomes is randomly silenced in each cell.
#' This function infers which X is inactive by fitting an EM model to allelic
#' read counts at heterozygous SNPs on the X chromosome.
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
#' the \code{inactive_x} column.
#'
#' @param x SNPData object containing X chromosome SNP data with donor
#'   assignments and heterozygosity information
#' @param n_inits Number of random initialisations for the EM algorithm.
#'   The run with the highest log-likelihood is returned. Default 10.
#' @param confidence_threshold Posterior probability threshold for hard
#'   assignment. Cells with \code{post_X1 >= confidence_threshold} are
#'   assigned "X1", those with \code{post_X1 <= 1 - confidence_threshold}
#'   are assigned "X2", and the remainder receive \code{NA}. Default 0.95.
#' @param refit_after_filter Logical; if TRUE, re-run the EM algorithm after
#'   filtering genes with inconsistent allelic patterns. Provides sharper
#'   posteriors on the cleaned gene set. Default FALSE.
#'
#' @return SNPData object with an additional \code{inactive_x} column in
#'   barcode metadata, with values "X1" or "X2" indicating the inferred
#'   inactive X chromosome state. Cells that do not meet the confidence
#'   threshold receive \code{NA}.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Assign inactive X chromosome to cells
#' snp_data <- assign_inactive_x(snp_data)
#'
#' # View results
#' get_barcode_info(snp_data) %>%
#'   count(donor, inactive_x)
#' }
setGeneric("assign_inactive_x", function(x, n_inits = 10, confidence_threshold = 0.95, refit_after_filter = FALSE) {
    standardGeneric("assign_inactive_x")
})

#' @rdname assign_inactive_x
#' @include SNPData-class.R
setMethod(
    "assign_inactive_x",
    signature(x = "SNPData"),
    function(x, n_inits = 10, confidence_threshold = 0.95, refit_after_filter = FALSE) {
        unique_donors <- sort(unique(get_barcode_info(x)$donor))
        result <- furrr::future_map(
            unique_donors,
            function(d) {
                tryCatch(
                    .assign_inactive_x_single_donor(
                        filter_samples(x, donor == d),
                        n_inits,
                        confidence_threshold,
                        refit_after_filter
                    ),
                    error = function(e) {
                        logger::log_warn("Failed to assign inactive X for donor {d}: {conditionMessage(e)}")
                        tibble::tibble(cell_id = character(), inactive_x = character())
                    }
                )
            },
            .options = furrr::furrr_options(packages = "snplet", seed = TRUE)
        ) %>%
            dplyr::bind_rows()
        add_barcode_metadata(x, result)
    }
)

#' Assign inactive X chromosome to cells by clonotype
#'
#' Identifies which X chromosome is inactive in female cells based on allelic
#' imbalance at heterozygous SNPs, aggregating counts to the clonotype level
#' before running the EM model and projecting assignments back to individual cells.
#'
#' @details
#' This function is similar to \code{\link{assign_inactive_x}} but runs the
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
#' @param x SNPData object containing X chromosome SNP data with donor
#'   assignments, clonotype information, and heterozygosity information
#' @param n_inits Number of random initialisations for the EM algorithm.
#'   The run with the highest log-likelihood is returned. Default 10.
#' @param confidence_threshold Posterior probability threshold for hard
#'   assignment. Clonotypes with \code{post_X1 >= confidence_threshold} are
#'   assigned "X1", those with \code{post_X1 <= 1 - confidence_threshold}
#'   are assigned "X2", and the remainder receive \code{NA}. Default 0.95.
#' @param refit_after_filter Logical; if TRUE, re-run the EM algorithm after
#'   filtering genes with inconsistent allelic patterns. Provides sharper
#'   posteriors on the cleaned gene set. Default FALSE.
#'
#' @return SNPData object with an additional \code{inactive_x} column in
#'   barcode metadata, with values "X1" or "X2" indicating the inferred
#'   inactive X chromosome state. Cells from clonotypes that do not meet the
#'   confidence threshold receive \code{NA}.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Assign inactive X chromosome to cells by clonotype
#' snp_data <- assign_inactive_x_by_clonotype(snp_data)
#'
#' # View results
#' get_barcode_info(snp_data) %>%
#'   count(donor, clonotype, inactive_x)
#' }
setGeneric(
    "assign_inactive_x_by_clonotype",
    function(x, n_inits = 10, confidence_threshold = 0.95, refit_after_filter = FALSE) {
        standardGeneric("assign_inactive_x_by_clonotype")
    }
)

#' @rdname assign_inactive_x_by_clonotype
#' @include SNPData-class.R
setMethod(
    "assign_inactive_x_by_clonotype",
    signature(x = "SNPData"),
    function(x, n_inits = 10, confidence_threshold = 0.95, refit_after_filter = FALSE) {
        barcode_info <- get_barcode_info(x)
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
        unique_donors <- sort(unique(barcode_info$donor))
        result <- furrr::future_map(
            unique_donors,
            function(d) {
                tryCatch(
                    .assign_inactive_x_single_donor_by_clonotype(
                        filter_samples(x, donor == d),
                        n_inits,
                        confidence_threshold,
                        refit_after_filter
                    ),
                    error = function(e) {
                        logger::log_warn("Failed to assign inactive X for donor {d}: {conditionMessage(e)}")
                        tibble::tibble(cell_id = character(), inactive_x = character())
                    }
                )
            },
            .options = furrr::furrr_options(packages = "snplet", seed = TRUE)
        ) %>%
            dplyr::bind_rows()
        add_barcode_metadata(x, result)
    }
)

#' Fit X-chromosome inactivation model
#'
#' Runs the beta-binomial EM algorithm for each donor and returns a full fit
#' object containing per-cell assignments, per-gene haplotypes, and the final
#' allele count matrices used in the model. Use this when diagnostic outputs
#' are needed; for barcode annotation only, use \code{\link{assign_inactive_x}}.
#'
#' @param x SNPData object containing X chromosome SNP data with donor
#'   assignments and heterozygosity information
#' @param n_inits Number of random initialisations for the EM algorithm.
#'   Default 10.
#' @param confidence_threshold Posterior probability threshold for hard
#'   assignment. Default 0.95.
#' @param refit_after_filter Logical; if TRUE, re-run the EM algorithm after
#'   filtering genes with inconsistent allelic patterns. Provides sharper
#'   posteriors on the cleaned gene set. Default FALSE.
#'
#' @return An \code{xci_fit} object (named list by donor). Use
#'   \code{\link{xci_assignments}}, \code{\link{xci_haplotypes}}, and
#'   \code{\link{plot_inactive_x_assignment_heatmap}} to extract results.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' fit <- fit_inactive_x(snp_data)
#' xci_assignments(fit)
#' xci_haplotypes(fit)
#' plot_inactive_x_assignment_heatmap(fit, donor = "donor1")
#' }
setGeneric("fit_inactive_x", function(x, n_inits = 10, confidence_threshold = 0.95, refit_after_filter = FALSE) {
    standardGeneric("fit_inactive_x")
})

#' @rdname fit_inactive_x
#' @include SNPData-class.R
setMethod(
    "fit_inactive_x",
    signature(x = "SNPData"),
    function(x, n_inits = 10, confidence_threshold = 0.95, refit_after_filter = FALSE) {
        unique_donors <- sort(unique(get_barcode_info(x)$donor))

        result <- furrr::future_map(
            unique_donors,
            function(d) {
                tryCatch(
                    .fit_xci_donor(filter_samples(x, donor == d), n_inits, confidence_threshold, refit_after_filter),
                    error = function(e) {
                        logger::log_warn("Failed to fit XCI for donor {d}: {conditionMessage(e)}")
                        NULL
                    }
                )
            },
            .options = furrr::furrr_options(packages = "snplet", seed = TRUE)
        ) %>%
            magrittr::set_names(unique_donors)

        structure(result, class = "xci_fit")
    }
)

#' Extract cell assignments from an XCI fit
#'
#' @param fit An \code{xci_fit} object from \code{\link{fit_inactive_x}}
#'
#' @return A tibble with columns \code{donor}, \code{cell_id}, \code{post_X1},
#'   \code{post_X2}, and \code{assignment}.
#'
#' @export
xci_assignments <- function(fit) {
    stopifnot(inherits(fit, "xci_fit"))
    purrr::keep(fit, ~ !is.null(.x)) %>%
        purrr::map(~ .x$assignments) %>%
        dplyr::bind_rows(.id = "donor")
}

#' Extract SNP haplotypes from an XCI fit
#'
#' Returns the inferred phase and escape fraction for each SNP in the final
#' gene set used by the EM model.
#'
#' @param fit An \code{xci_fit} object from \code{\link{fit_inactive_x}}
#'
#' @return A tibble with columns \code{donor}, \code{snp_id}, \code{gene_name},
#'   \code{allele_on_x1} ("REF" or "ALT"), and \code{pi_g} (estimated escape
#'   fraction from the inactive chromosome).
#'
#' @export
xci_haplotypes <- function(fit) {
    stopifnot(inherits(fit, "xci_fit"))
    purrr::keep(fit, ~ !is.null(.x)) %>%
        purrr::map(~ .x$haplotypes) %>%
        dplyr::bind_rows(.id = "donor")
}

#' Plot assignment heatmap from an XCI fit
#'
#' Visualizes the REF allele fraction at the final EM-selected genes for a
#' single donor. Cells are ordered by assignment (X1 then X2 then unassigned)
#' and annotated at the top. Rows are the genes that survived both the outlier
#' filter and the post-convergence escapee filter.
#'
#' @param fit An \code{xci_fit} object from \code{\link{fit_inactive_x}}
#' @param donor Character string specifying which donor to visualize
#'
#' @return A ComplexHeatmap object
#'
#' @export
plot_inactive_x_assignment_heatmap <- function(fit, donor) {
    stopifnot(inherits(fit, "xci_fit"))
    donor_fit <- fit[[donor]]
    if (is.null(donor_fit)) {
        stop(paste0("No fit available for donor: ", donor))
    }

    ref_mat <- donor_fit$ref_mat
    alt_mat <- donor_fit$alt_mat
    assignments <- donor_fit$assignments

    # REF allele fraction; NA where uncovered
    cov_mat <- ref_mat + alt_mat
    frac_mat <- ref_mat / cov_mat
    frac_mat[cov_mat == 0] <- NA

    # Drop rows that are entirely NA — they carry no signal and break hclust
    row_has_data <- rowSums(!is.na(frac_mat)) > 0
    frac_mat <- frac_mat[row_has_data, , drop = FALSE]

    # Order cells: X1 → X2 → unassigned
    idx <- order(assignments$post_X1, decreasing = TRUE)
    frac_mat <- frac_mat[, assignments$cell_id[idx], drop = FALSE]

    col_ann <- ComplexHeatmap::HeatmapAnnotation(
        assignment = assignments$assignment[idx],
        posterior_X1 = assignments$post_X1[idx],
        col = list(
            assignment = c(X1 = "#4e79a7", X2 = "#f28e2b", unassigned = "grey70"),
            posterior_X1 = circlize::colorRamp2(c(0, 0.1, 0.9, 1), c("#e9a3c9", "white", "white", "#a1d76a"))
        ),
        show_annotation_name = TRUE
    )

    # Impute NAs with 0.5 (balanced) only for row clustering so that sparse
    # coverage does not produce NaN distances; the display matrix keeps NAs.
    frac_for_cluster <- frac_mat
    frac_for_cluster[is.na(frac_for_cluster)] <- 0.5
    row_dend <- stats::hclust(stats::dist(as.matrix(frac_for_cluster)), method = "ward.D2")

    ComplexHeatmap::Heatmap(
        as.matrix(frac_mat),
        top_annotation = col_ann,
        cluster_columns = FALSE,
        cluster_rows = row_dend,
        show_row_names = FALSE,
        show_column_names = FALSE,
        name = "REF fraction",
        na_col = "grey95"
    )
}

.fit_xci_donor <- function(snp_data, n_inits = 10, confidence_threshold = 0.95, refit_after_filter = FALSE) {
    donor <- unique(get_barcode_info(snp_data)$donor)
    logger::log_info("Fitting XCI model for donor {donor}")

    snp_data <- .filter_to_informative_het_snps(snp_data)

    ref_mat <- ref_count(snp_data)
    alt_mat <- alt_count(snp_data)
    cell_ids <- colnames(ref_mat)
    snp_info <- get_snp_info(snp_data)

    if (nrow(ref_mat) == 0 || ncol(ref_mat) == 0) {
        logger::log_warn("Insufficient data for donor {donor}")
        return(NULL)
    }

    xci_result <- .infer_xci(
        ref_mat,
        alt_mat,
        n_inits = n_inits,
        confidence_threshold = confidence_threshold,
        refit_after_filter = refit_after_filter
    )

    assignments <- xci_result$post %>%
        dplyr::mutate(cell_id = cell_ids[cell], post_X2 = 1 - post_X1) %>%
        dplyr::select(cell_id, post_X1, post_X2, assignment)

    snp_info_filtered <- snp_info[xci_result$gene_keep, ]
    haplotypes <- tibble::tibble(
        snp_id = snp_info_filtered$snp_id,
        gene_name = snp_info_filtered$gene_name,
        allele_on_x1 = ifelse(xci_result$h_g == 0, "REF", "ALT"),
        pi_g = xci_result$pi_g
    )

    ref_mat_filtered <- ref_mat[xci_result$gene_keep, , drop = FALSE]
    alt_mat_filtered <- alt_mat[xci_result$gene_keep, , drop = FALSE]
    rownames(ref_mat_filtered) <- snp_info_filtered$snp_id
    rownames(alt_mat_filtered) <- snp_info_filtered$snp_id
    colnames(ref_mat_filtered) <- cell_ids
    colnames(alt_mat_filtered) <- cell_ids

    list(
        donor = donor,
        assignments = assignments,
        haplotypes = haplotypes,
        ref_mat = ref_mat_filtered,
        alt_mat = alt_mat_filtered
    )
}

.assign_inactive_x_single_donor <- function(
    snp_data,
    n_inits = 10,
    confidence_threshold = 0.95,
    refit_after_filter = FALSE
) {
    fit <- .fit_xci_donor(snp_data, n_inits, confidence_threshold, refit_after_filter)
    if (is.null(fit)) {
        return(tibble::tibble(cell_id = character(), inactive_x = character()))
    }
    fit$assignments %>%
        dplyr::filter(assignment != "unassigned") %>%
        dplyr::select(cell_id, inactive_x = assignment)
}

.assign_inactive_x_single_donor_by_clonotype <- function(
    snp_data,
    n_inits = 10,
    confidence_threshold = 0.95,
    refit_after_filter = FALSE
) {
    donor <- unique(get_barcode_info(snp_data)$donor)
    logger::log_info("Assigning inactive X for donor {donor} at clonotype level")

    snp_data <- .filter_to_informative_het_snps(snp_data)

    barcode_info <- get_barcode_info(snp_data)
    has_clonotype <- !is.na(barcode_info$clonotype)
    if (!any(has_clonotype)) {
        stop(glue::glue("No cells with non-NA clonotype values for donor {donor}"))
    }
    snp_data <- snp_data[, has_clonotype]
    clonotypes <- get_barcode_info(snp_data)$clonotype

    ref_mat <- groupedRowSums(ref_count(snp_data), clonotypes)
    alt_mat <- groupedRowSums(alt_count(snp_data), clonotypes)
    clonotype_ids <- colnames(ref_mat)

    if (nrow(ref_mat) == 0 || ncol(ref_mat) == 0) {
        logger::log_warn("Skipping inactive X assignment for donor {donor}: insufficient data")
        return(tibble::tibble(cell_id = character(), inactive_x = character()))
    }

    xci_result <- .infer_xci(
        ref_mat,
        alt_mat,
        n_inits = n_inits,
        confidence_threshold = confidence_threshold,
        refit_after_filter = refit_after_filter
    )

    clonotype_assignment <- xci_result$post %>%
        dplyr::filter(assignment != "unassigned") %>%
        dplyr::mutate(clonotype = clonotype_ids[cell]) %>%
        dplyr::select(clonotype, inactive_x = assignment)

    # Project clonotype assignments back to cells
    barcode_info <- get_barcode_info(snp_data)
    dplyr::select(barcode_info, cell_id, clonotype) %>%
        dplyr::inner_join(clonotype_assignment, by = "clonotype") %>%
        dplyr::select(cell_id, inactive_x)
}

.filter_to_informative_het_snps <- function(snp_data) {
    het_snp_ids <- snp_data %>%
        donor_het_status_df() %>%
        dplyr::filter(zygosity == "het") %>%
        dplyr::pull(snp_id)

    snp_data <- snp_data %>%
        filter_snps(snp_id %in% het_snp_ids)

    top_snp_per_gene <- get_snp_info(snp_data) %>%
        dplyr::arrange(dplyr::desc(coverage)) %>%
        dplyr::slice_head(n = 1, by = "gene_name")

    snp_data %>%
        filter_snps(snp_id %in% top_snp_per_gene$snp_id)
}

.infer_xci <- function(
    ref_mat,
    alt_mat,
    n_inits = 10,
    confidence_threshold = 0.95,
    min_cells = 10,
    min_cov = 1,
    refit_after_filter = FALSE
) {
    passes_outlier_filter <- .filter_outlier_genes(ref_mat, alt_mat, min_cells, min_cov)
    ref_mat <- ref_mat[passes_outlier_filter, , drop = FALSE]
    alt_mat <- alt_mat[passes_outlier_filter, , drop = FALSE]

    dat <- .pivot_counts_to_long(ref_mat, alt_mat, min_cov)
    n_genes <- nrow(ref_mat)

    fits <- lapply(seq_len(n_inits), function(s) .run_em(dat, n_genes, init_seed = s))
    best <- fits[[which.max(sapply(fits, `[[`, "ll"))]]

    # Post-convergence escapee filter: genes with LLR <= 0 are inconsistent with
    # current cell assignments; MAD filter on pi_g removes outlier escape fractions.
    passes_escapee_filter <- .filter_escapee_genes(dat, n_genes, best$h_g, best$pi_g, best$rho, best$post)
    escaped_gene_indices <- which(passes_escapee_filter)

    if (refit_after_filter && !all(passes_escapee_filter)) {
        # Re-run EM on the cleaned gene set for sharper posteriors.
        n_removed <- sum(!passes_escapee_filter)
        logger::log_info("Removing {n_removed} escape genes after initial EM; re-running")

        gene_indices <- which(passes_escapee_filter)
        dat <- dat %>%
            dplyr::filter(gene %in% gene_indices) %>%
            dplyr::mutate(gene = match(gene, gene_indices))
        n_genes <- length(gene_indices)

        fits <- lapply(seq_len(n_inits), function(s) .run_em(dat, n_genes, init_seed = s))
        best <- fits[[which.max(sapply(fits, `[[`, "ll"))]]

        # After refit, all genes in dat passed the escapee filter, so mark them all TRUE
        passes_escapee_filter <- rep(TRUE, n_genes)
    }

    best$post <- best$post %>%
        dplyr::mutate(
            assignment = dplyr::case_when(
                post_X1 >= confidence_threshold ~ "X1",
                post_X1 <= 1 - confidence_threshold ~ "X2",
                TRUE ~ "unassigned"
            )
        )

    # gene_keep: logical of length nrow(original ref_mat), TRUE = gene survived
    # both the outlier filter and the post-convergence escapee filter
    gene_keep <- passes_outlier_filter
    if (refit_after_filter) {
        # Map filtered gene indices back to original space
        gene_keep[escaped_gene_indices] <- passes_escapee_filter
    } else {
        # Direct assignment works when no refit (passes_escapee_filter has same length as filtered genes)
        gene_keep[passes_outlier_filter] <- passes_escapee_filter
    }

    best$h_g  <- best$h_g[passes_escapee_filter]
    best$pi_g <- best$pi_g[passes_escapee_filter]

    c(best, list(gene_keep = gene_keep))
}

.compute_gene_llr <- function(dat, post, h_g, pi_g, rho) {
    p_if_X1 <- ifelse(h_g[dat$gene] == 0, pi_g[dat$gene], 1 - pi_g[dat$gene])
    p_if_X2 <- 1 - p_if_X1

    ll_X1 <- .loglik_obs(dat$ref, dat$n, p_if_X1, rho)
    ll_X2 <- .loglik_obs(dat$ref, dat$n, p_if_X2, rho)

    dat %>%
        dplyr::left_join(post %>% dplyr::select(cell, post_X1), by = "cell") %>%
        # Expected LLR per observation: weighted by P(X1-inactive) for each cell.
        # Positive contribution means the observation is consistent with current assignments.
        dplyr::mutate(llr = post_X1 * (ll_X1 - ll_X2)) %>%
        dplyr::group_by(gene) %>%
        dplyr::summarise(llr = sum(llr), .groups = "drop")
}

.filter_escapee_genes <- function(dat, n_genes, h_g, pi_g, rho, post, mad_threshold = 2) {
    gene_llr <- .compute_gene_llr(dat, post, h_g, pi_g, rho)

    # Genes absent from dat have no observations above min_cov — exclude them
    keep_llr <- rep(FALSE, n_genes)
    keep_llr[gene_llr$gene[gene_llr$llr > 0]] <- TRUE # LLR > 0: data supports current assignment

    # Secondary pass: remove genes with unusually high escape fraction relative
    # to the rest of the sample. Guard against zero MAD (all pi_g identical).
    pi_mad <- stats::mad(pi_g)
    keep_pi <- if (pi_mad > 0) {
        # Robust z-score relative to the sample's own escape distribution
        (pi_g - stats::median(pi_g)) / pi_mad <= mad_threshold
    } else {
        rep(TRUE, n_genes) # all pi_g identical — no outliers possible
    }

    keep_llr & keep_pi
}

.run_em <- function(dat, n_genes, max_iter = 50, tol = 1e-4, init_seed = 1) {
    h_g <- withr::with_seed(init_seed, sample(0:1, n_genes, replace = TRUE)) # random phase initialisation
    pi_g <- rep(0.05, n_genes) # start conservatively: assume 5% escape fraction
    rho <- 0.05 # fixed beta-binomial overdispersion

    ll_prev <- -Inf
    for (iter in seq_len(max_iter)) {
        post <- .e_step(dat, h_g, pi_g, rho)
        h_g <- .m_step_phase(dat, post, pi_g, rho, h_g)
        pi_g <- .m_step_pi(dat, post, h_g)
        # log-likelihood: sum of log(sigmoid(lor)) for each cell
        # Numerically stable: log(sigmoid(lor)) = pmin(0, lor) - log(1 + exp(-|lor|))
        ll_current <- sum(pmin(0, post$lor) - log1p(exp(-abs(post$lor))))
        if (abs(ll_current - ll_prev) < tol) {
            break
        }
        ll_prev <- ll_current
    }
    list(post = post, h_g = h_g, pi_g = pi_g, rho = rho, ll = ll_current)
}

.e_step <- function(dat, h_g, pi_g, rho, prior = 0.5) {
    # p(REF | X1 inactive): if h=0 (X1 carries REF), REF is silenced → pi_g (small)
    #                        if h=1 (X1 carries ALT), REF is active  → 1 - pi_g (large)
    p_if_X1 <- ifelse(h_g[dat$gene] == 0, pi_g[dat$gene], 1 - pi_g[dat$gene])
    p_if_X2 <- 1 - p_if_X1 # p(REF | X2 inactive) is the complement by symmetry

    ll_X1 <- .loglik_obs(dat$ref, dat$n, p_if_X1, rho)
    ll_X2 <- .loglik_obs(dat$ref, dat$n, p_if_X2, rho)

    logit_prior <- log(prior / (1 - prior)) # log(1) = 0 for equal prior
    cell_ll <- dat %>%
        dplyr::mutate(ll_X1 = ll_X1, ll_X2 = ll_X2) %>%
        dplyr::group_by(cell) %>%
        dplyr::summarise(sum_X1 = sum(ll_X1), sum_X2 = sum(ll_X2), .groups = "drop")

    # log-odds ratio: Σ_g [loglik(X1) - loglik(X2)] + logit(prior)
    lor <- cell_ll$sum_X1 - cell_ll$sum_X2 + logit_prior
    cell_ll$post_X1 <- 1 / (1 + exp(-lor)) # sigmoid converts LOR to posterior
    cell_ll$lor <- lor
    cell_ll
}

.m_step_phase <- function(dat, post, pi_g, rho, h_g) {
    counts_with_posterior <- dat %>%
        dplyr::left_join(post %>% dplyr::select(cell, post_X1), by = "cell")

    # Expected log-likelihood under each phase: E_q[log p(ref | h, pi_g)]
    # h=0: X1-inactive cells (weight post_X1) see REF fraction pi_g (silenced);
    #       X2-inactive cells (weight 1-post_X1) see REF fraction 1 - pi_g (active)
    # h=1: REF and ALT roles are swapped
    ll_by_gene <- counts_with_posterior %>%
        dplyr::group_by(gene) %>%
        dplyr::summarise(
            ll_h0 = sum(
                post_X1 *
                    VGAM::dbetabinom(ref, n, pi_g[gene[1]], rho, log = TRUE) +
                    (1 - post_X1) * VGAM::dbetabinom(ref, n, 1 - pi_g[gene[1]], rho, log = TRUE)
            ),
            ll_h1 = sum(
                post_X1 *
                    VGAM::dbetabinom(ref, n, 1 - pi_g[gene[1]], rho, log = TRUE) +
                    (1 - post_X1) * VGAM::dbetabinom(ref, n, pi_g[gene[1]], rho, log = TRUE)
            ),
            .groups = "drop"
        )

    h_g_new <- h_g
    h_g_new[ll_by_gene$gene] <- as.integer(ll_by_gene$ll_h1 > ll_by_gene$ll_h0)
    h_g_new
}

.m_step_pi <- function(dat, post, h_g, pi_bounds = c(0.001, 0.499)) {
    counts_with_posterior <- dat %>%
        dplyr::left_join(post %>% dplyr::select(cell, post_X1), by = "cell") %>%
        dplyr::mutate(
            # Soft-assigned inactive-allele read count:
            # h=0 (X1 carries REF): X1-inactive cells contribute REF, X2-inactive contribute ALT
            # h=1 (X1 carries ALT): roles reversed
            xi_ref_count = ifelse(
                h_g[gene] == 0,
                post_X1 * ref + (1 - post_X1) * alt,
                post_X1 * alt + (1 - post_X1) * ref
            ),
            xi_total = n
        )

    # MLE: inactive-allele fraction = (soft inactive reads) / (total reads)
    pi_new <- counts_with_posterior %>%
        dplyr::group_by(gene) %>%
        dplyr::summarise(pi = sum(xi_ref_count) / sum(xi_total), .groups = "drop")

    pi_g_new <- rep(0.05, length(h_g))
    # Clamp to (0.001, 0.499) so pi_g stays interpretable as a minor fraction
    pi_g_new[pi_new$gene] <- pmax(pi_bounds[1], pmin(pi_bounds[2], pi_new$pi))
    pi_g_new
}

.loglik_obs <- function(ref, n, p, rho) {
    # Beta-binomial: rho > 0 adds overdispersion relative to binomial; rho = 0 reduces to binomial
    VGAM::dbetabinom(ref, size = n, prob = p, rho = rho, log = TRUE)
}

.filter_outlier_genes <- function(ref_mat, alt_mat, min_cells = 10, min_cov = 1, mad_threshold = 2) {
    # Coverage per cell-gene pair
    n_mat <- ref_mat + alt_mat
    covered <- n_mat >= min_cov

    # Per gene: how many cells have sufficient coverage, and how many of those are REF-majority
    n_expressing <- rowSums(covered)
    ref_majority <- rowSums(covered & (ref_mat > alt_mat))

    # Drop genes seen in too few cells — not enough information to estimate allelic skew
    passes_count_filter <- n_expressing >= min_cells

    # Allelic skew: fraction of covered cells where REF > ALT, folded to [0.5, 1]
    # so that genes skewed toward either allele score equally high
    skew <- ref_majority[passes_count_filter] / n_expressing[passes_count_filter]
    skew <- pmax(skew, 1 - skew)

    # Robust z-score: genes with unusually extreme skew relative to the rest are
    # likely systematic (e.g. mapping bias, escape from XCI) rather than informative.
    # Guard against zero MAD (all skew values identical).
    skew_mad <- stats::mad(skew)
    passes_skew_filter <- if (skew_mad > 0) {
        z <- (skew - stats::median(skew)) / skew_mad
        abs(z) <= mad_threshold
    } else {
        rep(TRUE, length(skew)) # all skew identical — no outliers possible
    }

    keep <- rep(FALSE, nrow(ref_mat))
    keep[passes_count_filter] <- passes_skew_filter
    keep
}

.pivot_counts_to_long <- function(ref_mat, alt_mat, min_cov = 3) {
    n_mat <- ref_mat + alt_mat
    # Matrix::which handles sparse lgCMatrix; base which() does not dispatch S4
    idx <- Matrix::which(n_mat >= min_cov, arr.ind = TRUE)
    tibble::tibble(
        gene = idx[, 1],
        cell = idx[, 2],
        ref = ref_mat[idx],
        alt = alt_mat[idx],
        n = n_mat[idx]
    )
}

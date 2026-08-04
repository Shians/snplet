#' @keywords internal
.filter_to_informative_het_snps <- function(snp_data, donor = NULL, zygosity_source = NULL) {
    # Suppress the generic per-filter logs from filter_snps; we emit a single
    # consolidated, donor-labelled summary instead.
    old_threshold <- logger::log_threshold()
    logger::log_threshold(logger::WARN)
    on.exit(logger::log_threshold(old_threshold), add = TRUE)

    n_start <- nrow(get_snp_info(snp_data))

    # Restrict to X-chromosome SNPs. XCI is an X-chromosome phenomenon, so SNPs
    # on other chromosomes carry no signal and would only add noise. Use the
    # canonical (UCSC) chromosome name set at construction so the match is robust
    # to the input naming style.
    snp_info <- get_snp_info(snp_data)
    if (!"chrom_canonical" %in% colnames(snp_info)) {
        stop(
            "No canonical chromosome names available. Ensure the SNPData object was built with a 'chrom' column so chrX SNPs can be selected."
        )
    }
    snp_data <- snp_data %>%
        filter_snps(chrom_canonical == "chrX")
    n_chrx <- nrow(get_snp_info(snp_data))

    het_snp_ids <- snp_data %>%
        donor_het_status_df(zygosity_source = zygosity_source) %>%
        dplyr::filter(zygosity == "het") %>%
        dplyr::pull(snp_id)

    snp_data <- snp_data %>%
        filter_snps(snp_id %in% het_snp_ids)
    n_het <- nrow(get_snp_info(snp_data))

    top_snp_per_gene <- get_snp_info(snp_data) %>%
        dplyr::arrange(dplyr::desc(coverage)) %>%
        dplyr::slice_head(n = 1, by = "gene_name")

    snp_data <- snp_data %>%
        filter_snps(snp_id %in% top_snp_per_gene$snp_id)
    n_genes <- nrow(get_snp_info(snp_data))

    logger::log_info(
        "[{donor}] het-SNP filter: {n_start} SNPs -> {n_chrx} chrX -> {n_het} het -> {n_genes} genes (top SNP per gene)"
    )

    snp_data
}

#' @keywords internal
.infer_xci <- function(
    ref_mat,
    alt_mat,
    n_inits = 10,
    confidence_threshold = 0.95,
    min_cells = 10,
    min_cov = 1,
    donor = NULL
) {
    passes_outlier_filter <- .filter_outlier_genes(ref_mat, alt_mat, min_cells, min_cov)
    ref_mat <- ref_mat[passes_outlier_filter, , drop = FALSE]
    alt_mat <- alt_mat[passes_outlier_filter, , drop = FALSE]

    dat <- .pivot_counts_to_long(ref_mat, alt_mat, min_cov)
    n_genes <- nrow(ref_mat)

    logger::log_info(
        "[{donor}] Running EM: {n_inits} random restarts over {n_genes} genes, {nrow(dat)} observations"
    )
    fits <- purrr::map(seq_len(n_inits), function(s) {
        fit <- .run_em(dat, n_genes, init_seed = s)
        logger::log_debug("[{donor}] EM restart {s}/{n_inits} done (logLik = {round(fit$ll, 2)})")
        fit
    })
    best <- fits[[which.max(purrr::map_dbl(fits, "ll"))]]
    logger::log_info("[{donor}] EM complete: best logLik = {round(best$ll, 2)}")

    # Post-convergence escapee filter: genes with LLR <= 0 are inconsistent with
    # current cell assignments; MAD filter on pi_g removes outlier escape fractions.
    # This decides which genes get reported as xci_informative (drove active-X
    # calling) — it does not change best$post, so escapee genes' contribution to
    # the cell-assignment posterior above is already baked in either way. Computed
    # against the pre-relabel h_g/post (relabelling below is a symmetric X1<->X2
    # swap and does not change which genes pass).
    passes_escapee_filter <- .filter_escapee_genes(dat, n_genes, best$h_g, best$pi_g, best$rho, best$post)

    best$post <- best$post %>%
        dplyr::mutate(
            assignment = dplyr::case_when(
                post_X1_inactive >= confidence_threshold ~ "X1",
                post_X1_inactive <= 1 - confidence_threshold ~ "X2",
                TRUE ~ "unassigned"
            )
        )

    best <- .canonicalise_x1_label(best)

    # gene_keep: logical of length nrow(original ref_mat), TRUE = gene survived
    # both the outlier filter and the post-convergence escapee filter.
    gene_keep <- passes_outlier_filter
    gene_keep[passes_outlier_filter] <- passes_escapee_filter

    c(best, list(gene_keep = gene_keep))
}

#' Canonicalise which mixture component is called "X1"
#'
#' The EM's component naming is arbitrary (fixed only by the random phase
#' initialisation), so the same biological haplotype could end up "X1" in one
#' run and "X2" in the next. Relabels so X1 always denotes the active-X
#' majority (assignment "X2" == X1 active, per \code{\link{.active_from_inactive}}).
#'
#' @keywords internal
.canonicalise_x1_label <- function(fit) {
    if (sum(fit$post$assignment == "X1") > sum(fit$post$assignment == "X2")) {
        fit$post <- fit$post %>%
            dplyr::mutate(
                post_X1_inactive = 1 - post_X1_inactive,
                assignment = dplyr::case_when(
                    assignment == "X1" ~ "X2",
                    assignment == "X2" ~ "X1",
                    TRUE ~ assignment
                )
            )
        fit$h_g <- 1L - fit$h_g
    }
    fit
}

#' @keywords internal
.compute_gene_llr <- function(dat, post, h_g, pi_g, rho) {
    p_if_X1 <- ifelse(h_g[dat$gene] == 0, pi_g[dat$gene], 1 - pi_g[dat$gene])
    p_if_X2 <- 1 - p_if_X1

    ll_X1 <- .loglik_obs(dat$ref, dat$n, p_if_X1, rho)
    ll_X2 <- .loglik_obs(dat$ref, dat$n, p_if_X2, rho)

    dat %>%
        dplyr::left_join(post %>% dplyr::select(cell, post_X1_inactive), by = "cell") %>%
        # Expected LLR per observation: weighted by P(X1-inactive) for each cell.
        # Positive contribution means the observation is consistent with current assignments.
        dplyr::mutate(llr = post_X1_inactive * (ll_X1 - ll_X2)) %>%
        dplyr::group_by(gene) %>%
        dplyr::summarise(llr = sum(llr), .groups = "drop")
}

#' @keywords internal
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

#' @keywords internal
.run_em <- function(dat, n_genes, max_iter = 50, tol = 1e-4, init_seed = 1) {
    h_g <- withr::with_seed(init_seed, sample(0:1, n_genes, replace = TRUE)) # random phase initialisation
    pi_g <- rep(0.05, n_genes) # start conservatively: assume 5% escape fraction
    rho <- 0.05 # fixed beta-binomial overdispersion

    # Unique kernel arguments are fixed across iterations; index them once.
    dedup <- .build_ll_dedup(dat)

    ll_prev <- -Inf
    for (iter in seq_len(max_iter)) {
        # Both escape orientations, computed once and shared by E- and M-steps.
        ll <- .betabinom_ll_both(dedup, pi_g, rho)
        post <- .e_step(dat, h_g, ll)
        h_g <- .m_step_phase(dat, post, h_g, ll)
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

#' @keywords internal
.e_step <- function(dat, h_g, ll, prior = 0.5) {
    # ll$L0 = loglik at pi_g[gene] (REF silenced), ll$L1 = loglik at 1 - pi_g[gene].
    # p(REF | X1 inactive): if h=0 (X1 carries REF), REF is silenced → pi_g → L0
    #                        if h=1 (X1 carries ALT), REF is active  → 1 - pi_g → L1
    h_row <- h_g[dat$gene] == 0
    # Per-observation X1-vs-X2 log-likelihood contrast. X2 is the complement
    # orientation, so the contrast is +(L0 - L1) when h=0 and -(L0 - L1) when h=1.
    obs_lor <- ifelse(h_row, ll$L0 - ll$L1, ll$L1 - ll$L0)

    logit_prior <- log(prior / (1 - prior)) # log(1) = 0 for equal prior
    # Σ_g [loglik(X1) - loglik(X2)] per cell, via rowsum (faster than group_by).
    cell_lor <- rowsum(obs_lor, dat$cell)
    lor <- as.numeric(cell_lor) + logit_prior

    tibble::tibble(
        cell = as.integer(rownames(cell_lor)),
        post_X1_inactive = 1 / (1 + exp(-lor)), # sigmoid converts LOR to posterior
        lor = lor
    )
}

#' @keywords internal
.m_step_phase <- function(dat, post, h_g, ll) {
    # Expected log-likelihood under each phase: E_q[log p(ref | h, pi_g)], built
    # from the shared orientation pair. ll$L0 = loglik at pi_g (REF silenced),
    # ll$L1 = loglik at 1 - pi_g (REF active).
    # h=0: X1-inactive cells (weight post_X1_inactive) see REF fraction pi_g → L0;
    #       X2-inactive cells (weight 1-post_X1_inactive) see 1 - pi_g → L1.
    # h=1: the two orientations swap.
    # Scatter cell posteriors into a lookup indexed by cell id. post$cell is a
    # sorted subset of cell ids (cells with no covered gene are absent), so a
    # positional index would misalign — index by id instead.
    post_lookup <- numeric(max(dat$cell))
    post_lookup[post$cell] <- post$post_X1_inactive
    post_X1_inactive <- post_lookup[dat$cell]

    obs_h0 <- post_X1_inactive * ll$L0 + (1 - post_X1_inactive) * ll$L1
    obs_h1 <- post_X1_inactive * ll$L1 + (1 - post_X1_inactive) * ll$L0

    ll_h0 <- rowsum(obs_h0, dat$gene)
    ll_h1 <- rowsum(obs_h1, dat$gene)
    genes <- as.integer(rownames(ll_h0))

    h_g_new <- h_g
    h_g_new[genes] <- as.integer(ll_h1 > ll_h0)
    h_g_new
}

#' @keywords internal
.m_step_pi <- function(dat, post, h_g, pi_bounds = c(0.001, 0.499)) {
    counts_with_posterior <- dat %>%
        dplyr::left_join(post %>% dplyr::select(cell, post_X1_inactive), by = "cell") %>%
        dplyr::mutate(
            # Soft-assigned inactive-allele read count:
            # h=0 (X1 carries REF): X1-inactive cells contribute REF, X2-inactive contribute ALT
            # h=1 (X1 carries ALT): roles reversed
            xi_ref_count = ifelse(
                h_g[gene] == 0,
                post_X1_inactive * ref + (1 - post_X1_inactive) * alt,
                post_X1_inactive * alt + (1 - post_X1_inactive) * ref
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

#' Re-estimate phase and escape fraction for every gene against a frozen cell assignment
#'
#' \code{.infer_xci} discards phase (\code{h_g}) and escape fraction (\code{pi_g}) for any gene
#' that fails the outlier or escapee filters, because a gene uninformative for calling active-X
#' state has no business influencing that call. But once the call is made, those same genes can
#' still be phased against it: given the cells' active-X state as fixed, ground truth, fitting a
#' gene's own phase and escape fraction is a well-posed one-gene problem, unlike the joint EM
#' where a weakly-discriminative gene's own fit is nearly indifferent between phase orientations.
#' Alternates the existing M-steps (\code{\link{.m_step_phase}}, \code{\link{.m_step_pi}}) with no
#' E-step, since the point is to quantify escape without letting escapee data influence which cell
#' is called X1- or X2-active (that would be circular).
#'
#' @param ref_mat,alt_mat Full per-donor gene x cell count matrices (one row per het SNP/gene),
#'   unfiltered by the outlier or escapee filters.
#' @param post Frozen per-cell posterior from the informative-gene fit (\code{xci_result$post}),
#'   with \code{cell} and \code{post_X1_inactive} columns.
#' @param rho Beta-binomial overdispersion from the informative-gene fit (\code{xci_result$rho}).
#'
#' @return List with \code{h_g} (integer phase, 0 = REF on X1) and \code{pi_g} (numeric escape
#'   fraction), both length \code{nrow(ref_mat)}. Genes with no covered cell are \code{NA} in
#'   both, since they cannot be phased at all.
#'
#' @keywords internal
.rephase_all_genes <- function(ref_mat, alt_mat, post, rho, min_cov = 1, max_iter = 50, tol = 1e-4) {
    n_genes <- nrow(ref_mat)
    dat <- .pivot_counts_to_long(ref_mat, alt_mat, min_cov)
    # Only cells the informative-gene fit actually scored carry a posterior; an unscored cell
    # would default to 0 in .m_step_phase's post lookup, silently reading as "certainly
    # X2-active" and corrupting phase for genes uniquely covered there.
    dat <- dplyr::filter(dat, cell %in% post$cell)

    h_g <- rep(0L, n_genes)
    pi_g <- rep(0.05, n_genes)
    if (nrow(dat) > 0) {
        dedup <- .build_ll_dedup(dat)
        for (iter in seq_len(max_iter)) {
            ll <- .betabinom_ll_both(dedup, pi_g, rho)
            h_g <- .m_step_phase(dat, post, h_g, ll)
            pi_g_new <- .m_step_pi(dat, post, h_g)
            converged <- max(abs(pi_g_new - pi_g)) < tol
            pi_g <- pi_g_new
            if (converged) {
                break
            }
        }
    }

    # Genes absent from dat have no covered cell and cannot be phased at all, unlike the 0.05
    # filler .m_step_pi otherwise supplies for bookkeeping.
    covered <- seq_len(n_genes) %in% unique(dat$gene)
    list(
        h_g = ifelse(covered, h_g, NA_integer_),
        pi_g = ifelse(covered, pi_g, NA_real_)
    )
}

#' @keywords internal
.loglik_obs <- function(ref, n, p, rho) {
    # Beta-binomial: rho > 0 adds overdispersion relative to binomial; rho = 0 reduces to binomial
    VGAM::dbetabinom(ref, size = n, prob = p, rho = rho, log = TRUE)
}

#' Beta-binomial log-density without the binomial coefficient
#'
#' Inlined equivalent of \code{VGAM::dbetabinom(ref, n, p, rho, log = TRUE)}
#' with the \code{lchoose(n, ref)} term dropped. That term depends only on the
#' data (not on \code{p}), so it cancels in every comparison the EM makes: the
#' E-step log-odds ratio is a difference of log-likelihoods at the same
#' \code{(ref, n)}, and the M-phase compares the two phase orientations at the
#' same rows. Dropping it avoids a per-row \code{lchoose} on every iteration
#' and matches the VGAM result up to that additive constant.
#'
#' The beta-binomial shape parameters relate to the mean \code{p} and the
#' correlation \code{rho} by \eqn{a = p (1 - rho) / rho} and
#' \eqn{b = (1 - p)(1 - rho) / rho}.
#'
#' @keywords internal
.betabinom_ll_kernel <- function(ref, n, p, rho) {
    scale <- (1 - rho) / rho
    a <- p * scale
    b <- (1 - p) * scale
    lbeta(ref + a, n - ref + b) - lbeta(a, b)
}

#' Precompute the deduplication index for the log-likelihood kernel
#'
#' Read depths are small integers, so the kernel arguments \code{(ref, n, gene)}
#' take only a few thousand distinct values across the hundreds of thousands of
#' observations in \code{dat} (often ~40x fewer). Because \code{dat} is constant
#' across EM iterations, we identify the unique argument triples once here and
#' reuse the mapping every iteration, evaluating the kernel on the unique rows
#' and scattering the result back.
#'
#' @return A list with \code{idx} (row -> unique-row index) and the unique-row
#'   slices \code{ref}, \code{n}, \code{gene}.
#'
#' @keywords internal
.build_ll_dedup <- function(dat) {
    idx <- as.integer(factor(paste(dat$ref, dat$n, dat$gene)))
    keep <- !duplicated(idx)
    ord <- order(idx[keep])
    list(
        idx = idx,
        ref = dat$ref[keep][ord],
        n = dat$n[keep][ord],
        gene = dat$gene[keep][ord]
    )
}

#' Log-likelihood of each observation under both escape orientations
#'
#' Returns, for every row of \code{dat}, the beta-binomial log-likelihood of
#' the REF count under \code{pi_g[gene]} (\code{L0}) and under
#' \code{1 - pi_g[gene]} (\code{L1}). These are the only two probabilities the
#' EM ever evaluates, so the E-step and both M-phase orientations are built
#' from this single pair rather than repeated \code{dbetabinom} calls.
#'
#' The kernel is evaluated only on the unique \code{(ref, n, gene)} rows
#' identified by \code{dedup} and scattered back to full length, avoiding the
#' large redundancy in the raw observations.
#'
#' @keywords internal
.betabinom_ll_both <- function(dedup, pi_g, rho) {
    p0 <- pi_g[dedup$gene]
    L0u <- .betabinom_ll_kernel(dedup$ref, dedup$n, p0, rho)
    L1u <- .betabinom_ll_kernel(dedup$ref, dedup$n, 1 - p0, rho)
    list(
        L0 = L0u[dedup$idx],
        L1 = L1u[dedup$idx]
    )
}

#' @keywords internal
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

#' @keywords internal
.pivot_counts_to_long <- function(ref_mat, alt_mat, min_cov = 1) {
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

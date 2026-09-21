# Per-donor X-chromosome inactivation (XCI) inference.
#
# .filter_to_informative_het_snps() narrows a donor's SNPData down to one
# heterozygous chrX SNP per gene. .infer_xci() then fits a beta-binomial
# mixture model by EM over those genes' REF/ALT counts, jointly calling each
# gene's phase (which allele sits on the active X) and escape fraction, and
# each cell's active-X identity; it takes the best of several random
# restarts and flags genes that turned out uninformative for the call.
# .rephase_all_genes() is a sibling driver that re-fits phase and escape for
# every gene, including ones .infer_xci() excluded, against the cell
# assignment .infer_xci() already froze — so excluded genes are still scored
# without influencing that assignment. Both drivers share the EM machinery
# (.run_em() and its E/M-steps) and the beta-binomial kernel at the bottom of
# the file.

.filter_to_informative_het_snps <- function(snp_data, donor = NULL) {
    # Takes a single donor's SNPData and filters it down to chrX SNPs, then to
    # heterozygous SNPs, then to one SNP per gene (the highest-coverage one).
    # The donor label is used only for the log message. Returns the filtered
    # SNPData.

    # Suppress the generic per-filter logs from filter_snps; we emit a single
    # consolidated, donor-labelled summary instead.
    old_threshold <- logger::log_threshold()
    logger::log_threshold(logger::WARN)
    on.exit(logger::log_threshold(old_threshold), add = TRUE)

    n_start <- nrow(snp_info(snp_data))

    # Restrict to X-chromosome SNPs. XCI is an X-chromosome phenomenon, so SNPs
    # on other chromosomes carry no signal and would only add noise. Use the
    # canonical (UCSC) chromosome name set at construction so the match is robust
    # to the input naming style.
    snp_info <- snp_info(snp_data)
    if (!"chrom_canonical" %in% colnames(snp_info)) {
        stop(
            "No canonical chromosome names available. Ensure the SNPData object was built with a 'chrom' column so chrX SNPs can be selected."
        )
    }
    snp_data <- snp_data %>%
        filter_snps(chrom_canonical == "chrX")
    n_chrx <- nrow(snp_info(snp_data))

    het_snp_ids <- snp_data %>%
        donor_het_status_df() %>%
        dplyr::filter(zygosity == "het") %>%
        dplyr::pull(snp_id)

    snp_data <- snp_data %>%
        filter_snps(snp_id %in% het_snp_ids)
    n_het <- nrow(snp_info(snp_data))

    top_snp_per_gene <- snp_info(snp_data) %>%
        dplyr::arrange(dplyr::desc(coverage)) %>%
        dplyr::slice_head(n = 1, by = "gene_name")

    snp_data <- snp_data %>%
        filter_snps(snp_id %in% top_snp_per_gene$snp_id)
    n_genes <- nrow(snp_info(snp_data))

    logger::log_info(
        "[{donor}] het-SNP filter: {n_start} SNPs -> {n_chrx} chrX -> {n_het} het -> {n_genes} genes (top SNP per gene)"
    )

    snp_data
}

.infer_xci <- function(
    ref_mat,
    alt_mat,
    n_inits = 10,
    confidence_threshold = 0.95,
    min_cells = 10,
    min_cov = 1,
    donor = NULL
) {
    # Top-level per-donor driver, taking genes x cells REF/ALT count matrices.
    # Filters genes, then runs the EM (.run_em) from n_inits random restarts
    # and keeps the best by log-likelihood. Each cell's active-X call is then
    # labelled at confidence_threshold, and genes are flagged as informative
    # or not. Returns the best restart's list (post, h_g, pi_g, rho, prior,
    # ll; see .run_em), with post gaining an assignment column ("X1"/"X2"/
    # "unassigned") and a new gene_keep logical vector, length nrow(ref_mat),
    # that is TRUE for genes surviving both gene filters.
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

    # Post-convergence uninformative-gene filter: genes with LLR <= 0 are inconsistent
    # with current cell assignments; MAD filter on pi_g removes outlier escape fractions
    # (not incomplete escape in general — a gene with real but moderate escape keeps a
    # positive LLR and stays). This decides which genes get reported as xci_informative
    # (drove active-X calling) — it does not change best$post, so these genes'
    # contribution to the cell-assignment posterior above is already baked in either
    # way. Computed against the pre-relabel h_g/post (relabelling below is a symmetric
    # X1<->X2 swap and does not change which genes pass).
    passes_uninformative_filter <- .filter_uninformative_genes(dat, n_genes, best$h_g, best$pi_g, best$rho, best$post)

    # assignment names the *active* X, matching the stored active_x column.
    best$post <- best$post %>%
        dplyr::mutate(
            assignment = dplyr::case_when(
                post_X1_active >= confidence_threshold ~ "X1",
                post_X1_active <= 1 - confidence_threshold ~ "X2",
                TRUE ~ "unassigned"
            )
        )

    best <- .canonicalise_x1_label(best)

    # gene_keep: logical of length nrow(original ref_mat), TRUE = gene survived
    # both the outlier filter and the post-convergence uninformative-gene filter.
    gene_keep <- passes_outlier_filter
    gene_keep[passes_outlier_filter] <- passes_uninformative_filter

    c(best, list(gene_keep = gene_keep))
}

.canonicalise_x1_label <- function(fit) {
    # Takes a fit list as returned by .run_em. The EM's component labelling is
    # arbitrary, fixed only by the random phase initialisation, so this
    # relabels the fit so that X1 always denotes the active-X majority.
    # Returns the relabelled fit list.
    if (sum(fit$post$assignment == "X2") > sum(fit$post$assignment == "X1")) {
        fit$post <- fit$post %>%
            dplyr::mutate(
                post_X1_active = 1 - post_X1_active,
                assignment = dplyr::case_when(
                    assignment == "X1" ~ "X2",
                    assignment == "X2" ~ "X1",
                    TRUE ~ assignment
                )
            )
        fit$h_g <- 1L - fit$h_g
        fit$prior <- 1 - fit$prior
    }
    fit
}

.filter_uninformative_genes <- function(dat, n_genes, h_g, pi_g, rho, post, mad_threshold = 2) {
    # n_genes is the total gene count, including any already dropped by
    # .filter_outlier_genes. Drops genes with a non-positive log-likelihood
    # ratio (.compute_gene_llr), then drops genes with an outlying escape
    # fraction by robust z-score against mad_threshold. This only decides
    # xci_informative status and does not change the posterior. Returns a
    # logical vector of length n_genes.
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

.compute_gene_llr <- function(dat, post, h_g, pi_g, rho) {
    # Takes the long observation tibble (.pivot_counts_to_long), the per-cell
    # posterior (.e_step), and the fitted phase, escape fraction, and
    # overdispersion. For each gene, sums the posterior-weighted
    # log-likelihood-ratio of the fitted phase versus its flip. Returns a
    # tibble with one row per gene present in the data, giving the gene id and
    # its llr, where a positive value means the data supports the current
    # phase call.

    # Escape (pi_g) is a property of the silenced allele, so this contrast is
    # naturally expressed in inactive-X terms even though the posterior carried
    # alongside it is active-oriented.
    p_if_x1_inactive <- ifelse(h_g[dat$gene] == 0, pi_g[dat$gene], 1 - pi_g[dat$gene])
    p_if_x2_inactive <- 1 - p_if_x1_inactive

    ll_x1_inactive <- .loglik_obs(dat$ref, dat$n, p_if_x1_inactive, rho)
    ll_x2_inactive <- .loglik_obs(dat$ref, dat$n, p_if_x2_inactive, rho)

    dat %>%
        dplyr::left_join(post %>% dplyr::select(cell, post_X1_active), by = "cell") %>%
        # Expected LLR per observation: weighted by P(X1-inactive) for each cell.
        # Positive contribution means the observation is consistent with current assignments.
        dplyr::mutate(llr = (1 - post_X1_active) * (ll_x1_inactive - ll_x2_inactive)) %>%
        dplyr::group_by(gene) %>%
        dplyr::summarise(llr = sum(llr), .groups = "drop")
}

.filter_outlier_genes <- function(ref_mat, alt_mat, min_cells = 10, min_cov = 1, mad_threshold = 2) {
    # The pre-EM gene filter. Takes genes x cells REF/ALT count matrices,
    # drops genes covered in fewer than min_cells cells, then drops genes
    # whose REF/ALT skew, folded to the range [0.5, 1], is an outlier by
    # robust z-score against mad_threshold. Returns a logical vector of length
    # nrow(ref_mat).

    # Coverage per cell-gene pair
    n_mat <- ref_mat + alt_mat
    covered <- n_mat >= min_cov

    # Per gene: how many cells have sufficient coverage, and how many of those favour REF. A
    # tied cell (ref == alt) is split 0.5/0.5 rather than counted toward neither allele: the
    # strict ref > alt count silently pulls a gene's skew toward 0.5 by an amount that depends
    # on how many cells happen to tie, which is worst at low coverage — precisely where it is
    # most needed to tell moderate escape apart from a true 1:1 gene.
    n_expressing <- rowSums(covered)
    ref_majority <- rowSums(covered & (ref_mat > alt_mat)) + 0.5 * rowSums(covered & (ref_mat == alt_mat))

    # Drop genes seen in too few cells — not enough information to estimate allelic skew
    passes_count_filter <- n_expressing >= min_cells

    # Allelic skew: fraction of covered cells favouring REF, folded to [0.5, 1]
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

.pivot_counts_to_long <- function(ref_mat, alt_mat, min_cov = 1) {
    # Takes genes x cells REF/ALT count matrices and pivots them to one row
    # per cell-gene pair with at least min_cov total reads. Returns a tibble
    # with one row per covered pair, giving the gene and cell matrix indices,
    # the ref and alt counts, and their sum n.
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

.rephase_all_genes <- function(ref_mat, alt_mat, post, rho, min_cov = 1, max_iter = 50, tol = 1e-4) {
    # Re-fits phase and escape fraction for every gene against a frozen cell
    # assignment, so genes dropped from the EM as uninformative for calling
    # active-X are still phased and scored without influencing that call.
    # Alternates .m_step_phase and .m_step_pi with no E-step. Takes the full,
    # unfiltered per-donor count matrices, the frozen per-cell posterior from
    # the informative-gene fit (xci_result$post), and the overdispersion from
    # that fit (xci_result$rho). Returns a list with h_g (integer phase, 0
    # meaning REF is on X1) and pi_g (numeric escape fraction), both of length
    # nrow(ref_mat); both are NA for genes with no covered cell.
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

.run_em <- function(dat, n_genes, max_iter = 50, tol = 1e-4, init_seed = 1) {
    # Runs one EM restart from a random initial phase seeded by init_seed,
    # alternating the E-step (.e_step) with the four M-steps (.m_step_phase,
    # .m_step_pi, .m_step_rho, .m_step_prior) until pi_g, h_g and prior stop
    # changing by more than tol, or max_iter is reached. Returns a list
    # containing the final E-step tibble as post, the fitted phase per gene as
    # h_g, the fitted escape fraction per gene as pi_g, the fitted
    # overdispersion as rho, the fitted X1-active prior as prior, and the
    # final mixture log-likelihood as ll, used only to rank restarts.
    h_g <- withr::with_seed(init_seed, sample(0:1, n_genes, replace = TRUE)) # random phase initialisation
    pi_g <- rep(0.05, n_genes) # start conservatively: assume 5% escape fraction
    rho <- 0.05 # initial beta-binomial overdispersion, re-estimated by .m_step_rho below
    prior <- 0.5 # initial X1-active prior, re-estimated by .m_step_prior below

    # Unique kernel arguments are fixed across iterations; index them once.
    dedup <- .build_ll_dedup(dat)

    # Parameter-scale convergence (matches .rephase_all_genes): the observed-data
    # log-likelihood is a sum over every cell, so an absolute tolerance on it
    # never binds once cell counts get large -- max_iter would silently cap the
    # fit instead. pi_g, h_g and prior are what the fit is actually judged on,
    # and all three stabilise well before the log-likelihood sum stops moving
    # in its asymptotic tail.
    pi_g_prev <- pi_g
    h_g_prev <- h_g
    prior_prev <- prior
    for (iter in seq_len(max_iter)) {
        # Both escape orientations, computed once and shared by E- and M-steps.
        ll <- .betabinom_ll_both(dedup, pi_g, rho)
        post <- .e_step(dat, h_g, ll, prior = prior)
        h_g <- .m_step_phase(dat, post, h_g, ll)
        pi_g <- .m_step_pi(dat, post, h_g)
        rho <- .m_step_rho(dat, dedup, post, h_g, pi_g)
        prior <- .m_step_prior(post)
        if (max(abs(pi_g - pi_g_prev)) < tol && all(h_g == h_g_prev) && abs(prior - prior_prev) < tol) {
            break
        }
        pi_g_prev <- pi_g
        h_g_prev <- h_g
        prior_prev <- prior
    }
    # Observed-data mixture log-likelihood (see .e_step) at the final parameters,
    # used only to rank restarts in .infer_xci(); not part of the stopping rule.
    ll_final <- sum(post$ll_mix)
    list(post = post, h_g = h_g, pi_g = pi_g, rho = rho, prior = prior, ll = ll_final)
}

.e_step <- function(dat, h_g, ll, prior = 0.5) {
    # The E-step. Given the current phase per gene and the L0/L1 likelihoods
    # from .betabinom_ll_both, sums per-gene log-likelihoods for each cell to
    # get the log-odds that X1 is active, then converts that to a posterior
    # probability. prior is the prior probability that X1 is active. Returns a
    # tibble with one row per cell present in the data, giving the cell id,
    # the posterior post_X1_active, the underlying log-odds ratio lor, and the
    # mixture log-likelihood ll_mix, which is invariant to the X1/X2 labelling
    # and is used to rank EM restarts.

    # ll$L0 = loglik at pi_g[gene] (REF silenced), ll$L1 = loglik at 1 - pi_g[gene].
    # p(REF | X1 active): if h=0 (X1 carries REF), REF is expressed → 1 - pi_g → L1
    #                      if h=1 (X1 carries ALT), REF is silenced → pi_g     → L0
    # X2-active is the complement orientation, so its two cases are swapped.
    h_row <- h_g[dat$gene] == 0
    obs_x1_active <- ifelse(h_row, ll$L1, ll$L0)
    obs_x2_active <- ifelse(h_row, ll$L0, ll$L1)

    # Σ_g loglik per cell under each hypothesis, via rowsum (faster than group_by).
    cell_x1 <- rowsum(obs_x1_active, dat$cell)
    cell_x2 <- rowsum(obs_x2_active, dat$cell)

    logit_prior <- log(prior / (1 - prior)) # log(1) = 0 for equal prior
    lor <- as.numeric(cell_x1 - cell_x2) + logit_prior

    # Observed-data log-likelihood per cell, log(prior * P(obs | X1 active) +
    # (1 - prior) * P(obs | X2 active)), via log-sum-exp for stability. Unlike
    # the per-cell posterior, this is invariant to which mixture component is
    # labelled X1, so comparing it across EM restarts ranks them on fit quality
    # rather than on an arbitrary labelling. It is offset by a constant — the
    # lchoose(n, ref) term .betabinom_ll_kernel drops — which is identical across
    # restarts and labellings, so neither the argmax nor the convergence delta
    # is affected.
    a <- as.numeric(cell_x1) + log(prior)
    b <- as.numeric(cell_x2) + log1p(-prior)
    ll_mix <- pmax(a, b) + log1p(exp(-abs(a - b)))

    tibble::tibble(
        cell = as.integer(rownames(cell_x1)),
        post_X1_active = 1 / (1 + exp(-lor)), # sigmoid converts LOR to posterior
        lor = lor,
        ll_mix = ll_mix
    )
}

.check_post_covers_dat <- function(post_by_obs) {
    # Takes the per-observation posterior looked up by .m_step_phase,
    # .m_step_pi, or .m_step_rho. A cell that appears in the observations but
    # is absent from the posterior has no safe default, since any fill-in
    # value would read as a confident but wrong call. Stops with an error if
    # any value is NA; otherwise returns invisible(TRUE).
    if (anyNA(post_by_obs)) {
        stop(
            "Observations reference cells absent from `post`. Filter `dat` to the ",
            "cells the E-step scored before calling the M-steps."
        )
    }
    invisible(TRUE)
}

.m_step_phase <- function(dat, post, h_g, ll) {
    # The M-step for phase. Takes the observations, the per-cell posterior
    # (which must cover every cell id present in the observations, checked by
    # .check_post_covers_dat), the current phase per gene, and the L0/L1
    # likelihoods. For each gene, compares the expected log-likelihood of
    # phase 0 against phase 1 and picks whichever is higher. Returns an
    # integer vector of length n_genes giving the updated phase per gene
    # (0 means REF is on X1, 1 means ALT is on X1); genes with no covered cell
    # keep their previous value.

    # Expected log-likelihood under each phase: E_q[log p(ref | h, pi_g)], built
    # from the shared orientation pair. ll$L0 = loglik at pi_g (REF silenced),
    # ll$L1 = loglik at 1 - pi_g (REF active).
    # h=0: X1-active cells (weight post_X1_active) see REF fraction 1 - pi_g → L1;
    #       X2-active cells (weight 1-post_X1_active) see pi_g → L0.
    # h=1: the two orientations swap.
    # Scatter cell posteriors into a lookup indexed by cell id. post$cell is a
    # sorted subset of cell ids (cells with no covered gene are absent), so a
    # positional index would misalign — index by id instead. The fill is NA, not
    # 0, so an uncovered cell is rejected rather than silently read as a
    # confident call — see .check_post_covers_dat().
    post_lookup <- rep(NA_real_, max(dat$cell))
    post_lookup[post$cell] <- post$post_X1_active
    post_X1_active <- post_lookup[dat$cell]
    .check_post_covers_dat(post_X1_active)

    obs_h0 <- (1 - post_X1_active) * ll$L0 + post_X1_active * ll$L1
    obs_h1 <- (1 - post_X1_active) * ll$L1 + post_X1_active * ll$L0

    ll_h0 <- rowsum(obs_h0, dat$gene)
    ll_h1 <- rowsum(obs_h1, dat$gene)
    genes <- as.integer(rownames(ll_h0))

    h_g_new <- h_g
    h_g_new[genes] <- as.integer(ll_h1 > ll_h0)
    h_g_new
}

.m_step_pi <- function(dat, post, h_g, pi_bounds = c(0.001, 0.499)) {
    # The M-step for escape fraction. Takes the observations, the per-cell
    # posterior (covering every cell id present in the observations), and the
    # current phase per gene, which determines which allele is on the
    # inactive X in each cell. Computes the closed-form MLE as the
    # soft-assigned inactive-allele read count divided by total reads per
    # gene, clamped to pi_bounds so the result stays interpretable as a minor
    # fraction. Returns a numeric vector of length length(h_g) giving the
    # updated escape fraction per gene; genes with no covered cell default to
    # 0.05.

    # left_join leaves NA for any cell absent from post; reject rather than let
    # it propagate into pi_g (see .check_post_covers_dat).
    joined <- dat %>%
        dplyr::left_join(post %>% dplyr::select(cell, post_X1_active), by = "cell")
    .check_post_covers_dat(joined$post_X1_active)

    counts_with_posterior <- joined %>%
        dplyr::mutate(
            # Soft-assigned inactive-allele read count. pi_g stays the escape
            # (inactive-allele) fraction regardless of how the cell posterior is
            # oriented, so the silenced allele is the one on the *inactive* X:
            # h=0 (X1 carries REF): X2-active cells contribute REF, X1-active contribute ALT
            # h=1 (X1 carries ALT): roles reversed
            xi_ref_count = ifelse(
                h_g[gene] == 0,
                (1 - post_X1_active) * ref + post_X1_active * alt,
                (1 - post_X1_active) * alt + post_X1_active * ref
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

.m_step_prior <- function(post, prior_bounds = c(1e-3, 1 - 1e-3)) {
    # The M-step for the X1-active prior. Takes the per-cell posterior and
    # returns the closed-form MLE for a mixing weight: the average posterior
    # probability of X1 across cells, clamped to prior_bounds to keep
    # log(prior/(1-prior)) in .e_step finite even for a donor whose cells are
    # almost entirely on one X.
    prior <- mean(post$post_X1_active)
    pmax(prior_bounds[1], pmin(prior_bounds[2], prior))
}

.m_step_rho <- function(dat, dedup, post, h_g, pi_g, rho_bounds = c(1e-4, 0.999)) {
    # The M-step for overdispersion. Takes the observations, the
    # deduplication index from .build_ll_dedup, the per-cell posterior, and
    # the current phase and escape fraction per gene. Fits a single
    # dataset-wide rho by maximising the expected log-likelihood over
    # rho_bounds with everything else held fixed. Returns a numeric scalar
    # within rho_bounds.

    # Same posterior lookup/weighting pattern as .m_step_phase: expected
    # log-likelihood under each cell's current phase call, but as a function of
    # rho rather than h_g. A single dataset-wide scalar rather than per-gene --
    # per-gene estimation would trade off against pi_g on typically sparse
    # single-cell counts, especially for low-coverage genes.
    post_lookup <- rep(NA_real_, max(dat$cell))
    post_lookup[post$cell] <- post$post_X1_active
    post_X1_active <- post_lookup[dat$cell]
    .check_post_covers_dat(post_X1_active)

    neg_expected_ll <- function(rho) {
        ll <- .betabinom_ll_both(dedup, pi_g, rho)
        obs_h0 <- (1 - post_X1_active) * ll$L0 + post_X1_active * ll$L1
        obs_h1 <- (1 - post_X1_active) * ll$L1 + post_X1_active * ll$L0
        chosen <- ifelse(h_g[dat$gene] == 0, obs_h0, obs_h1)
        -sum(chosen)
    }
    stats::optimize(neg_expected_ll, interval = rho_bounds)$minimum
}

.loglik_obs <- function(ref, n, p, rho) {
    # Wraps VGAM::dbetabinom for use outside the EM's hot loop. ref and n are
    # the observed REF counts and totals, p is the expected REF fraction, and
    # rho is the overdispersion, where 0 reduces to the binomial. Returns a
    # numeric vector of log-densities the same length as ref.
    VGAM::dbetabinom(ref, size = n, prob = p, rho = rho, log = TRUE)
}

.betabinom_ll_kernel <- function(ref, n, p, rho) {
    # Computes the beta-binomial log-density without the binomial
    # coefficient. This is equivalent to VGAM::dbetabinom(ref, n, p, rho,
    # log = TRUE) minus the lchoose(n, ref) term, which cancels in every
    # comparison the EM makes and is dropped here to avoid a per-row lchoose
    # on every iteration. The shape parameters relate to p and rho by
    # a = p(1 - rho)/rho and b = (1 - p)(1 - rho)/rho.
    scale <- (1 - rho) / rho
    a <- p * scale
    b <- (1 - p) * scale
    lbeta(ref + a, n - ref + b) - lbeta(a, b)
}

.build_ll_dedup <- function(dat) {
    # Indexes the unique (ref, n, gene) argument triples in the observations
    # once, since read depths are small integers and the observations are
    # constant across EM iterations, often giving ~40x fewer unique rows than
    # observations. Returns a list with idx, mapping each row to its
    # unique-row index, and the unique-row slices ref, n, and gene.
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

.betabinom_ll_both <- function(dedup, pi_g, rho) {
    # Takes the deduplication index from .build_ll_dedup and the current
    # escape fraction and overdispersion, evaluates the kernel only on the
    # unique rows, and scatters the result back to full length. Returns a list
    # with L0, the log-likelihood of the REF count under pi_g[gene], and L1,
    # under 1 - pi_g[gene], each a numeric vector matching the length of the
    # original observations.
    p0 <- pi_g[dedup$gene]
    L0u <- .betabinom_ll_kernel(dedup$ref, dedup$n, p0, rho)
    L1u <- .betabinom_ll_kernel(dedup$ref, dedup$n, 1 - p0, rho)
    list(
        L0 = L0u[dedup$idx],
        L1 = L1u[dedup$idx]
    )
}

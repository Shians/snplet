# ==============================================================================
# Test Suite: X-chromosome inactivation fitting
# Description: assign_xci cell/clonotype modes, stored diagnostics, and
#              the SNPData accessor/heatmap methods.
# ==============================================================================

library(testthat)
library(Matrix)

# ------------------------------------------------------------------------------
# Test Data Setup
# ------------------------------------------------------------------------------

# Build a synthetic SNPData with a strong X-inactivation signal for one donor.
# Cells fall into two clonally-consistent groups. Within a group each gene
# expresses predominantly one allele; pooled across all cells each SNP looks
# heterozygous (balanced), so it survives the het-SNP filter, while per-cell
# the allele is skewed — exactly the XCI pattern the EM should recover.
make_xci_snpdata <- function(
    n_genes = 20,
    n_cells_per_group = 40,
    depth = 20,
    seed = 1,
    n_donors = 1,
    escapee = FALSE
) {
    withr::with_seed(seed, {
        n_cells <- 2 * n_cells_per_group
        # The X that stays active in each cell, i.e. the truth assign_xci() should
        # recover in its active_x column.
        group <- rep(c("X1", "X2"), each = n_cells_per_group)

        # Per gene, which allele sits on X1 (0 = REF sits on X1, so REF is the
        # expressed allele in X1-active cells).
        # When escapee = TRUE the last gene is an escapee: it stays balanced
        # (p_ref = 0.5) in every cell regardless of inactivation state, so it
        # carries no XCI signal and should be filtered out.
        allele_on_x1 <- sample(0:1, n_genes, replace = TRUE)
        is_escapee <- rep(FALSE, n_genes)
        if (escapee) {
            is_escapee[n_genes] <- TRUE
        }

        ref_mat <- matrix(0L, nrow = n_genes, ncol = n_cells)
        alt_mat <- matrix(0L, nrow = n_genes, ncol = n_cells)
        for (g in seq_len(n_genes)) {
            for (c in seq_len(n_cells)) {
                # X1-active cells express the X1 allele; probability of REF
                x1_active <- group[c] == "X1"
                p_ref <- if (is_escapee[g]) {
                    0.5
                } else if (allele_on_x1[g] == 0) {
                    if (x1_active) 0.95 else 0.05
                } else {
                    if (x1_active) 0.05 else 0.95
                }
                ref_mat[g, c] <- rbinom(1, depth, p_ref)
                alt_mat[g, c] <- depth - ref_mat[g, c]
            }
        }

        ref_mat <- Matrix(ref_mat, sparse = TRUE)
        alt_mat <- Matrix(alt_mat, sparse = TRUE)

        snp_info <- data.frame(
            chrom = "X",
            pos = seq_len(n_genes) * 1000L,
            ref = "A",
            alt = "G",
            gene_name = paste0("gene", seq_len(n_genes)),
            stringsAsFactors = FALSE
        )
        # Many small clonotypes per group so the clonotype-level EM has enough
        # units to survive the min_cells outlier filter. Clonotypes never span
        # groups, preserving clonal consistency of X-inactivation state.
        clono_within <- rep(seq_len(n_cells_per_group %/% 4 + 1), each = 4)[seq_len(n_cells_per_group)]
        clonotype <- paste0("clono_", group, "_", c(clono_within, clono_within))
        # Spread cells evenly across donors. Clonotypes stay donor-specific by
        # tagging the donor into the clonotype id, preserving clonal consistency.
        donor <- paste0("donor", (seq_len(n_cells) - 1) %% n_donors)
        clonotype <- paste0(donor, "_", clonotype)
        barcode_info <- data.frame(
            barcode = paste0("cell", seq_len(n_cells)),
            donor = donor,
            clonotype = clonotype,
            true_group = group,
            stringsAsFactors = FALSE
        )

        snpdata <- SNPData(
            ref_count = ref_mat,
            alt_count = alt_mat,
            snp_info = snp_info,
            barcode_info = barcode_info
        )
        # No genotype source is supplied, so establish an active zygosity
        # source via the binomial test before any assign_xci()/
        # donor_het_status_df() call, which now require one.
        snpdata <- infer_zygosity(snpdata)

        list(
            snpdata = snpdata,
            allele_on_x1 = allele_on_x1
        )
    })
}

# ==============================================================================

test_that("assign_xci recovers the true clonal split at cell level", {
    fixture <- make_xci_snpdata()
    stored <- assign_xci(fixture$snpdata, n_inits = 3)

    assignments <- xci_assignments(stored)
    # Verify assignments carry the stored annotation columns
    expect_true(all(c("cell_id", "active_x", "xci_post_X1_active") %in% colnames(assignments)))
    # Verify every cell received an assignment row
    expect_equal(nrow(assignments), ncol(fixture$snpdata))

    # Confirm the two inferred groups recover the true clonal split (up to a
    # label swap, since X1/X2 labels are exchangeable in the model). Compare
    # only cells that met the confidence threshold.
    truth <- barcode_info(fixture$snpdata)$true_group
    called <- !is.na(assignments$active_x)
    agree <- mean(assignments$active_x[called] == truth[called])
    # Confirm agreement with the true clonal split exceeds 90% either way
    expect_true(max(agree, 1 - agree) > 0.9)
})

test_that("assign_xci always labels the majority active-X group X1", {
    # Which physical haplotype the EM's random phase initialisation happens to
    # call "X1" is otherwise arbitrary and can flip between runs. Build an
    # imbalanced split (most cells true X2) across several data seeds so the
    # majority group is not the one the EM would land on by chance, and check
    # it is always relabelled X1 regardless.
    for (seed in 1:3) {
        fixture <- make_xci_snpdata(seed = seed)
        bi <- barcode_info(fixture$snpdata)
        minority_cells <- bi$barcode[bi$true_group == "X1"][1:10]
        imbalanced <- filter_barcodes(fixture$snpdata, true_group == "X2" | barcode %in% minority_cells)

        stored <- assign_xci(imbalanced, n_inits = 3)
        called <- xci_assignments(stored)$active_x
        called <- called[!is.na(called)]
        tab <- table(called)

        # Confirm the larger group is always reported as X1
        expect_equal(names(which.max(tab)), "X1")
    }
})

test_that("assign_xci labels active_x by the allele the cell expresses", {
    fixture <- make_xci_snpdata()
    stored <- assign_xci(fixture$snpdata, n_inits = 3)

    # The X1/X2 labels are exchangeable, so the recovery tests above pass under
    # either polarity. Phase pins the direction: allele_on_x1 and active_x share
    # the same X1 label, so a cell calling X1 active must show X1's allele. This
    # is the property the heatmap relies on, and it fails if the inactive-to-
    # active flip in .store_xci_fit is dropped or inverted.
    haplotypes <- xci_haplotypes(stored)
    snp_rows <- match(haplotypes$snp_id, snp_info(stored)$snp_id)
    ref <- as.matrix(ref_count(stored))[snp_rows, , drop = FALSE]
    alt <- as.matrix(alt_count(stored))[snp_rows, , drop = FALSE]

    assignments <- xci_assignments(stored)
    x1_active <- which(assignments$active_x %in% "X1")
    x2_active <- which(assignments$active_x %in% "X2")
    ref_on_x1 <- haplotypes$allele_on_x1 == "REF"

    # Pooled REF fraction over a set of SNPs x cells
    ref_fraction <- function(snps, cells) {
        r <- sum(ref[snps, cells])
        r / (r + sum(alt[snps, cells]))
    }

    # Verify X1-active cells express X1's allele: REF-heavy where X1 carries REF
    expect_gt(ref_fraction(ref_on_x1, x1_active), 0.8)
    # Verify the same cells are ALT-heavy where X1 instead carries ALT
    expect_lt(ref_fraction(!ref_on_x1, x1_active), 0.2)
    # Confirm X2-active cells mirror it, expressing X2's allele at both SNP sets
    expect_lt(ref_fraction(ref_on_x1, x2_active), 0.2)
    expect_gt(ref_fraction(!ref_on_x1, x2_active), 0.8)
})

test_that("xci_haplotypes reports phase and escape fraction per informative SNP", {
    fixture <- make_xci_snpdata()
    stored <- assign_xci(fixture$snpdata, n_inits = 3)

    haplotypes <- xci_haplotypes(stored)
    # Verify haplotype columns are present, including the per-donor donor column
    expect_true(all(c("snp_id", "gene_name", "donor", "allele_on_x1", "escape_fraction") %in% colnames(haplotypes)))
    # Check phase is reported as REF/ALT
    expect_true(all(haplotypes$allele_on_x1 %in% c("REF", "ALT")))
    # Confirm escape fraction is a valid minor fraction
    expect_true(all(haplotypes$escape_fraction > 0 & haplotypes$escape_fraction < 0.5))
})

test_that("xci_haplotypes emits one row per SNP and donor for a multi-donor fit", {
    fixture <- make_xci_snpdata(n_donors = 2)
    stored <- assign_xci(fixture$snpdata, n_inits = 3)

    haplotypes <- xci_haplotypes(stored)
    # Confirm both donors are represented in the unnested output
    expect_setequal(haplotypes$donor, c("donor0", "donor1"))
    # Each row keys on one donor, so no donor cell is comma-joined
    expect_false(any(grepl(",", haplotypes$donor)))
    # A SNP informative in both donors yields one row per donor
    per_snp <- table(haplotypes$snp_id)
    expect_true(all(per_snp <= 2))
    # Escape fraction parses back to a numeric minor fraction in every donor
    expect_type(haplotypes$escape_fraction, "double")
    expect_true(all(haplotypes$escape_fraction > 0 & haplotypes$escape_fraction < 0.5))
})

test_that("clonotype fit projects assignments back to cells consistently", {
    fixture <- make_xci_snpdata()

    # Reach into the per-donor engine to inspect the pre-storage clonotype fit.
    donor_data <- filter_samples(fixture$snpdata, donor == "donor0")
    fit <- snplet:::.fit_xci_donor(donor_data, n_inits = 3, by = "clonotype")

    # Verify the clonotype fit records its unit
    expect_equal(fit$unit, "clonotype")
    # Verify the fit carries a cell-level projection of assignments
    expect_true(!is.null(fit$cell_assignments))

    # The EM ran per clonotype, one assignment row each
    n_clonotypes <- dplyr::n_distinct(barcode_info(fixture$snpdata)$clonotype)
    # Confirm one assignment row exists per distinct clonotype
    expect_equal(nrow(fit$assignments), n_clonotypes)

    # Cells within a clonotype must share an assignment (clonal consistency)
    cell_assign <- fit$cell_assignments
    clono_map <- dplyr::select(barcode_info(fixture$snpdata), cell_id, clonotype)
    joined <- dplyr::inner_join(cell_assign, clono_map, by = "cell_id")
    per_clono <- tapply(joined$assignment, joined$clonotype, function(a) length(unique(a)))
    # Confirm every clonotype's cells were all projected to the same assignment
    expect_true(all(per_clono == 1))
})

test_that("assign_xci promotes diagnostics into SNPData slots and survives subsetting", {
    fixture <- make_xci_snpdata()
    stored <- assign_xci(fixture$snpdata, n_inits = 3)

    # Verify a SNPData object is returned
    expect_s4_class(stored, "SNPData")

    barcode_info <- barcode_info(stored)
    snp_info <- snp_info(stored)
    donor_snp_info <- donor_snp_info(stored)
    donor_info <- donor_info(stored)
    # Check barcode diagnostics were written
    expect_true(all(c("active_x", "xci_post_X1_active") %in% colnames(barcode_info)))
    # Check per-donor diagnostics were written to donor_info
    expect_true(all(c("xci_rho", "xci_median_pi_g") %in% colnames(donor_info)))
    # Confirm rho and the median escape fraction are valid, non-NA values
    expect_true(all(!is.na(donor_info$xci_rho) & donor_info$xci_rho >= 0 & donor_info$xci_rho < 1))
    expect_true(all(!is.na(donor_info$xci_median_pi_g)))
    # Check the per-SNP informative flag was NOT duplicated into snp_info; it
    # belongs solely to donor_snp_info since informativeness is donor-specific
    expect_false("xci_informative" %in% colnames(snp_info))
    # Check per-(SNP, donor) diagnostics were written to donor_snp_info
    expect_true(all(
        c("snp_id", "donor", "xci_informative", "allele_on_x1", "xci_escape_fraction") %in%
            colnames(donor_snp_info)
    ))
    # Confirm some (SNP, donor) pairs are flagged informative and none are NA
    expect_true(any(donor_snp_info$xci_informative))
    # Confirm the informative flag itself is never NA
    expect_false(any(is.na(donor_snp_info$xci_informative)))
    # donor_snp_info now carries one row per phased (SNP, donor) pair, not just
    # the informative subset that drove active-X calling: confirm phase/escape
    # fraction are never NA (every stored row was genuinely phased) and that at
    # least one row is phased but not informative (haplotype_expression()'s
    # broader default gene set).
    expect_false(any(is.na(donor_snp_info$allele_on_x1)))
    expect_false(any(is.na(donor_snp_info$xci_escape_fraction)))
    expect_true(any(!donor_snp_info$xci_informative))

    # Diagnostics must survive cell subsetting because they live in barcode_info
    subset_cells <- stored[, 1:10]
    # Verify the active_x column is retained after subsetting cells
    expect_true("active_x" %in% colnames(barcode_info(subset_cells)))
    # Verify barcode_info row count matches the subset size
    expect_equal(nrow(barcode_info(subset_cells)), 10)

    # And survive SNP subsetting because donor_snp_info rows are restricted
    # to the kept SNPs
    informative_snp_ids <- donor_snp_info$snp_id[donor_snp_info$xci_informative]
    subset_snps <- filter_snps(stored, snp_id %in% informative_snp_ids)
    # donor_snp_info rows survive SNP subsetting too, restricted to the kept SNPs
    expect_true(all(donor_snp_info(subset_snps)$snp_id %in% snp_info(subset_snps)$snp_id))
})

test_that("assign_xci diagnostics drop a donor's donor_snp_info rows once its cells are gone", {
    fixture <- make_xci_snpdata(n_donors = 2)
    stored <- assign_xci(fixture$snpdata, n_inits = 3)

    donor0_cells <- barcode_info(stored)$donor == "donor0"
    subset_obj <- stored[, donor0_cells]

    # Confirm donor1 no longer appears in donor_info
    expect_false("donor1" %in% donor_info(subset_obj)$donor)
    # Confirm donor1's donor_snp_info rows were dropped along with its cells
    expect_false("donor1" %in% donor_snp_info(subset_obj)$donor)
})

test_that("accessors and heatmap work on a stored SNPData object", {
    fixture <- make_xci_snpdata()
    stored <- assign_xci(fixture$snpdata, n_inits = 3)

    assignments <- xci_assignments(stored)
    # Verify stored-object assignments expose the annotation columns
    expect_true(all(c("cell_id", "active_x", "xci_post_X1_active") %in% colnames(assignments)))

    haplotypes <- xci_haplotypes(stored)
    # Verify stored-object haplotypes expose donor, phase and escape fraction
    expect_true(all(c("snp_id", "donor", "allele_on_x1", "escape_fraction") %in% colnames(haplotypes)))

    # Confirm the heatmap method runs on a stored object. It returns a drawn
    # HeatmapList (donor as title, unit-labelled column axis), so draw to a null
    # device to keep the test headless.
    grDevices::pdf(NULL)
    on.exit(grDevices::dev.off(), add = TRUE)
    hm <- plot_xci_heatmap(stored, donor = "donor0")
    # Verify the heatmap method returns a drawn HeatmapList object
    expect_s4_class(hm, "HeatmapList")
})

# Shared stored-object fixture for heatmap display-parameter tests below
create_stored_xci_fixture <- function() {
    fixture <- make_xci_snpdata()
    assign_xci(fixture$snpdata, n_inits = 3)
}

test_that("heatmap max_genes caps the number of displayed rows", {
    stored <- create_stored_xci_fixture()
    grDevices::pdf(NULL)
    on.exit(grDevices::dev.off(), add = TRUE)

    capped <- plot_xci_heatmap(stored, donor = "donor0", max_genes = 5)

    # Confirm max_genes caps the number of displayed rows
    expect_equal(nrow(capped@ht_list[[1]]@matrix), 5)
})

test_that("heatmap max_genes above the available count shows all retained genes", {
    stored <- create_stored_xci_fixture()
    grDevices::pdf(NULL)
    on.exit(grDevices::dev.off(), add = TRUE)

    full <- plot_xci_heatmap(stored, donor = "donor0")
    n_all <- nrow(full@ht_list[[1]]@matrix)
    big <- plot_xci_heatmap(stored, donor = "donor0", max_genes = n_all + 100)

    # Confirm requesting more genes than available shows all retained genes
    expect_equal(nrow(big@ht_list[[1]]@matrix), n_all)
})

test_that("heatmap show_unassigned = FALSE drops unassigned columns", {
    stored <- create_stored_xci_fixture()
    grDevices::pdf(NULL)
    on.exit(grDevices::dev.off(), add = TRUE)

    assigned_only <- plot_xci_heatmap(
        stored,
        donor = "donor0",
        show_unassigned = FALSE
    )
    n_assigned <- sum(!is.na(barcode_info(stored)$active_x))

    # Check show_unassigned = FALSE drops unassigned columns
    expect_equal(ncol(assigned_only@ht_list[[1]]@matrix), n_assigned)
})

test_that("heatmap accepts toggling gene names and row clustering", {
    stored <- create_stored_xci_fixture()
    grDevices::pdf(NULL)
    on.exit(grDevices::dev.off(), add = TRUE)

    variant <- plot_xci_heatmap(
        stored,
        donor = "donor0",
        show_gene_names = FALSE,
        cluster_rows = TRUE
    )

    # Ensure toggling gene names and clustering still returns a valid plot
    expect_s4_class(variant, "HeatmapList")
})

test_that("heatmap accepts toggling assignment boundary markers", {
    stored <- create_stored_xci_fixture()
    grDevices::pdf(NULL)
    on.exit(grDevices::dev.off(), add = TRUE)

    no_marks <- plot_xci_heatmap(
        stored,
        donor = "donor0",
        mark_boundaries = FALSE
    )

    # Ensure toggling the assignment boundary markers still returns a valid plot
    expect_s4_class(no_marks, "HeatmapList")
})

test_that("heatmap show_posterior toggles the posterior annotation row", {
    stored <- create_stored_xci_fixture()
    grDevices::pdf(NULL)
    on.exit(grDevices::dev.off(), add = TRUE)

    with_post <- plot_xci_heatmap(stored, donor = "donor0")
    without_post <- plot_xci_heatmap(
        stored,
        donor = "donor0",
        show_posterior = FALSE
    )
    with_names <- names(with_post@ht_list[[1]]@top_annotation@anno_list)
    without_names <- names(without_post@ht_list[[1]]@top_annotation@anno_list)

    # Confirm the posterior_X1_active annotation is present by default
    expect_true("posterior_X1_active" %in% with_names)
    # Confirm the posterior_X1_active annotation is dropped when show_posterior = FALSE
    expect_false("posterior_X1_active" %in% without_names)
    # Confirm the assignment annotation is retained either way
    expect_true("assignment" %in% without_names)
})

test_that("heatmap applies custom colour arguments to the plot body", {
    stored <- create_stored_xci_fixture()
    grDevices::pdf(NULL)
    on.exit(grDevices::dev.off(), add = TRUE)

    coloured <- plot_xci_heatmap(
        stored,
        donor = "donor0",
        ref_fraction_palette = c("#2166ac", "#f7f7f7", "#b2182b"),
        assignment_palette = c(X1 = "#1b9e77", X2 = "#d95f02", unassigned = "grey85"),
        posterior_palette = c("purple", "white", "orange"),
        na_fill = "white"
    )
    body_col <- coloured@ht_list[[1]]@matrix_color_mapping@col_fun

    # Ensure custom colour arguments are accepted and applied to the body
    expect_s4_class(coloured, "HeatmapList")
    # Confirm the supplied ramp anchors drive the REF fraction body colours
    expect_equal(body_col(0), "#2166ACFF")
    # Confirm the high end of the ramp maps to the supplied top colour
    expect_equal(body_col(1), "#B2182BFF")
    # Confirm na_col was overridden (normalised to hex with alpha)
    expect_equal(coloured@ht_list[[1]]@matrix_color_mapping@na_col, "#FFFFFFFF")
})

test_that("heatmap distinguishes no-coverage units from low-confidence unassigned", {
    fixture <- make_xci_snpdata()
    stored <- assign_xci(fixture$snpdata, n_inits = 3)

    # Force one cell to look like it was never scored (NA posterior): the model
    # gives it no coverage, distinct from a low-confidence unassigned call.
    bi <- barcode_info(stored)
    bi$active_x[1] <- NA
    bi$xci_post_X1_active[1] <- NA_real_
    stored <- add_barcode_metadata(
        stored,
        dplyr::select(bi, cell_id, active_x, xci_post_X1_active),
        join_by = "cell_id",
        overwrite = TRUE
    )

    grDevices::pdf(NULL)
    on.exit(grDevices::dev.off(), add = TRUE)

    # With show_no_coverage = TRUE the annotation carries a distinct level
    hm <- plot_xci_heatmap(
        stored,
        donor = "donor0",
        show_no_coverage = TRUE
    )
    ann <- hm@ht_list[[1]]@top_annotation@anno_list$assignment
    ann_colors <- ann@color_mapping@colors
    # Verify the "no coverage" level is exposed when show_no_coverage = TRUE
    expect_true("no coverage" %in% names(ann_colors))
    # Check "no coverage" uses the dark slate distinct from the default
    # low-confidence "unassigned" grey70 (#B3B3B3)
    expect_equal(toupper(ann_colors[["no coverage"]]), "#4D4D4DFF")

    # By default (show_no_coverage = FALSE) the no-coverage cell is dropped, so
    # the plotted matrix loses exactly that one column
    default_hm <- plot_xci_heatmap(stored, donor = "donor0")
    n_donor0 <- sum(barcode_info(stored)$donor == "donor0")
    expect_equal(ncol(default_hm@ht_list[[1]]@matrix), n_donor0 - 1)
})

test_that("clonotype fit requires clonotype information", {
    fixture <- make_xci_snpdata()
    # Strip clonotype information
    no_clono <- add_barcode_metadata(
        fixture$snpdata,
        data.frame(
            cell_id = barcode_info(fixture$snpdata)$cell_id,
            clonotype = NA_character_
        ),
        overwrite = TRUE
    )
    # Verify a clear error is raised when clonotypes are all NA
    expect_error(
        assign_xci_by_clonotype(no_clono),
        "clonotype"
    )
})

test_that("inlined beta-binomial kernel matches VGAM up to the dropped lchoose term", {
    withr::with_seed(42, {
        ref <- rpois(2000, 4)
        size <- ref + rpois(2000, 8) + 1L
        p <- runif(2000, 0.05, 0.45)
    })
    rho <- 0.05

    kernel <- snplet:::.betabinom_ll_kernel(ref, size, p, rho)
    vgam <- VGAM::dbetabinom(ref, size = size, prob = p, rho = rho, log = TRUE)

    # The kernel drops lchoose(size, ref); the difference must equal exactly that
    # constant, confirming the two agree on every EM comparison (where it cancels).
    expect_equal(vgam - kernel, lchoose(size, ref), tolerance = 1e-10)
})

test_that("deduplicated kernel evaluation matches the direct per-row computation", {
    withr::with_seed(3, {
        n_genes <- 20
        dat <- tibble::tibble(
            gene = sample(n_genes, 8000, replace = TRUE),
            cell = sample(1500, 8000, replace = TRUE),
            ref = rpois(8000, 4),
            n = rpois(8000, 8) + 2
        )
        dat$ref <- pmin(dat$ref, dat$n)
        pi_g <- runif(n_genes, 0.05, 0.45)
    })
    rho <- 0.05

    dedup <- snplet:::.build_ll_dedup(dat)
    both <- snplet:::.betabinom_ll_both(dedup, pi_g, rho)

    direct_L0 <- snplet:::.betabinom_ll_kernel(dat$ref, dat$n, pi_g[dat$gene], rho)
    direct_L1 <- snplet:::.betabinom_ll_kernel(dat$ref, dat$n, 1 - pi_g[dat$gene], rho)

    # Scatter-back must reproduce the direct evaluation exactly (same values)
    expect_identical(both$L0, direct_L0)
    # Confirm the L1 (flipped-probability) branch also matches the direct computation
    expect_identical(both$L1, direct_L1)
})

# A minimal (dat, post) pair where cell 3 appears in the observations but was
# never scored, so it carries no posterior. Both M-steps must refuse it: any
# numeric stand-in asserts a confident call and would bias phase and escape for
# genes that cell alone covers.
make_unscored_cell_case <- function() {
    list(
        dat = tibble::tibble(
            gene = c(1L, 1L, 2L, 2L),
            cell = c(1L, 2L, 2L, 3L),
            ref = c(8L, 1L, 7L, 9L),
            alt = c(1L, 8L, 2L, 1L),
            n = c(9L, 9L, 9L, 10L)
        ),
        post = tibble::tibble(cell = c(1L, 2L), post_X1_active = c(0.99, 0.01)),
        h_g = c(0L, 0L)
    )
}

test_that(".m_step_phase() rejects an observation whose cell was never scored", {
    case <- make_unscored_cell_case()
    ll <- snplet:::.betabinom_ll_both(snplet:::.build_ll_dedup(case$dat), c(0.05, 0.05), 0.05)

    # Verify the unscored cell is refused rather than silently defaulting to a
    # confident posterior
    expect_error(
        snplet:::.m_step_phase(case$dat, case$post, case$h_g, ll),
        "cells absent from `post`"
    )
})

test_that(".m_step_pi() rejects an observation whose cell was never scored", {
    case <- make_unscored_cell_case()

    # Verify the same precondition holds in the escape-fraction M-step, whose
    # join would otherwise leak NA into pi_g
    expect_error(
        snplet:::.m_step_pi(case$dat, case$post, case$h_g),
        "cells absent from `post`"
    )
})

test_that(".rephase_all_genes() phases against scored cells and leaves the rest NA", {
    # Genes 1-2 are covered in the scored cells 1-3; gene 3 is covered only in
    # cell 4, which the fit never scored, so it cannot be phased at all.
    ref_mat <- matrix(c(8L, 1L, 7L, 0L, 7L, 2L, 6L, 0L, 0L, 0L, 0L, 9L), nrow = 3, byrow = TRUE)
    alt_mat <- matrix(c(1L, 8L, 2L, 0L, 2L, 7L, 3L, 0L, 0L, 0L, 0L, 1L), nrow = 3, byrow = TRUE)
    post <- tibble::tibble(cell = 1:3, post_X1_active = c(0.99, 0.01, 0.98))

    rephased <- snplet:::.rephase_all_genes(ref_mat, alt_mat, post, rho = 0.05)

    # Confirm the guard drops the unscored cell instead of tripping the M-step
    # precondition, so genes backed by scored cells are still phased
    expect_false(any(is.na(rephased$h_g[1:2])))
    # Confirm a gene covered only in an unscored cell is reported unphaseable
    expect_true(is.na(rephased$h_g[3]))
    # Confirm its escape fraction is likewise NA rather than a bookkeeping filler
    expect_true(is.na(rephased$pi_g[3]))
})

test_that("assign_xci_by_clonotype recovers the clonal split and stores diagnostics", {
    fixture <- make_xci_snpdata()
    stored <- assign_xci_by_clonotype(fixture$snpdata, n_inits = 3)

    # Verify a stored SNPData object is returned
    expect_s4_class(stored, "SNPData")

    barcode_info <- barcode_info(stored)
    # Check the fit recorded that it ran on clonotype units
    expect_true(all(barcode_info$xci_fit_unit == "clonotype"))

    assignments <- xci_assignments(stored)
    # Verify every cell received an assignment row via the cell projection
    expect_equal(nrow(assignments), ncol(fixture$snpdata))

    # Confirm the inferred groups recover the true clonal split up to a label
    # swap (labels are exchangeable in the model)
    truth <- barcode_info$true_group
    called <- !is.na(assignments$active_x)
    agree <- mean(assignments$active_x[called] == truth[called])
    expect_true(max(agree, 1 - agree) > 0.9)
})

test_that("confidence_threshold controls how many cells are hard-assigned", {
    fixture <- make_xci_snpdata()

    strict <- assign_xci(fixture$snpdata, n_inits = 3, confidence_threshold = 0.999)
    lax <- assign_xci(fixture$snpdata, n_inits = 3, confidence_threshold = 0.6)

    n_called_strict <- sum(!is.na(barcode_info(strict)$active_x))
    n_called_lax <- sum(!is.na(barcode_info(lax)$active_x))

    # Verify a looser threshold assigns at least as many cells as a strict one
    expect_gte(n_called_lax, n_called_strict)

    # Confirm cells below the strict threshold receive NA rather than a call
    post <- barcode_info(strict)$xci_post_X1_active
    unassigned <- is.na(barcode_info(strict)$active_x)
    borderline <- post > 1 - 0.999 & post < 0.999
    # Every cell whose posterior sits inside the strict band must be unassigned
    expect_true(all(unassigned[borderline]))
})

test_that("assign_xci drops an escapee gene from the informative set but still phases it", {
    fixture <- make_xci_snpdata(escapee = TRUE)
    stored <- assign_xci(fixture$snpdata, n_inits = 3)

    # The clonal split must still be recovered
    assignments <- xci_assignments(stored)
    truth <- barcode_info(fixture$snpdata)$true_group
    called <- !is.na(assignments$active_x)
    agree <- mean(assignments$active_x[called] == truth[called])
    # Confirm agreement with the true clonal split exceeds 90%
    expect_true(max(agree, 1 - agree) > 0.9)

    # The escapee gene (last gene, balanced in every cell) carries no XCI signal
    # and must be filtered out of the informative (active-X-calling) set
    snp_info <- snp_info(stored)
    donor_snp_info <- donor_snp_info(stored)
    escapee_snp <- snp_info$snp_id[snp_info$gene_name == paste0("gene", 20)]
    escapee_informative <- donor_snp_info$xci_informative[donor_snp_info$snp_id == escapee_snp]
    # Confirm the escapee gene is not marked informative
    expect_false(any(escapee_informative))

    # Confirm genuine XCI genes are still retained
    expect_true(any(donor_snp_info$xci_informative))

    # Though excluded from calling, the escapee gene is still phased: its
    # balanced-in-every-cell signal should read as escape fraction near the
    # 0.499 ceiling, not be discarded as NA.
    escapee_escape_fraction <- donor_snp_info$xci_escape_fraction[donor_snp_info$snp_id == escapee_snp]
    # Confirm the escapee gene's escape fraction was estimated, not dropped
    expect_false(any(is.na(escapee_escape_fraction)))
    # Confirm it reads as strongly escaping (near-complete symmetry, not mild)
    expect_true(all(escapee_escape_fraction > 0.4))

    # xci_median_pi_g summarises only the informative gene set, so the
    # escapee's near-0.5 escape fraction must not pull the donor's recorded
    # background level up toward it.
    donor <- unique(donor_snp_info$donor)
    median_pi_g <- donor_info(stored)$xci_median_pi_g[donor_info(stored)$donor == donor]
    # Confirm the recorded background escape level excludes the escapee gene
    expect_true(median_pi_g < 0.4)
})

test_that("test_escape recommended pipeline flags an injected escapee gene against non-escapees", {
    # Full pipeline as documented in test_escape()'s examples: collapse
    # haplotype_expression()'s SNP x active-X-group rows to one row per
    # (donor, gene), attach the donor's empirical null, then BH-correct
    # per donor via group_modify.
    #
    # by_snp = TRUE is required here, not incidental: the injected escapee is
    # near-complete, so its dominant allele does not flip and no gene-level
    # election can pick it. The default grain would therefore never present it to
    # test_escape at all -- which is the escape-detection blind spot this test
    # exists to guard. Each gene here carries a single SNP, so summing the rows
    # pools only the two active-X groups and cannot over-count a molecule.
    fixture <- make_xci_snpdata(escapee = TRUE)
    stored <- assign_xci(fixture$snpdata, n_inits = 3)

    hap_by_gene <- haplotype_expression(stored, by_snp = TRUE) %>%
        dplyr::summarise(
            active_count = sum(active_count),
            inactive_count = sum(inactive_count),
            .by = c(donor, gene_name)
        ) %>%
        dplyr::left_join(dplyr::select(donor_info(stored), donor, xci_median_pi_g, xci_rho), by = "donor")

    escape_result <- hap_by_gene %>%
        dplyr::group_by(donor) %>%
        dplyr::group_modify(~ test_escape(.x, p = .x$xci_median_pi_g, rho = .x$xci_rho)) %>%
        dplyr::ungroup()

    # Verify the collapse actually produced one row per gene, not per SNP
    expect_equal(nrow(escape_result), length(unique(hap_by_gene$gene_name)))

    escapee_row <- escape_result[escape_result$gene_name == paste0("gene", 20), ]
    other_rows <- escape_result[escape_result$gene_name != paste0("gene", 20), ]
    # Confirm the injected escapee is flagged significant against the
    # per-donor empirical null
    expect_true(escapee_row$adj_p_val < 0.05)
    # Confirm the escapee's p-value is smaller than every genuine XCI gene's
    expect_true(escapee_row$p_val < min(other_rows$p_val))
})

test_that("haplotype_expression() surfaces an escapee gene that assign_xci excluded from calling", {
    fixture <- make_xci_snpdata(escapee = TRUE)
    stored <- assign_xci(fixture$snpdata, n_inits = 3)

    snp_info <- snp_info(stored)
    escapee_snp <- snp_info$snp_id[snp_info$gene_name == paste0("gene", 20)]

    # Verify the per-SNP grain reports the escapee at all, even though assign_xci
    # excluded it from calling and no gene-level election can pick it
    hap_default <- haplotype_expression(stored, by_snp = TRUE)
    expect_true(escapee_snp %in% hap_default$snp_id)
    # Confirm the escapee still reads as elevated inactive-haplotype expression
    # in at least one active-X group, even though assign_xci excluded it from
    # calling
    escapee_rows <- dplyr::filter(hap_default, snp_id == escapee_snp)
    expect_true(any(escapee_rows$escapes))

    # Both groups are reported, including one that reverses the phase. This gene
    # escapes nearly completely, so its two groups straddle 0.5 by noise and which
    # side each lands on depends on the EM's arbitrary X1/X2 labelling; dropping
    # the reversed group would let that coin flip decide whether a real escapee is
    # visible at all, so it is flagged instead.
    expect_equal(nrow(escapee_rows), 2L)
    expect_true(any(escapee_rows$phase_contradiction))
    # This escapee's two groups fall on OPPOSITE sides of 0.5 (one reversed, one
    # not), which is what breaks the flip and sets same_allele_dominant -- not the
    # magnitude of its escape. A gene escaping heavily but symmetrically would
    # still flip and still be elected. The consequence here is that this
    # particular escapee is absent from the default gene-level grain entirely.
    expect_true(all(escapee_rows$same_allele_dominant))
    expect_false(paste0("gene", 20) %in% haplotype_expression(stored)$gene_name)

    # Verify xci_informative_only = TRUE restricts back to the calling-informative set
    hap_informative_only <- haplotype_expression(stored, xci_informative_only = TRUE, by_snp = TRUE)
    expect_false(escapee_snp %in% hap_informative_only$snp_id)
})

test_that("assign_xci fits each donor independently in a multi-donor object", {
    fixture <- make_xci_snpdata(n_donors = 2)
    stored <- assign_xci(fixture$snpdata, n_inits = 3)

    barcode_info <- barcode_info(stored)
    # Verify both donors received assignments
    expect_setequal(unique(barcode_info$donor), c("donor0", "donor1"))

    # Confirm the clonal split is recovered within each donor (labels are
    # exchangeable per donor, so compare after resolving the swap per donor)
    agreement_rate <- function(donor) {
        rows <- barcode_info$donor == donor & !is.na(barcode_info$active_x)
        mean(barcode_info$active_x[rows] == barcode_info$true_group[rows])
    }
    agree_donor0 <- agreement_rate("donor0")
    agree_donor1 <- agreement_rate("donor1")

    # Confirm the clonal split is recovered within donor0
    expect_true(max(agree_donor0, 1 - agree_donor0) > 0.9)
    # Confirm the clonal split is recovered within donor1
    expect_true(max(agree_donor1, 1 - agree_donor1) > 0.9)
})

test_that("zygosity_source<- switches which stored zygosity call the het-SNP filter trusts", {
    fixture <- make_xci_snpdata(n_genes = 3, n_cells_per_group = 20, seed = 42)
    snp_info <- snp_info(fixture$snpdata)
    target_snp <- snp_info$snp_id[1]

    # fixture$snpdata already carries binomial-derived calls (infer_zygosity() inside
    # make_xci_snpdata()). Add a competing Vireo call for the same pair claiming it's
    # homozygous, contradicting what the binomial test genuinely found. The widened
    # (snp_id, donor, zygosity_source) key lets both coexist without conflict.
    snp_data <- add_donor_snp_metadata(
        fixture$snpdata,
        data.frame(snp_id = target_snp, donor = "donor0", zygosity = "hom", zygosity_source = "vireo_gt")
    )

    # Confirm the fixture's active source is still "binomial" (unaffected by adding a
    # second, non-active source)
    expect_equal(zygosity_source(snp_data), "binomial")
    trusting_binomial <- assign_xci(snp_data, n_inits = 3)
    informative_binomial <- donor_snp_info(trusting_binomial) %>%
        dplyr::filter(snp_id == target_snp, xci_informative) %>%
        nrow()
    # Verify the binomial call (correctly heterozygous) lets the SNP into the fit
    expect_equal(informative_binomial, 1)

    zygosity_source(snp_data) <- "vireo_gt"
    trusting_vireo <- assign_xci(snp_data, n_inits = 3)
    informative_vireo <- donor_snp_info(trusting_vireo) %>%
        dplyr::filter(snp_id == target_snp, xci_informative) %>%
        nrow()
    # Verify switching the active source to the (contradicting) Vireo "hom" call
    # excludes the SNP from the het-SNP filter
    expect_equal(informative_vireo, 0)
})

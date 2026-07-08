# ==============================================================================
# Test Suite: X-chromosome inactivation fitting
# Description: assign_inactive_x cell/clonotype modes, stored diagnostics, and
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
        group <- rep(c("X1", "X2"), each = n_cells_per_group)

        # Per gene, which allele sits on X1 (0 = REF active when X2 inactive).
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
                # X1-inactive cells silence the X1 allele; probability of REF
                x1_inactive <- group[c] == "X1"
                p_ref <- if (is_escapee[g]) {
                    0.5
                } else if (allele_on_x1[g] == 0) {
                    if (x1_inactive) 0.05 else 0.95
                } else {
                    if (x1_inactive) 0.95 else 0.05
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

        list(
            snpdata = SNPData(
                ref_count = ref_mat,
                alt_count = alt_mat,
                snp_info = snp_info,
                barcode_info = barcode_info
            ),
            allele_on_x1 = allele_on_x1
        )
    })
}

# ==============================================================================

test_that("assign_inactive_x recovers the true clonal split at cell level", {
    fixture <- make_xci_snpdata()
    stored <- assign_inactive_x(fixture$snpdata, n_inits = 3)

    assignments <- xci_assignments(stored)
    # Verify assignments carry the stored annotation columns
    expect_true(all(c("cell_id", "inactive_x", "xci_post_X1") %in% colnames(assignments)))
    # Verify every cell received an assignment row
    expect_equal(nrow(assignments), ncol(fixture$snpdata))

    # Confirm the two inferred groups recover the true clonal split (up to a
    # label swap, since X1/X2 labels are exchangeable in the model). Compare
    # only cells that met the confidence threshold.
    truth <- get_barcode_info(fixture$snpdata)$true_group
    called <- !is.na(assignments$inactive_x)
    agree <- mean(assignments$inactive_x[called] == truth[called])
    # Confirm agreement with the true clonal split exceeds 90% either way
    expect_true(max(agree, 1 - agree) > 0.9)
})

test_that("xci_haplotypes reports phase and escape fraction per informative SNP", {
    fixture <- make_xci_snpdata()
    stored <- assign_inactive_x(fixture$snpdata, n_inits = 3)

    haplotypes <- xci_haplotypes(stored)
    # Verify haplotype columns are present
    expect_true(all(c("snp_id", "gene_name", "allele_on_x1", "escape_fraction") %in% colnames(haplotypes)))
    # Check phase is reported as REF/ALT
    expect_true(all(haplotypes$allele_on_x1 %in% c("REF", "ALT")))
    # Confirm escape fraction is a valid minor fraction
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
    n_clonotypes <- dplyr::n_distinct(get_barcode_info(fixture$snpdata)$clonotype)
    # Confirm one assignment row exists per distinct clonotype
    expect_equal(nrow(fit$assignments), n_clonotypes)

    # Cells within a clonotype must share an assignment (clonal consistency)
    cell_assign <- fit$cell_assignments
    clono_map <- dplyr::select(get_barcode_info(fixture$snpdata), cell_id, clonotype)
    joined <- dplyr::inner_join(cell_assign, clono_map, by = "cell_id")
    per_clono <- tapply(joined$assignment, joined$clonotype, function(a) length(unique(a)))
    # Confirm every clonotype's cells were all projected to the same assignment
    expect_true(all(per_clono == 1))
})

test_that("assign_inactive_x promotes diagnostics into SNPData slots and survives subsetting", {
    fixture <- make_xci_snpdata()
    stored <- assign_inactive_x(fixture$snpdata, n_inits = 3)

    # Verify a SNPData object is returned
    expect_s4_class(stored, "SNPData")

    barcode_info <- get_barcode_info(stored)
    snp_info <- get_snp_info(stored)
    # Check barcode diagnostics were written
    expect_true(all(c("inactive_x", "xci_post_X1") %in% colnames(barcode_info)))
    # Check SNP diagnostics were written
    expect_true(all(c("xci_informative", "xci_allele_on_x1", "xci_escape_fraction") %in% colnames(snp_info)))
    # Confirm some SNPs are flagged informative and none are NA
    expect_true(any(snp_info$xci_informative))
    # Confirm the informative flag itself is never NA
    expect_false(any(is.na(snp_info$xci_informative)))

    # Diagnostics must survive cell subsetting because they live in barcode_info
    subset_cells <- stored[, 1:10]
    # Verify the inactive_x column is retained after subsetting cells
    expect_true("inactive_x" %in% colnames(get_barcode_info(subset_cells)))
    # Verify barcode_info row count matches the subset size
    expect_equal(nrow(get_barcode_info(subset_cells)), 10)

    # And survive SNP subsetting because they live in snp_info
    subset_snps <- filter_snps(stored, xci_informative)
    expect_true(all(get_snp_info(subset_snps)$xci_informative))
})

test_that("accessors and heatmap work on a stored SNPData object", {
    fixture <- make_xci_snpdata()
    stored <- assign_inactive_x(fixture$snpdata, n_inits = 3)

    assignments <- xci_assignments(stored)
    # Verify stored-object assignments expose the annotation columns
    expect_true(all(c("cell_id", "inactive_x", "xci_post_X1") %in% colnames(assignments)))

    haplotypes <- xci_haplotypes(stored)
    # Verify stored-object haplotypes expose phase and escape fraction
    expect_true(all(c("snp_id", "allele_on_x1", "escape_fraction") %in% colnames(haplotypes)))

    # Confirm the heatmap method runs on a stored object. It returns a drawn
    # HeatmapList (donor as title, unit-labelled column axis), so draw to a null
    # device to keep the test headless.
    grDevices::pdf(NULL)
    on.exit(grDevices::dev.off(), add = TRUE)
    hm <- plot_inactive_x_assignment_heatmap(stored, donor = "donor0")
    # Verify the heatmap method returns a drawn HeatmapList object
    expect_s4_class(hm, "HeatmapList")
})

# Shared stored-object fixture for heatmap display-parameter tests below
create_stored_xci_fixture <- function() {
    fixture <- make_xci_snpdata()
    assign_inactive_x(fixture$snpdata, n_inits = 3)
}

test_that("heatmap max_genes caps the number of displayed rows", {
    stored <- create_stored_xci_fixture()
    grDevices::pdf(NULL)
    on.exit(grDevices::dev.off(), add = TRUE)

    capped <- plot_inactive_x_assignment_heatmap(stored, donor = "donor0", max_genes = 5)

    # Confirm max_genes caps the number of displayed rows
    expect_equal(nrow(capped@ht_list[[1]]@matrix), 5)
})

test_that("heatmap max_genes above the available count shows all retained genes", {
    stored <- create_stored_xci_fixture()
    grDevices::pdf(NULL)
    on.exit(grDevices::dev.off(), add = TRUE)

    full <- plot_inactive_x_assignment_heatmap(stored, donor = "donor0")
    n_all <- nrow(full@ht_list[[1]]@matrix)
    big <- plot_inactive_x_assignment_heatmap(stored, donor = "donor0", max_genes = n_all + 100)

    # Confirm requesting more genes than available shows all retained genes
    expect_equal(nrow(big@ht_list[[1]]@matrix), n_all)
})

test_that("heatmap show_unassigned = FALSE drops unassigned columns", {
    stored <- create_stored_xci_fixture()
    grDevices::pdf(NULL)
    on.exit(grDevices::dev.off(), add = TRUE)

    assigned_only <- plot_inactive_x_assignment_heatmap(
        stored,
        donor = "donor0",
        show_unassigned = FALSE
    )
    n_assigned <- sum(!is.na(get_barcode_info(stored)$inactive_x))

    # Check show_unassigned = FALSE drops unassigned columns
    expect_equal(ncol(assigned_only@ht_list[[1]]@matrix), n_assigned)
})

test_that("heatmap accepts toggling gene names and row clustering", {
    stored <- create_stored_xci_fixture()
    grDevices::pdf(NULL)
    on.exit(grDevices::dev.off(), add = TRUE)

    variant <- plot_inactive_x_assignment_heatmap(
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

    no_marks <- plot_inactive_x_assignment_heatmap(
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

    with_post <- plot_inactive_x_assignment_heatmap(stored, donor = "donor0")
    without_post <- plot_inactive_x_assignment_heatmap(
        stored,
        donor = "donor0",
        show_posterior = FALSE
    )
    with_names <- names(with_post@ht_list[[1]]@top_annotation@anno_list)
    without_names <- names(without_post@ht_list[[1]]@top_annotation@anno_list)

    # Confirm the posterior_X1 annotation is present by default
    expect_true("posterior_X1" %in% with_names)
    # Confirm the posterior_X1 annotation is dropped when show_posterior = FALSE
    expect_false("posterior_X1" %in% without_names)
    # Confirm the assignment annotation is retained either way
    expect_true("assignment" %in% without_names)
})

test_that("heatmap applies custom colour arguments to the plot body", {
    stored <- create_stored_xci_fixture()
    grDevices::pdf(NULL)
    on.exit(grDevices::dev.off(), add = TRUE)

    coloured <- plot_inactive_x_assignment_heatmap(
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
    stored <- assign_inactive_x(fixture$snpdata, n_inits = 3)

    # Force one cell to look like it was never scored (NA posterior): the model
    # gives it no coverage, distinct from a low-confidence unassigned call.
    bi <- get_barcode_info(stored)
    bi$inactive_x[1] <- NA
    bi$xci_post_X1[1] <- NA_real_
    stored <- add_barcode_metadata(
        stored,
        dplyr::select(bi, cell_id, inactive_x, xci_post_X1),
        join_by = "cell_id",
        overwrite = TRUE
    )

    grDevices::pdf(NULL)
    on.exit(grDevices::dev.off(), add = TRUE)

    # With show_no_coverage = TRUE the annotation carries a distinct level
    hm <- plot_inactive_x_assignment_heatmap(
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
    default_hm <- plot_inactive_x_assignment_heatmap(stored, donor = "donor0")
    n_donor0 <- sum(get_barcode_info(stored)$donor == "donor0")
    expect_equal(ncol(default_hm@ht_list[[1]]@matrix), n_donor0 - 1)
})

test_that("clonotype fit requires clonotype information", {
    fixture <- make_xci_snpdata()
    # Strip clonotype information
    no_clono <- add_barcode_metadata(
        fixture$snpdata,
        data.frame(
            cell_id = get_barcode_info(fixture$snpdata)$cell_id,
            clonotype = NA_character_
        ),
        overwrite = TRUE
    )
    # Verify a clear error is raised when clonotypes are all NA
    expect_error(
        assign_inactive_x_by_clonotype(no_clono),
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

test_that("assign_inactive_x_by_clonotype recovers the clonal split and stores diagnostics", {
    fixture <- make_xci_snpdata()
    stored <- assign_inactive_x_by_clonotype(fixture$snpdata, n_inits = 3)

    # Verify a stored SNPData object is returned
    expect_s4_class(stored, "SNPData")

    barcode_info <- get_barcode_info(stored)
    # Check the fit recorded that it ran on clonotype units
    expect_true(all(barcode_info$xci_fit_unit == "clonotype"))

    assignments <- xci_assignments(stored)
    # Verify every cell received an assignment row via the cell projection
    expect_equal(nrow(assignments), ncol(fixture$snpdata))

    # Confirm the inferred groups recover the true clonal split up to a label
    # swap (labels are exchangeable in the model)
    truth <- barcode_info$true_group
    called <- !is.na(assignments$inactive_x)
    agree <- mean(assignments$inactive_x[called] == truth[called])
    expect_true(max(agree, 1 - agree) > 0.9)
})

test_that("confidence_threshold controls how many cells are hard-assigned", {
    fixture <- make_xci_snpdata()

    strict <- assign_inactive_x(fixture$snpdata, n_inits = 3, confidence_threshold = 0.999)
    lax <- assign_inactive_x(fixture$snpdata, n_inits = 3, confidence_threshold = 0.6)

    n_called_strict <- sum(!is.na(get_barcode_info(strict)$inactive_x))
    n_called_lax <- sum(!is.na(get_barcode_info(lax)$inactive_x))

    # Verify a looser threshold assigns at least as many cells as a strict one
    expect_gte(n_called_lax, n_called_strict)

    # Confirm cells below the strict threshold receive NA rather than a call
    post <- get_barcode_info(strict)$xci_post_X1
    unassigned <- is.na(get_barcode_info(strict)$inactive_x)
    borderline <- post > 1 - 0.999 & post < 0.999
    # Every cell whose posterior sits inside the strict band must be unassigned
    expect_true(all(unassigned[borderline]))
})

test_that("refit_after_filter returns a valid fit and drops escapee genes", {
    fixture <- make_xci_snpdata(escapee = TRUE)

    # Verify the refit path returns a valid stored object
    stored <- assign_inactive_x(fixture$snpdata, n_inits = 3, refit_after_filter = TRUE)
    expect_s4_class(stored, "SNPData")

    # The clonal split must still be recovered after refitting
    assignments <- xci_assignments(stored)
    truth <- get_barcode_info(fixture$snpdata)$true_group
    called <- !is.na(assignments$inactive_x)
    agree <- mean(assignments$inactive_x[called] == truth[called])
    # Confirm agreement with the true clonal split still exceeds 90% after refitting
    expect_true(max(agree, 1 - agree) > 0.9)

    # The escapee gene (last gene, balanced in every cell) carries no XCI signal
    # and must be filtered out of the informative set by the refit pass
    snp_info <- get_snp_info(stored)
    escapee_snp <- snp_info$snp_id[snp_info$gene_name == paste0("gene", 20)]
    escapee_informative <- snp_info$xci_informative[snp_info$snp_id == escapee_snp]
    # Confirm the escapee gene is not marked informative
    expect_false(isTRUE(escapee_informative))

    # Confirm genuine XCI genes are still retained
    expect_true(any(snp_info$xci_informative))
})

test_that("assign_inactive_x fits each donor independently in a multi-donor object", {
    fixture <- make_xci_snpdata(n_donors = 2)
    stored <- assign_inactive_x(fixture$snpdata, n_inits = 3)

    barcode_info <- get_barcode_info(stored)
    # Verify both donors received assignments
    expect_setequal(unique(barcode_info$donor), c("donor0", "donor1"))

    # Confirm the clonal split is recovered within each donor (labels are
    # exchangeable per donor, so compare after resolving the swap per donor)
    agreement_rate <- function(donor) {
        rows <- barcode_info$donor == donor & !is.na(barcode_info$inactive_x)
        mean(barcode_info$inactive_x[rows] == barcode_info$true_group[rows])
    }
    agree_donor0 <- agreement_rate("donor0")
    agree_donor1 <- agreement_rate("donor1")

    # Confirm the clonal split is recovered within donor0
    expect_true(max(agree_donor0, 1 - agree_donor0) > 0.9)
    # Confirm the clonal split is recovered within donor1
    expect_true(max(agree_donor1, 1 - agree_donor1) > 0.9)
})

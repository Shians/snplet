# ==============================================================================
# Test Suite: X-chromosome inactivation fitting
# Description: fit_inactive_x cell/clonotype modes, stored diagnostics, and
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
    seed = 1
) {
    withr::with_seed(seed, {
        n_cells <- 2 * n_cells_per_group
        group <- rep(c("X1", "X2"), each = n_cells_per_group)

        # Per gene, which allele sits on X1 (0 = REF active when X2 inactive)
        allele_on_x1 <- sample(0:1, n_genes, replace = TRUE)

        ref_mat <- matrix(0L, nrow = n_genes, ncol = n_cells)
        alt_mat <- matrix(0L, nrow = n_genes, ncol = n_cells)
        for (g in seq_len(n_genes)) {
            for (c in seq_len(n_cells)) {
                # X1-inactive cells silence the X1 allele; probability of REF
                x1_inactive <- group[c] == "X1"
                p_ref <- if (allele_on_x1[g] == 0) {
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
        barcode_info <- data.frame(
            barcode = paste0("cell", seq_len(n_cells)),
            donor = "donor0",
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

test_that("fit_inactive_x returns an xci_fit with cell-level assignments", {
    fixture <- make_xci_snpdata()
    fit <- fit_inactive_x(fixture$snpdata, n_inits = 3)

    # Verify the fit is an xci_fit object keyed by donor
    expect_s3_class(fit, "xci_fit")
    expect_true("donor0" %in% names(fit))

    assignments <- xci_assignments(fit)
    # Verify assignments carry the expected columns
    expect_true(all(c("donor", "cell_id", "post_X1", "post_X2", "assignment") %in% colnames(assignments)))
    # Verify every cell received an assignment row
    expect_equal(nrow(assignments), ncol(fixture$snpdata))

    # Confirm the two inferred groups recover the true clonal split (up to a
    # label swap, since X1/X2 labels are exchangeable in the model).
    truth <- get_barcode_info(fixture$snpdata)$true_group
    agree <- mean(assignments$assignment == truth)
    expect_true(max(agree, 1 - agree) > 0.9)
})

test_that("xci_haplotypes reports phase and escape fraction per informative SNP", {
    fixture <- make_xci_snpdata()
    fit <- fit_inactive_x(fixture$snpdata, n_inits = 3)

    haplotypes <- xci_haplotypes(fit)
    # Verify haplotype columns are present
    expect_true(all(c("donor", "snp_id", "gene_name", "allele_on_x1", "escape_fraction") %in% colnames(haplotypes)))
    # Check phase is reported as REF/ALT
    expect_true(all(haplotypes$allele_on_x1 %in% c("REF", "ALT")))
    # Confirm escape fraction is a valid minor fraction
    expect_true(all(haplotypes$escape_fraction > 0 & haplotypes$escape_fraction < 0.5))
})

test_that("fit_inactive_x by clonotype projects assignments back to cells", {
    fixture <- make_xci_snpdata()
    fit <- fit_inactive_x(fixture$snpdata, n_inits = 3, by = "clonotype")

    # Verify the clonotype fit records its unit and carries a cell projection
    expect_equal(fit[["donor0"]]$unit, "clonotype")
    expect_true(!is.null(fit[["donor0"]]$cell_assignments))

    # The EM ran per clonotype, one assignment row each
    n_clonotypes <- dplyr::n_distinct(get_barcode_info(fixture$snpdata)$clonotype)
    expect_equal(nrow(fit[["donor0"]]$assignments), n_clonotypes)

    # Cells within a clonotype must share an assignment (clonal consistency)
    cell_assign <- fit[["donor0"]]$cell_assignments
    clono_map <- dplyr::select(get_barcode_info(fixture$snpdata), cell_id, clonotype)
    joined <- dplyr::inner_join(cell_assign, clono_map, by = "cell_id")
    per_clono <- tapply(joined$assignment, joined$clonotype, function(a) length(unique(a)))
    expect_true(all(per_clono == 1))
})

test_that("store = TRUE promotes diagnostics into SNPData slots and survives subsetting", {
    fixture <- make_xci_snpdata()
    stored <- fit_inactive_x(fixture$snpdata, n_inits = 3, store = TRUE)

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
    expect_false(any(is.na(snp_info$xci_informative)))

    # Diagnostics must survive cell subsetting because they live in barcode_info
    subset_cells <- stored[, 1:10]
    expect_true("inactive_x" %in% colnames(get_barcode_info(subset_cells)))
    expect_equal(nrow(get_barcode_info(subset_cells)), 10)

    # And survive SNP subsetting because they live in snp_info
    subset_snps <- filter_snps(stored, xci_informative)
    expect_true(all(get_snp_info(subset_snps)$xci_informative))
})

test_that("accessors and heatmap work on a stored SNPData object", {
    fixture <- make_xci_snpdata()
    stored <- fit_inactive_x(fixture$snpdata, n_inits = 3, store = TRUE)

    assignments <- xci_assignments(stored)
    # Verify stored-object assignments expose the annotation columns
    expect_true(all(c("cell_id", "inactive_x", "xci_post_X1") %in% colnames(assignments)))

    haplotypes <- xci_haplotypes(stored)
    # Verify stored-object haplotypes expose phase and escape fraction
    expect_true(all(c("snp_id", "allele_on_x1", "escape_fraction") %in% colnames(haplotypes)))

    # Confirm the heatmap method runs on a stored object
    hm <- plot_inactive_x_assignment_heatmap(stored, donor = "donor0")
    expect_s4_class(hm, "Heatmap")
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
        fit_inactive_x(no_clono, by = "clonotype"),
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
    expect_identical(both$L1, direct_L1)
})

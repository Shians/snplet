# ==============================================================================
# Test Suite: Phased haplotype expression
# Description: haplotype_expression() splits ALT/REF into active/inactive
#              haplotype reads per donor using stored XCI phase, flags the
#              same-allele-dominant paradox, and honours per-donor phase.
# ==============================================================================

library(testthat)
library(Matrix)

# ------------------------------------------------------------------------------
# Test Data Setup
# ------------------------------------------------------------------------------

# Build a SNPData with XCI diagnostics injected directly (rather than fitting the
# EM) so the counts and stored phase are fully controlled and deterministic.
#
# 8 cells: two donors, each with two X2-active and two X1-active cells.
# Columns:  d0_x2a d0_x2a d0_x1a d0_x1a | d1_x2a d1_x2a d1_x1a d1_x1a
#
# Three SNPs, each depth 10 per cell:
#   snp1 - clean XCI in donor0 only. X1 carries REF, so X2-active cells express
#          ALT and X1-active cells express REF: the dominant physical allele
#          flips between groups.
#   snp2 - escaper in donor0. REF dominates in BOTH groups (biologically
#          impossible under clean XCI), with elevated inactive-haplotype reads.
#   snp3 - informative in BOTH donors with OPPOSITE phase (X1=REF in donor0,
#          X1=ALT in donor1). The scalar phase column stores donor0's value, so a
#          fit that flattened phase would mis-split donor1; clean XCI counts in
#          both donors let the test catch that.
make_hap_fixture <- function() {
    # rows = snp1, snp2, snp3; cols = the 8 cells above
    ref <- rbind(
        c(0L, 0L, 10L, 10L, 0L, 0L, 0L, 0L),
        c(7L, 7L, 8L, 8L, 0L, 0L, 0L, 0L),
        c(0L, 0L, 10L, 10L, 10L, 10L, 0L, 0L)
    )
    alt <- rbind(
        c(10L, 10L, 0L, 0L, 0L, 0L, 0L, 0L),
        c(3L, 3L, 2L, 2L, 0L, 0L, 0L, 0L),
        c(10L, 10L, 0L, 0L, 0L, 0L, 10L, 10L)
    )

    snp_info <- data.frame(
        chrom = "X",
        pos = c(1000L, 2000L, 3000L),
        ref = "A",
        alt = "G",
        gene_name = c("gene1", "gene2", "gene3"),
        stringsAsFactors = FALSE
    )
    barcode_info <- data.frame(
        barcode = paste0("cell", 1:8),
        donor = rep(c("donor0", "donor1"), each = 4),
        stringsAsFactors = FALSE
    )

    obj <- SNPData(
        ref_count = Matrix(ref, sparse = TRUE),
        alt_count = Matrix(alt, sparse = TRUE),
        snp_info = snp_info,
        barcode_info = barcode_info
    )

    snp_ids <- get_snp_info(obj)$snp_id
    cell_ids <- get_barcode_info(obj)$cell_id

    # One row per informative (snp_id, donor) pair: snp1 and snp2 are
    # informative in donor0 only; snp3 is informative in both donors, with
    # opposite phase.
    donor_snp_diag <- data.frame(
        snp_id = c(snp_ids[1], snp_ids[2], snp_ids[3], snp_ids[3]),
        donor = c("donor0", "donor0", "donor0", "donor1"),
        xci_informative = TRUE,
        allele_on_x1 = c("REF", "REF", "REF", "ALT"),
        stringsAsFactors = FALSE
    )
    barcode_diag <- data.frame(
        cell_id = cell_ids,
        active_x = c("X2", "X2", "X1", "X1", "X2", "X2", "X1", "X1"),
        stringsAsFactors = FALSE
    )

    obj <- add_donor_snp_metadata(obj, donor_snp_diag, join_by = c("snp_id", "donor"), overwrite = TRUE)
    obj <- add_barcode_metadata(obj, barcode_diag, join_by = "cell_id", overwrite = TRUE)

    list(obj = obj, snp_ids = snp_ids)
}

# Pull the single row for one donor / SNP / active-X group.
hap_row <- function(res, donor, snp_id, active_x) {
    dplyr::filter(res, donor == !!donor, snp_id == !!snp_id, active_x == !!active_x)
}

# ==============================================================================

test_that("haplotype_expression() errors when no XCI diagnostics are stored", {
    ref <- Matrix(matrix(1L, 2, 2), sparse = TRUE)
    alt <- Matrix(matrix(1L, 2, 2), sparse = TRUE)
    snp_info <- data.frame(chrom = "X", pos = c(1L, 2L), ref = "A", alt = "G")
    barcode_info <- data.frame(barcode = c("c1", "c2"))
    obj <- SNPData(ref_count = ref, alt_count = alt, snp_info = snp_info, barcode_info = barcode_info)

    # Verify the method refuses to run without stored phase/assignments
    expect_error(haplotype_expression(obj), "Run assign_xci")
})

test_that("haplotype_expression() carries a donor column keyed to the assignments", {
    fixture <- make_hap_fixture()

    res <- haplotype_expression(fixture$obj)

    # Verify the donor column is present
    expect_true("donor" %in% colnames(res))
    # Confirm both donors appear in the output
    expect_setequal(res$donor, c("donor0", "donor1"))
})

test_that("haplotype_expression() reports a clean XCI SNP as flipping alleles", {
    fixture <- make_hap_fixture()
    snp1 <- fixture$snp_ids[1]

    res <- haplotype_expression(fixture$obj)
    x1 <- hap_row(res, "donor0", snp1, "X1")
    x2 <- hap_row(res, "donor0", snp1, "X2")

    # In X1-active cells the active X1 haplotype (REF) dominates
    expect_equal(x1$dominant_allele, "REF")
    # In X2-active cells the active X2 haplotype (ALT) dominates
    expect_equal(x2$dominant_allele, "ALT")
    # Inactive-haplotype reads are absent, so escape fraction is zero both ways
    expect_equal(x1$escape_fraction, 0)
    expect_equal(x2$escape_fraction, 0)
    # Confirm the dominant allele flips, so the SNP is not flagged paradoxical
    expect_false(unique(c(x1$same_allele_dominant, x2$same_allele_dominant)))
})

test_that("haplotype_expression() flags an escaper where the same allele dominates both groups", {
    fixture <- make_hap_fixture()
    snp2 <- fixture$snp_ids[2]

    res <- haplotype_expression(fixture$obj)
    x1 <- hap_row(res, "donor0", snp2, "X1")
    x2 <- hap_row(res, "donor0", snp2, "X2")

    # REF dominates in both active-X groups
    expect_equal(x1$dominant_allele, "REF")
    expect_equal(x2$dominant_allele, "REF")
    # Elevated inactive-haplotype fraction quantifies the escape. X1 carries REF,
    # so the X1-active group reads 4/20 from the silenced X2 and the X2-active
    # group reads 14/20 from the silenced X1.
    expect_equal(x1$escape_fraction, 0.2)
    expect_equal(x2$escape_fraction, 0.7)
    # Confirm the biologically impossible same-allele dominance is flagged
    expect_true(x1$same_allele_dominant)
    expect_true(x2$same_allele_dominant)
})

test_that("haplotype_expression() splits active/inactive counts by the stored phase", {
    fixture <- make_hap_fixture()
    snp2 <- fixture$snp_ids[2]

    res <- haplotype_expression(fixture$obj)
    x2 <- hap_row(res, "donor0", snp2, "X2")

    # X1 carries REF; in X2-active cells REF is the silenced (inactive) allele:
    # pooled ref = 14 (inactive), pooled alt = 6 (active X2 haplotype)
    expect_equal(x2$inactive_count, 14)
    expect_equal(x2$active_count, 6)
    expect_equal(x2$coverage, 20)
})

test_that("haplotype_expression() applies each donor's own phase, not the flattened scalar", {
    fixture <- make_hap_fixture()
    snp3 <- fixture$snp_ids[3]

    res <- haplotype_expression(fixture$obj)
    d1_x1 <- hap_row(res, "donor1", snp3, "X1")
    d1_x2 <- hap_row(res, "donor1", snp3, "X2")

    # snp3 is clean XCI in donor1 under its own phase (X1 = ALT). Had the scalar
    # donor0 phase (REF) been used instead, the split would invert and escape
    # would read ~1 rather than ~0.
    expect_equal(d1_x1$escape_fraction, 0)
    expect_equal(d1_x2$escape_fraction, 0)
    # Confirm the dominant physical allele still flips between groups in donor1
    expect_equal(d1_x1$dominant_allele, "ALT")
    expect_equal(d1_x2$dominant_allele, "REF")
})

test_that("haplotype_expression() omits SNPs a donor's model never retained", {
    fixture <- make_hap_fixture()
    snp1 <- fixture$snp_ids[1]

    res <- haplotype_expression(fixture$obj)

    # snp1 is informative in donor0 only, so no donor1 rows should exist for it
    expect_equal(nrow(dplyr::filter(res, donor == "donor1", snp_id == snp1)), 0L)
})

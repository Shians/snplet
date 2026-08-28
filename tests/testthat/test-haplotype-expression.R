# ==============================================================================
# Test Suite: Phased haplotype expression
# Description: haplotype_expression() splits ALT/REF into active/inactive
#              haplotype reads per donor using stored XCI phase, elects one
#              representative SNP per gene by default, flags the
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

    snp_ids <- snp_info(obj)$snp_id
    cell_ids <- barcode_info(obj)$cell_id

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

    res <- haplotype_expression(fixture$obj, by_snp = TRUE)

    # Verify the donor column is present
    expect_true("donor" %in% colnames(res))
    # Confirm both donors appear in the output
    expect_setequal(res$donor, c("donor0", "donor1"))
})

test_that("haplotype_expression() reports a clean XCI SNP as flipping alleles", {
    fixture <- make_hap_fixture()
    snp1 <- fixture$snp_ids[1]

    res <- haplotype_expression(fixture$obj, by_snp = TRUE, by_active_x = TRUE)
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

test_that("haplotype_expression() flags a group that reverses the phase without dropping it", {
    fixture <- make_hap_fixture()
    snp2 <- fixture$snp_ids[2]

    res <- haplotype_expression(fixture$obj, by_snp = TRUE, by_active_x = TRUE)
    x1 <- hap_row(res, "donor0", snp2, "X1")
    x2 <- hap_row(res, "donor0", snp2, "X2")

    # snp2's X2-active group has its inactive (X1) haplotype outnumbering its
    # active (X2) one: escape_fraction 0.7, a phase contradiction rather than
    # partial leakage
    expect_equal(x2$escape_fraction, 0.7)
    expect_true(x2$phase_contradiction)
    # Verify the group is reported rather than excluded, so a near-complete
    # escapee cannot be hidden by which side of 0.5 the EM's arbitrary X1/X2
    # labelling happens to place it
    expect_equal(nrow(x2), 1L)
    # Confirm the partner group is unflagged: at 0.2 it is ordinary leakage
    expect_equal(x1$escape_fraction, 0.2)
    expect_false(x1$phase_contradiction)
})

test_that("same_allele_dominant tracks the groups straddling 0.5, not high escape", {
    # One flagged group does not imply the flip failed: what breaks the flip is the
    # two groups disagreeing about which side of 0.5 they sit on. When BOTH exceed
    # 0.5 each is dominated by the other X's allele -- different physical alleles --
    # so the SNP still flips and stays electable. A fixture with only a one-sided
    # case would pass a wrong assertion here, so both cases are built explicitly.
    #
    # The symmetric case is constructed, not observed: assign_xci() would have
    # inverted such a SNP's phase (escape > 0.5 is unidentifiable from expression
    # alone, so its phase step always makes the silenced allele the minority).
    # The phase is injected by hand here precisely to pin the branch, which real
    # fits never reach -- so a future change cannot quietly redefine the flip test
    # as a magnitude test and still pass.
    #
    # X1 carries REF; 2 cells per group, depth 100.
    #   one_sided - escape 0.20 in the X1-active group, 0.80 in the X2-active one
    #   symmetric - escape 0.80 in both
    ref <- rbind(
        c(80L, 80L, 80L, 80L),
        c(20L, 20L, 80L, 80L)
    )
    alt <- rbind(
        c(20L, 20L, 20L, 20L),
        c(80L, 80L, 20L, 20L)
    )
    snp_info <- data.frame(
        chrom = "X",
        pos = c(1000L, 2000L),
        ref = "A",
        alt = "G",
        gene_name = c("one_sided", "symmetric"),
        stringsAsFactors = FALSE
    )
    barcode_info <- data.frame(barcode = paste0("cell", 1:4), donor = "donor0", stringsAsFactors = FALSE)
    obj <- SNPData(
        ref_count = Matrix(ref, sparse = TRUE),
        alt_count = Matrix(alt, sparse = TRUE),
        snp_info = snp_info,
        barcode_info = barcode_info
    )
    obj <- add_donor_snp_metadata(
        obj,
        data.frame(
            snp_id = snp_info(obj)$snp_id,
            donor = "donor0",
            xci_informative = TRUE,
            allele_on_x1 = "REF",
            stringsAsFactors = FALSE
        ),
        join_by = c("snp_id", "donor"),
        overwrite = TRUE
    )
    obj <- add_barcode_metadata(
        obj,
        data.frame(
            cell_id = barcode_info(obj)$cell_id,
            active_x = c("X1", "X1", "X2", "X2"),
            stringsAsFactors = FALSE
        ),
        join_by = "cell_id",
        overwrite = TRUE
    )

    res <- haplotype_expression(obj, by_snp = TRUE, by_active_x = TRUE)
    one_sided <- dplyr::filter(res, gene_name == "one_sided")
    symmetric <- dplyr::filter(res, gene_name == "symmetric")

    # Verify only one of one_sided's groups crosses 0.5, and the flip fails
    expect_equal(sum(one_sided$phase_contradiction), 1L)
    expect_true(all(one_sided$same_allele_dominant))

    # Verify both of symmetric's groups cross 0.5, yet the dominant allele still
    # flips -- high escape alone is not a phase contradiction
    expect_equal(sum(symmetric$phase_contradiction), 2L)
    expect_false(any(symmetric$same_allele_dominant))
    expect_setequal(symmetric$dominant_allele, c("ALT", "REF"))

    # Confirm the election follows same_allele_dominant, so the symmetric escapee
    # is still elected at gene level while the one-sided SNP is dropped
    elected <- haplotype_expression(obj)
    expect_setequal(elected$gene_name, "symmetric")
})


test_that("phase_contradiction is mechanical, with no per-gene exemption", {
    fixture <- make_hap_fixture()
    snp2 <- fixture$snp_ids[2]
    snp_info(fixture$obj)$gene_name[snp_info(fixture$obj)$snp_id == snp2] <- "XIST"

    res <- haplotype_expression(fixture$obj, by_snp = TRUE, by_active_x = TRUE)
    x2 <- hap_row(res, "donor0", snp2, "X2")

    # The reversed group is reported rather than excluded, as for any gene
    expect_equal(nrow(x2), 1L)
    expect_equal(x2$escape_fraction, 0.7)
    # A former carve-out suppressed this flag for XIST alone. It never fired for
    # the case it was written for -- assign_xci() cannot store a phase that puts
    # XIST above 0.5 -- and where it could fire it made two genes with identical
    # counts report different flags. The flag is now purely a function of the
    # counts.
    expect_true(x2$phase_contradiction)
})

test_that("phase_likely_inverted labels curated genes without altering their counts", {
    fixture <- make_hap_fixture()
    snp1 <- fixture$snp_ids[1]
    snp_info(fixture$obj)$gene_name[snp_info(fixture$obj)$snp_id == snp1] <- "XIST"

    res <- haplotype_expression(fixture$obj, by_snp = TRUE)
    plain <- haplotype_expression(fixture$obj, by_snp = TRUE, inverted_phase_genes = character(0))

    # Verify the curated gene is labelled and other genes are not
    expect_true(all(dplyr::filter(res, gene_name == "XIST")$phase_likely_inverted))
    expect_false(any(dplyr::filter(res, gene_name != "XIST")$phase_likely_inverted))

    # Ensure the label is inert: counts and escape fractions are untouched, so
    # the flag asserts prior biology rather than correcting anything measured
    expect_equal(res$active_count, plain$active_count)
    expect_equal(res$inactive_count, plain$inactive_count)
    expect_equal(res$escape_fraction, plain$escape_fraction)
    # Confirm the list is configurable rather than hard-coded
    expect_false(any(plain$phase_likely_inverted))
})

test_that("phase_likely_inverted survives the gene-level election", {
    fixture <- make_hap_fixture()
    snp1 <- fixture$snp_ids[1]
    snp_info(fixture$obj)$gene_name[snp_info(fixture$obj)$snp_id == snp1] <- "XIST"

    res <- haplotype_expression(fixture$obj)

    # The flag is a property of the gene, not a selection criterion, so unlike
    # same_allele_dominant it must be carried into the elected output
    expect_true("phase_likely_inverted" %in% colnames(res))
    expect_true(all(dplyr::filter(res, gene_name == "XIST")$phase_likely_inverted))
})

test_that("haplotype_expression() splits active/inactive counts by the stored phase", {
    fixture <- make_hap_fixture()
    snp2 <- fixture$snp_ids[2]
    # snp2 carries partial leakage (16:4), which is what makes it a useful test of
    # where the minority reads land -- snp1's 20:0 would pass even if leakage were
    # discarded entirely.
    res <- haplotype_expression(fixture$obj, by_snp = TRUE, by_active_x = TRUE)
    x1 <- hap_row(res, "donor0", snp2, "X1")

    # X1 carries REF; in X1-active cells REF is the active haplotype:
    # pooled ref = 16 (active), pooled alt = 4 (inactive X2 haplotype)
    expect_equal(x1$active_count, 16)
    expect_equal(x1$inactive_count, 4)
    expect_equal(x1$coverage, 20)
})

test_that("haplotype_expression() applies each donor's own phase, not the flattened scalar", {
    fixture <- make_hap_fixture()
    snp3 <- fixture$snp_ids[3]

    res <- haplotype_expression(fixture$obj, by_snp = TRUE, by_active_x = TRUE)
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

    res <- haplotype_expression(fixture$obj, by_snp = TRUE)

    # snp1 is informative in donor0 only, so no donor1 rows should exist for it
    expect_equal(nrow(dplyr::filter(res, donor == "donor1", snp_id == snp1)), 0L)
})

# ------------------------------------------------------------------------------
# Gene-level election (the default grain)
# ------------------------------------------------------------------------------

# A fixture with two SNPs per gene, so the election between them is exercised.
# The existing fixture has one SNP per gene and cannot distinguish "picked the
# right SNP" from "picked the only SNP".
#
# 4 cells, one donor: cells 1-2 X1-active, cells 3-4 X2-active. X1 carries REF
# at every SNP.
#   geneA/snpA_hi  - REF dominant in BOTH groups (tied 5:5 in the X2 group, so
#                    its escape fraction is exactly 0.5 and it survives the
#                    phase-contradiction filter). Coverage 40: the best-covered
#                    SNP of the gene, and the one a naive coverage-only rule
#                    would elect.
#   geneA/snpA_lo  - clean flip, coverage 12. Must win despite lower coverage.
#   geneB/snpB_low - clean flip, coverage 12.
#   geneB/snpB_hi  - clean flip, coverage 20. Must win on coverage.
#   geneC/snpC     - REF dominant in both groups, the gene's only SNP, so geneC
#                    has no valid representative and must be dropped.
make_gene_collapse_fixture <- function() {
    ref <- rbind(
        c(8L, 8L, 5L, 5L),
        c(3L, 3L, 0L, 0L),
        c(3L, 3L, 0L, 0L),
        c(5L, 5L, 0L, 0L),
        c(8L, 8L, 5L, 5L)
    )
    alt <- rbind(
        c(2L, 2L, 5L, 5L),
        c(0L, 0L, 3L, 3L),
        c(0L, 0L, 3L, 3L),
        c(0L, 0L, 5L, 5L),
        c(2L, 2L, 5L, 5L)
    )

    snp_info <- data.frame(
        chrom = "X",
        pos = c(1000L, 2000L, 3000L, 4000L, 5000L),
        ref = "A",
        alt = "G",
        gene_name = c("geneA", "geneA", "geneB", "geneB", "geneC"),
        stringsAsFactors = FALSE
    )
    barcode_info <- data.frame(
        barcode = paste0("cell", 1:4),
        donor = "donor0",
        stringsAsFactors = FALSE
    )

    obj <- SNPData(
        ref_count = Matrix(ref, sparse = TRUE),
        alt_count = Matrix(alt, sparse = TRUE),
        snp_info = snp_info,
        barcode_info = barcode_info
    )

    snp_ids <- snp_info(obj)$snp_id
    donor_snp_diag <- data.frame(
        snp_id = snp_ids,
        donor = "donor0",
        xci_informative = TRUE,
        allele_on_x1 = "REF",
        stringsAsFactors = FALSE
    )
    barcode_diag <- data.frame(
        cell_id = barcode_info(obj)$cell_id,
        active_x = c("X1", "X1", "X2", "X2"),
        stringsAsFactors = FALSE
    )

    obj <- add_donor_snp_metadata(obj, donor_snp_diag, join_by = c("snp_id", "donor"), overwrite = TRUE)
    obj <- add_barcode_metadata(obj, barcode_diag, join_by = "cell_id", overwrite = TRUE)

    list(obj = obj, snp_ids = snp_ids)
}

test_that("haplotype_expression() elects the flipping SNP over a better-covered non-flipping one", {
    fixture <- make_gene_collapse_fixture()

    res <- haplotype_expression(fixture$obj)
    gene_a <- dplyr::filter(res, gene_name == "geneA")

    # Verify geneA is represented by the flipping SNP, not the higher-coverage
    # SNP whose dominant allele fails to flip between the groups
    expect_setequal(gene_a$snp_id, fixture$snp_ids[2])
    # Confirm the gene collapses to a single row, its two active-X groups pooled
    expect_equal(nrow(gene_a), 1L)
    # Check the counts are the elected SNP's own two groups summed (6 + 6), not a
    # sum over the gene's SNPs
    expect_equal(gene_a$coverage, 12)
})

test_that("haplotype_expression() breaks ties between electable SNPs on total coverage", {
    fixture <- make_gene_collapse_fixture()

    res <- haplotype_expression(fixture$obj)
    gene_b <- dplyr::filter(res, gene_name == "geneB")

    # Both of geneB's SNPs flip, so the higher-coverage one wins
    expect_setequal(gene_b$snp_id, fixture$snp_ids[4])
    # Verify the coverage is the elected SNP's alone (10 + 10 over its two
    # groups), confirming no summing across the gene's SNPs occurred
    expect_equal(gene_b$coverage, 20)
})

test_that("haplotype_expression() drops a gene with no flipping SNP", {
    fixture <- make_gene_collapse_fixture()

    res <- haplotype_expression(fixture$obj)

    # Confirm geneC, whose only SNP fails the flip test, is absent
    expect_false("geneC" %in% res$gene_name)
    # Ensure only the two representable genes survive, one row each
    expect_setequal(res$gene_name, c("geneA", "geneB"))
})

test_that("haplotype_expression() drops the selection-criterion columns from gene output", {
    fixture <- make_gene_collapse_fixture()

    res <- haplotype_expression(fixture$obj)

    # every surviving SNP flips by construction, so a uniformly FALSE column
    # would read as evidence rather than as the selection criterion it is
    expect_false("same_allele_dominant" %in% colnames(res))
    # Confirm the default grain is one row per (donor, gene)
    expect_equal(nrow(dplyr::distinct(res, donor, gene_name)), nrow(res))
    # Check the per-group columns are gone with the groups they described
    expect_false(any(c("active_x", "dominant_allele", "phase_contradiction") %in% colnames(res)))
})

test_that("haplotype_expression() requires gene annotation unless by_snp", {
    fixture <- make_gene_collapse_fixture()
    snp_info(fixture$obj)$gene_name <- NULL

    # Electing a representative per gene is impossible without gene annotation, and
    # the error must point at the argument that returns per-SNP output instead
    expect_error(haplotype_expression(fixture$obj), "requires gene annotation")
})

test_that("haplotype_expression() reports every phased SNP under by_snp", {
    fixture <- make_gene_collapse_fixture()

    per_snp <- haplotype_expression(fixture$obj, by_snp = TRUE, by_active_x = TRUE)

    # Verify by_snp reports every phased SNP, including the non-flipping ones no
    # gene could elect -- this is the surface for inspecting what the default omits
    expect_true(all(fixture$snp_ids %in% per_snp$snp_id))
    # Confirm the per-SNP grain carries both selection-criterion columns
    expect_true(all(c("same_allele_dominant", "phase_contradiction") %in% colnames(per_snp)))
})

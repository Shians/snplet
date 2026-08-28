# ==============================================================================
# Test Suite: Molecule-pooled haplotype expression
# Description: haplotype_expression_by_molecule() pools a gene's molecules so
#              each counts once per gene regardless of how many SNPs it spans,
#              unlike the per-SNP haplotype_expression(). pool_blocks controls
#              whether every phase block of a gene is counted or only its
#              largest.
# ==============================================================================

library(testthat)
library(Matrix)

# ------------------------------------------------------------------------------
# Test Data Setup
# ------------------------------------------------------------------------------

# One donor, 4 cells (2 X1-active, 2 X2-active). GENE1 has three heterozygous
# SNPs: snpA/snpB share phase_block 1 (the dominant block, more molecules),
# snpC sits alone in phase_block 2 (stranded). allele_on_x1 = "REF" for all
# three (i.e. REF reads = X1, ALT reads = X2), simulating post-add_molecule_phase()
# state without needing a real BAM or EM fit.
make_molecule_hap_fixture <- function() {
    ref <- matrix(0L, nrow = 3, ncol = 4)
    alt <- matrix(0L, nrow = 3, ncol = 4)
    snp_info <- data.frame(
        chrom = "chrX",
        pos = c(1000L, 1010L, 2000L),
        ref = "A",
        alt = "G",
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
    names(snp_ids) <- c("snpA", "snpB", "snpC")

    donor_snp_diag <- data.frame(
        snp_id = snp_ids,
        donor = "donor0",
        xci_informative = c(TRUE, FALSE, FALSE),
        allele_on_x1 = "REF",
        allele_on_x1_molecule = c(NA, "REF", "REF"),
        phase_block = c(1L, 1L, 2L),
        phase_source = c(NA, "read_backed", "read_backed"),
        phase_conflict = FALSE,
        zygosity = "het",
        zygosity_source = "vireo_gt",
        stringsAsFactors = FALSE
    )
    obj <- add_donor_snp_metadata(obj, donor_snp_diag, join_by = c("snp_id", "donor"), overwrite = TRUE)
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

    snp_gene_map <- data.frame(
        snp_id = snp_ids,
        gene_name = "GENE1",
        gene_strand = "+",
        ambiguous = FALSE,
        stringsAsFactors = FALSE
    )

    snp_gene_map(obj) <- snp_gene_map

    list(obj = obj, snp_ids = snp_ids, snp_gene_map = snp_gene_map)
}

# haplotype_expression_by_molecule() reads its molecule calls off the object,
# where add_molecule_phase() would have left them; these fixtures skip the BAM
# extraction and attach the calls directly.
with_molecule_calls <- function(obj, molecule_calls) {
    attr(obj, "molecule_calls") <- molecule_calls
    obj
}

# ==============================================================================

test_that("haplotype_expression_by_molecule() errors when no XCI diagnostics are stored", {
    ref <- Matrix(matrix(1L, 2, 2), sparse = TRUE)
    alt <- Matrix(matrix(1L, 2, 2), sparse = TRUE)
    snp_info <- data.frame(chrom = "chrX", pos = c(1L, 2L), ref = "A", alt = "G")
    barcode_info <- data.frame(barcode = c("c1", "c2"), donor = "donor0")
    obj <- SNPData(ref_count = ref, alt_count = alt, snp_info = snp_info, barcode_info = barcode_info)

    # Verify the function refuses to run before assign_xci() has stored a fit
    expect_error(haplotype_expression_by_molecule(obj), "Run assign_xci")
})

test_that("haplotype_expression_by_molecule() errors when molecule phase has not been added", {
    fixture <- make_molecule_hap_fixture()
    # Strip the phase_block/allele_on_x1_molecule columns to simulate a SNPData
    # that only ever had assign_xci() run, not add_molecule_phase()
    obj <- fixture$obj
    donor_snp_info <- donor_snp_info(obj)
    donor_snp_info$phase_block <- NULL
    obj@donor_snp_info <- donor_snp_info

    # Verify a clear error names the missing prerequisite step
    expect_error(
        haplotype_expression_by_molecule(obj),
        "Run add_molecule_phase"
    )
})

test_that("haplotype_expression_by_molecule() errors when the object carries no molecule calls", {
    fixture <- make_molecule_hap_fixture()

    # The fixture has stored phase but no "molecule_calls" attribute, as an
    # object subset or rebuilt after add_molecule_phase() would be.
    # Verify the error names the step that attaches them rather than returning
    # an empty result
    expect_error(
        haplotype_expression_by_molecule(fixture$obj),
        "carries no molecule calls"
    )
})

test_that("haplotype_expression_by_molecule() counts a multi-SNP molecule once, not once per SNP", {
    fixture <- make_molecule_hap_fixture()
    snpA <- fixture$snp_ids[["snpA"]]
    snpB <- fixture$snp_ids[["snpB"]]

    # A single molecule covering both snpA and snpB, both REF (i.e. X1).
    molecule_calls <- tibble::tribble(
        ~donor,
        ~barcode,
        ~umi,
        ~snp_id,
        ~allele,
        ~transcript_strand,
        "donor0",
        "cell1",
        "u1",
        snpA,
        "REF",
        "+",
        "donor0",
        "cell1",
        "u1",
        snpB,
        "REF",
        "+"
    )

    result <- haplotype_expression_by_molecule(with_molecule_calls(fixture$obj, molecule_calls), by_active_x = TRUE)
    x1_row <- dplyr::filter(result, gene_name == "GENE1", active_x == "X1")

    # Verify one molecule spanning two SNPs of the dominant block contributes
    # exactly one count, not two -- the entire point of pooling at the
    # molecule level rather than reporting per-SNP rows.
    expect_equal(x1_row$coverage, 1)
    expect_equal(x1_row$active_count, 1)
})

test_that("haplotype_expression_by_molecule() reports non-dominant-block molecules as stranded, not dropped", {
    fixture <- make_molecule_hap_fixture()
    snpA <- fixture$snp_ids[["snpA"]]
    snpB <- fixture$snp_ids[["snpB"]]
    snpC <- fixture$snp_ids[["snpC"]]

    # Two molecules in the dominant block (snpA/snpB), one molecule alone in
    # the smaller block (snpC).
    molecule_calls <- tibble::tribble(
        ~donor,
        ~barcode,
        ~umi,
        ~snp_id,
        ~allele,
        ~transcript_strand,
        "donor0",
        "cell1",
        "u1",
        snpA,
        "REF",
        "+",
        "donor0",
        "cell1",
        "u1",
        snpB,
        "REF",
        "+",
        "donor0",
        "cell2",
        "u2",
        snpA,
        "ALT",
        "+",
        "donor0",
        "cell2",
        "u2",
        snpB,
        "ALT",
        "+",
        "donor0",
        "cell3",
        "u3",
        snpC,
        "REF",
        "+"
    )

    result <- haplotype_expression_by_molecule(with_molecule_calls(fixture$obj, molecule_calls))

    # Verify the dominant block (more molecules) is the one pooled
    expect_equal(unique(result$phase_block_used), 1L)
    # Check the lone molecule in the other block is counted as stranded, not silently lost
    expect_equal(unique(result$n_stranded_molecules), 1L)
})

test_that("haplotype_expression_by_molecule() drops a molecule with a tied haplotype vote", {
    fixture <- make_molecule_hap_fixture()
    snpA <- fixture$snp_ids[["snpA"]]
    snpB <- fixture$snp_ids[["snpB"]]

    # cell1's molecule covers one REF (X1) and one ALT (X2) call -- a tie.
    molecule_calls <- tibble::tribble(
        ~donor,
        ~barcode,
        ~umi,
        ~snp_id,
        ~allele,
        ~transcript_strand,
        "donor0",
        "cell1",
        "u1",
        snpA,
        "REF",
        "+",
        "donor0",
        "cell1",
        "u1",
        snpB,
        "ALT",
        "+"
    )

    result <- haplotype_expression_by_molecule(with_molecule_calls(fixture$obj, molecule_calls))

    # Verify a tied molecule contributes to neither group's coverage
    expect_equal(sum(result$coverage), 0)
})

test_that("haplotype_expression_by_molecule() reports both active_x groups under by_active_x", {
    fixture <- make_molecule_hap_fixture()
    snpA <- fixture$snp_ids[["snpA"]]
    snpB <- fixture$snp_ids[["snpB"]]

    # Only an X1-active cell has molecules; no X2-active cell contributes any.
    molecule_calls <- tibble::tribble(
        ~donor,
        ~barcode,
        ~umi,
        ~snp_id,
        ~allele,
        ~transcript_strand,
        "donor0",
        "cell1",
        "u1",
        snpA,
        "REF",
        "+",
        "donor0",
        "cell1",
        "u1",
        snpB,
        "REF",
        "+"
    )

    result <- haplotype_expression_by_molecule(
        with_molecule_calls(fixture$obj, molecule_calls),
        by_active_x = TRUE
    )

    # Verify both groups appear even though one has no supporting molecules
    expect_equal(nrow(result), 2)
    x2_row <- dplyr::filter(result, active_x == "X2")
    # Check the zero-coverage group reads as coverage zero, not NA or missing
    expect_equal(x2_row$coverage, 0)
    expect_true(is.na(x2_row$escape_fraction))

    # Confirm the default pools that empty group into the gene's single row
    # rather than leaving a half-covered row to be tested on its own
    pooled <- haplotype_expression_by_molecule(with_molecule_calls(fixture$obj, molecule_calls))
    expect_equal(nrow(pooled), 1L)
    expect_equal(pooled$coverage, sum(result$coverage))
})

test_that("haplotype_expression_by_molecule() errors on missing molecule_calls columns", {
    fixture <- make_molecule_hap_fixture()
    bad_calls <- tibble::tibble(donor = "donor0", barcode = "cell1", snp_id = "x")

    # Verify a clear error names the missing columns
    expect_error(
        haplotype_expression_by_molecule(with_molecule_calls(fixture$obj, bad_calls)),
        "missing required column"
    )
})

test_that("haplotype_expression_by_molecule() errors when the object carries no SNP-to-gene map", {
    fixture <- make_molecule_hap_fixture()
    empty_calls <- tibble::tibble(
        donor = character(),
        barcode = character(),
        umi = character(),
        snp_id = character(),
        allele = character(),
        transcript_strand = character()
    )
    obj <- fixture$obj
    # As an object imported with an annotation lacking a strand column would be
    obj@snp_gene_map <- .empty_snp_gene_map()

    # Verify the error names how to supply the map rather than reporting no genes
    expect_error(
        haplotype_expression_by_molecule(with_molecule_calls(obj, empty_calls)),
        "no SNP-to-gene map"
    )
})

test_that("snp_gene_map<- rejects a map that does not match the object", {
    fixture <- make_molecule_hap_fixture()
    obj <- fixture$obj

    # Check a map missing the strand/ambiguous columns is refused
    expect_error(
        snp_gene_map(obj) <- tibble::tibble(snp_id = "x", gene_name = "GENE1"),
        "missing required column"
    )

    # Verify a map built against a different SNP set is refused rather than
    # silently contributing rows for SNPs the object does not hold
    expect_error(
        snp_gene_map(obj) <- tibble::tibble(
            snp_id = "chrX:999:A:G",
            gene_name = "GENE1",
            gene_strand = "+",
            ambiguous = FALSE
        ),
        "not present in snp_info"
    )
})

# ==============================================================================
# Test: strand resolution of ambiguous (multi-gene) SNPs
# ==============================================================================

# A single het SNP overlapping two opposite-strand genes, as assign_snp_genes()
# would report it: one row per candidate, both ambiguous = TRUE.
make_ambiguous_gene_fixture <- function() {
    ref <- matrix(0L, nrow = 1, ncol = 1)
    alt <- matrix(0L, nrow = 1, ncol = 1)
    snp_info <- data.frame(chrom = "chrX", pos = 5000L, ref = "A", alt = "G", stringsAsFactors = FALSE)
    barcode_info <- data.frame(barcode = "cell1", donor = "donor0", stringsAsFactors = FALSE)
    obj <- SNPData(
        ref_count = Matrix(ref, sparse = TRUE),
        alt_count = Matrix(alt, sparse = TRUE),
        snp_info = snp_info,
        barcode_info = barcode_info
    )
    snp_id <- snp_info(obj)$snp_id

    donor_snp_diag <- data.frame(
        snp_id = snp_id,
        donor = "donor0",
        xci_informative = TRUE,
        allele_on_x1 = "REF",
        allele_on_x1_molecule = NA_character_,
        phase_block = 1L,
        phase_source = NA_character_,
        phase_conflict = FALSE,
        zygosity = "het",
        zygosity_source = "vireo_gt",
        stringsAsFactors = FALSE
    )
    obj <- add_donor_snp_metadata(obj, donor_snp_diag, join_by = c("snp_id", "donor"), overwrite = TRUE)
    obj <- add_barcode_metadata(
        obj,
        data.frame(cell_id = barcode_info(obj)$cell_id, active_x = "X1", stringsAsFactors = FALSE),
        join_by = "cell_id",
        overwrite = TRUE
    )

    snp_gene_map <- data.frame(
        snp_id = c(snp_id, snp_id),
        gene_name = c("GENE_PLUS", "GENE_MINUS"),
        gene_strand = c("+", "-"),
        ambiguous = c(TRUE, TRUE),
        stringsAsFactors = FALSE
    )

    snp_gene_map(obj) <- snp_gene_map

    list(obj = obj, snp_id = snp_id, snp_gene_map = snp_gene_map)
}

test_that("haplotype_expression_by_molecule() attributes an ambiguous SNP by molecule strand", {
    fixture <- make_ambiguous_gene_fixture()

    # u1 is a "+" strand molecule (matches GENE_PLUS); u2 is "-" (matches GENE_MINUS).
    molecule_calls <- tibble::tribble(
        ~donor,
        ~barcode,
        ~umi,
        ~snp_id,
        ~allele,
        ~transcript_strand,
        "donor0",
        "cell1",
        "u1",
        fixture$snp_id,
        "REF",
        "+",
        "donor0",
        "cell1",
        "u2",
        fixture$snp_id,
        "REF",
        "-"
    )

    result <- haplotype_expression_by_molecule(with_molecule_calls(fixture$obj, molecule_calls), by_active_x = TRUE)

    # Verify each molecule is attributed only to the gene matching its own strand
    plus_row <- dplyr::filter(result, gene_name == "GENE_PLUS", active_x == "X1")
    minus_row <- dplyr::filter(result, gene_name == "GENE_MINUS", active_x == "X1")
    expect_equal(plus_row$active_count, 1)
    expect_equal(minus_row$active_count, 1)
    # Check neither gene double-counts the other strand's molecule
    expect_equal(plus_row$coverage, 1)
    expect_equal(minus_row$coverage, 1)
})

test_that("haplotype_expression_by_molecule() drops an ambiguous SNP's molecule when strand is unknown", {
    fixture <- make_ambiguous_gene_fixture()

    # A molecule with no resolvable transcript strand cannot be attributed to
    # either candidate gene, and must be excluded rather than guessed.
    molecule_calls <- tibble::tribble(
        ~donor,
        ~barcode,
        ~umi,
        ~snp_id,
        ~allele,
        ~transcript_strand,
        "donor0",
        "cell1",
        "u1",
        fixture$snp_id,
        "REF",
        NA_character_
    )

    expect_error(
        haplotype_expression_by_molecule(with_molecule_calls(fixture$obj, molecule_calls)),
        "No molecule calls"
    )
})

test_that("test_escape() prefers molecule counts once the object carries them", {
    fixture <- make_molecule_hap_fixture()
    snpA <- fixture$snp_ids[["snpA"]]
    snpB <- fixture$snp_ids[["snpB"]]

    # One molecule spanning both of the dominant block's SNPs
    molecule_calls <- tibble::tribble(
        ~donor,
        ~barcode,
        ~umi,
        ~snp_id,
        ~allele,
        ~transcript_strand,
        "donor0",
        "cell1",
        "u1",
        snpA,
        "REF",
        "+",
        "donor0",
        "cell1",
        "u1",
        snpB,
        "REF",
        "+"
    )
    obj <- with_molecule_calls(fixture$obj, molecule_calls)
    obj <- add_donor_metadata(
        obj,
        data.frame(donor = "donor0", xci_median_pi_g = 0.05, xci_rho = 0.05, stringsAsFactors = FALSE),
        join_by = "donor",
        overwrite = TRUE
    )

    result <- test_escape(obj)

    # Verify the read-backed counts are used, which count the two-SNP molecule
    # once rather than once per SNP
    expect_equal(unique(result$count_source), "molecule")
    expect_equal(result$coverage, 1)
    # Confirm the gene is tested once, not once per active-X group
    expect_equal(nrow(result), 1L)
})

test_that("haplotype_expression_by_molecule() counts other blocks' molecules when pool_blocks is TRUE", {
    fixture <- make_molecule_hap_fixture()
    snpA <- fixture$snp_ids[["snpA"]]
    snpB <- fixture$snp_ids[["snpB"]]
    snpC <- fixture$snp_ids[["snpC"]]

    # Two molecules in the dominant block (snpA/snpB), one alone in block 2.
    molecule_calls <- tibble::tribble(
        ~donor,
        ~barcode,
        ~umi,
        ~snp_id,
        ~allele,
        ~transcript_strand,
        "donor0",
        "cell1",
        "u1",
        snpA,
        "REF",
        "+",
        "donor0",
        "cell1",
        "u1",
        snpB,
        "REF",
        "+",
        "donor0",
        "cell2",
        "u2",
        snpA,
        "ALT",
        "+",
        "donor0",
        "cell2",
        "u2",
        snpB,
        "ALT",
        "+",
        "donor0",
        "cell3",
        "u3",
        snpC,
        "REF",
        "+"
    )

    result <- haplotype_expression_by_molecule(
        with_molecule_calls(fixture$obj, molecule_calls),
        pool_blocks = TRUE
    )

    # Verify the block-2 molecule is counted rather than stranded
    expect_equal(result$coverage, 3)
    # Confirm both of the gene's blocks were counted
    expect_equal(result$n_blocks_pooled, 2L)
    # Check the block-2 molecule is still reported in the stranded tally
    expect_equal(result$n_stranded_molecules, 1L)
})

test_that("haplotype_expression_by_molecule() counts only the largest block when pool_blocks is FALSE", {
    fixture <- make_molecule_hap_fixture()
    snpA <- fixture$snp_ids[["snpA"]]
    snpB <- fixture$snp_ids[["snpB"]]
    snpC <- fixture$snp_ids[["snpC"]]

    molecule_calls <- tibble::tribble(
        ~donor,
        ~barcode,
        ~umi,
        ~snp_id,
        ~allele,
        ~transcript_strand,
        "donor0",
        "cell1",
        "u1",
        snpA,
        "REF",
        "+",
        "donor0",
        "cell1",
        "u1",
        snpB,
        "REF",
        "+",
        "donor0",
        "cell2",
        "u2",
        snpA,
        "ALT",
        "+",
        "donor0",
        "cell2",
        "u2",
        snpB,
        "ALT",
        "+",
        "donor0",
        "cell3",
        "u3",
        snpC,
        "REF",
        "+"
    )

    result <- haplotype_expression_by_molecule(
        with_molecule_calls(fixture$obj, molecule_calls),
        pool_blocks = FALSE
    )

    # Verify only the dominant block's two molecules are counted
    expect_equal(result$coverage, 2)
    # Confirm a single block was counted
    expect_equal(result$n_blocks_pooled, 1L)
})

test_that("haplotype_expression_by_molecule() reports molecules from a backwards-oriented block", {
    fixture <- make_molecule_hap_fixture()
    snpA <- fixture$snp_ids[["snpA"]]
    snpC <- fixture$snp_ids[["snpC"]]

    # Block 1 (3 molecules) reads entirely active. Block 2's two molecules sit
    # in an X1-active cell yet both call the X2 allele, so that block's
    # inactive fraction is 1.0 -- impossible for genuine escape, which is
    # bounded at 0.5, and so the signature of a block oriented backwards.
    molecule_calls <- tibble::tribble(
        ~donor,
        ~barcode,
        ~umi,
        ~snp_id,
        ~allele,
        ~transcript_strand,
        "donor0",
        "cell1",
        "u1",
        snpA,
        "REF",
        "+",
        "donor0",
        "cell2",
        "u2",
        snpA,
        "REF",
        "+",
        "donor0",
        "cell3",
        "u3",
        snpA,
        "ALT",
        "+",
        "donor0",
        "cell1",
        "u4",
        snpC,
        "ALT",
        "+",
        "donor0",
        "cell1",
        "u5",
        snpC,
        "ALT",
        "+"
    )

    result <- haplotype_expression_by_molecule(
        with_molecule_calls(fixture$obj, molecule_calls),
        pool_blocks = TRUE
    )

    # Verify both of the backwards block's molecules are flagged as discordant
    expect_equal(result$discordant_block_molecules, 2L)
    # Confirm they are counted rather than dropped
    expect_equal(result$coverage, 5)
    # Check they land on the inactive haplotype, inflating escape as expected
    expect_equal(result$inactive_count, 2)
})

test_that("haplotype_expression_by_molecule() reports no discordant molecules when blocks agree", {
    fixture <- make_molecule_hap_fixture()
    snpA <- fixture$snp_ids[["snpA"]]
    snpC <- fixture$snp_ids[["snpC"]]

    molecule_calls <- tibble::tribble(
        ~donor,
        ~barcode,
        ~umi,
        ~snp_id,
        ~allele,
        ~transcript_strand,
        "donor0",
        "cell1",
        "u1",
        snpA,
        "REF",
        "+",
        "donor0",
        "cell2",
        "u2",
        snpA,
        "REF",
        "+",
        "donor0",
        "cell3",
        "u3",
        snpA,
        "ALT",
        "+",
        "donor0",
        "cell1",
        "u4",
        snpC,
        "REF",
        "+"
    )

    result <- haplotype_expression_by_molecule(
        with_molecule_calls(fixture$obj, molecule_calls),
        pool_blocks = TRUE
    )

    # Verify a block agreeing with its gene's orientation is not flagged
    expect_equal(result$discordant_block_molecules, 0L)
})

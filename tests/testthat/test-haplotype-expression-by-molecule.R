# ==============================================================================
# Test Suite: Molecule-pooled haplotype expression
# Description: haplotype_expression_by_molecule() pools molecules within a
#              gene's dominant phase block so a molecule counts once per gene
#              regardless of how many SNPs it spans, unlike the per-SNP
#              haplotype_expression().
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
    snp_ids <- get_snp_info(obj)$snp_id
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
            cell_id = get_barcode_info(obj)$cell_id,
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

    list(obj = obj, snp_ids = snp_ids, snp_gene_map = snp_gene_map)
}

# ==============================================================================

test_that("haplotype_expression_by_molecule() errors when no XCI diagnostics are stored", {
    ref <- Matrix(matrix(1L, 2, 2), sparse = TRUE)
    alt <- Matrix(matrix(1L, 2, 2), sparse = TRUE)
    snp_info <- data.frame(chrom = "chrX", pos = c(1L, 2L), ref = "A", alt = "G")
    barcode_info <- data.frame(barcode = c("c1", "c2"), donor = "donor0")
    obj <- SNPData(ref_count = ref, alt_count = alt, snp_info = snp_info, barcode_info = barcode_info)
    empty_calls <- tibble::tibble(
        donor = character(),
        barcode = character(),
        umi = character(),
        snp_id = character(),
        allele = character()
    )
    empty_map <- tibble::tibble(snp_id = character(), gene_name = character())

    # Verify the function refuses to run before assign_xci() has stored a fit
    expect_error(haplotype_expression_by_molecule(obj, empty_calls, empty_map), "Run assign_xci")
})

test_that("haplotype_expression_by_molecule() errors when molecule phase has not been added", {
    fixture <- make_molecule_hap_fixture()
    # Strip the phase_block/allele_on_x1_molecule columns to simulate a SNPData
    # that only ever had assign_xci() run, not add_molecule_phase()
    obj <- fixture$obj
    donor_snp_info <- get_donor_snp_info(obj)
    donor_snp_info$phase_block <- NULL
    obj@donor_snp_info <- donor_snp_info

    empty_calls <- tibble::tibble(
        donor = character(),
        barcode = character(),
        umi = character(),
        snp_id = character(),
        allele = character()
    )

    # Verify a clear error names the missing prerequisite step
    expect_error(
        haplotype_expression_by_molecule(obj, empty_calls, fixture$snp_gene_map),
        "Run add_molecule_phase"
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

    result <- haplotype_expression_by_molecule(fixture$obj, molecule_calls, fixture$snp_gene_map)
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

    result <- haplotype_expression_by_molecule(fixture$obj, molecule_calls, fixture$snp_gene_map)

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

    result <- haplotype_expression_by_molecule(fixture$obj, molecule_calls, fixture$snp_gene_map)

    # Verify a tied molecule contributes to neither group's coverage
    expect_equal(sum(result$coverage), 0)
})

test_that("haplotype_expression_by_molecule() reports both active_x groups even with zero coverage", {
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

    result <- haplotype_expression_by_molecule(fixture$obj, molecule_calls, fixture$snp_gene_map)

    # Verify both groups appear even though one has no supporting molecules
    expect_equal(nrow(result), 2)
    x2_row <- dplyr::filter(result, active_x == "X2")
    # Check the zero-coverage group reads as coverage zero, not NA or missing
    expect_equal(x2_row$coverage, 0)
    expect_true(is.na(x2_row$escape_fraction))
})

test_that("haplotype_expression_by_molecule() errors on missing molecule_calls columns", {
    fixture <- make_molecule_hap_fixture()
    bad_calls <- tibble::tibble(donor = "donor0", barcode = "cell1", snp_id = "x")

    # Verify a clear error names the missing columns
    expect_error(
        haplotype_expression_by_molecule(fixture$obj, bad_calls, fixture$snp_gene_map),
        "missing required column"
    )
})

test_that("haplotype_expression_by_molecule() errors on missing snp_gene_map columns", {
    fixture <- make_molecule_hap_fixture()
    empty_calls <- tibble::tibble(
        donor = character(),
        barcode = character(),
        umi = character(),
        snp_id = character(),
        allele = character(),
        transcript_strand = character()
    )
    bad_map <- tibble::tibble(snp_id = "x")

    # Verify a clear error names the missing columns
    expect_error(
        haplotype_expression_by_molecule(fixture$obj, empty_calls, bad_map),
        "missing required column"
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
    snp_id <- get_snp_info(obj)$snp_id

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
        data.frame(cell_id = get_barcode_info(obj)$cell_id, active_x = "X1", stringsAsFactors = FALSE),
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

    result <- haplotype_expression_by_molecule(fixture$obj, molecule_calls, fixture$snp_gene_map)

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
        haplotype_expression_by_molecule(fixture$obj, molecule_calls, fixture$snp_gene_map),
        "No molecule calls"
    )
})

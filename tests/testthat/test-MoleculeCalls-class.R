# ==============================================================================
# Test Suite: MoleculeCalls S4 Class
# Description: Tests for the MoleculeCalls container phase_from_molecules()
#              stores in a SNPData object's molecules slot
# ==============================================================================

library(testthat)

# ------------------------------------------------------------------------------
# Test Data Setup
# ------------------------------------------------------------------------------

test_molecule_calls <- tibble::tibble(
    library_id = "lib_A",
    barcode = c("AAA", "AAA"),
    umi = c("u1", "u1"),
    snp_id = c("snp_1", "snp_2"),
    allele = c("REF", "ALT"),
    transcript_strand = "+"
)

test_snp_gene_map <- tibble::tibble(
    snp_id = c("snp_1", "snp_2"),
    gene_name = "GENE1",
    gene_strand = "+",
    ambiguous = FALSE
)

# ==============================================================================

test_that("MoleculeCalls() with no arguments is an empty container", {
    molecules <- MoleculeCalls()

    # Verify the empty default is a valid object
    expect_s4_class(molecules, "MoleculeCalls")
    # Check each table starts with no rows but its full set of columns
    expect_equal(nrow(molecule_calls(molecules)), 0L)
    expect_named(snp_gene_map(molecules), c("snp_id", "gene_name", "gene_strand", "ambiguous"))
    expect_equal(nrow(bam_calibration(molecules)), 0L)
})

test_that("MoleculeCalls() stores its tables and returns them through the accessors", {
    molecules <- MoleculeCalls(calls = test_molecule_calls, snp_gene_map = test_snp_gene_map)

    # Confirm the accessors return what was stored
    expect_equal(molecule_calls(molecules), test_molecule_calls)
    expect_equal(snp_gene_map(molecules), test_snp_gene_map)
})

test_that("MoleculeCalls() rejects calls missing a required column", {
    bad_calls <- dplyr::select(test_molecule_calls, -library_id)

    # Ensure a call table without its key column is refused at construction
    expect_error(MoleculeCalls(calls = bad_calls), "calls is missing required column\\(s\\): library_id")
})

test_that("MoleculeCalls() rejects a gene map missing a required column", {
    bad_map <- dplyr::select(test_snp_gene_map, -gene_strand)

    # Ensure a map without strand is refused, since ambiguous SNPs need it
    expect_error(MoleculeCalls(snp_gene_map = bad_map), "snp_gene_map is missing required column\\(s\\): gene_strand")
})

test_that("show() summarises calls, molecules and SNPs", {
    molecules <- MoleculeCalls(calls = test_molecule_calls, bam_files = list(lib_A = "a.bam"))

    output <- capture.output(show(molecules))

    # Confirm the summary counts one molecule spanning two SNPs
    expect_match(output, "2 calls from 1 molecules at 2 SNPs", all = FALSE)
    # Check the BAM provenance is listed
    expect_match(output, "BAM files \\[lib_A\\]: a.bam", all = FALSE)
})

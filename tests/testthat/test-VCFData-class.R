# ==============================================================================
# Test Suite: VCFData S4 Class
# Description: Tests for VCFData class constructor, accessors, and read_vcf
# ==============================================================================

library(testthat)

# ------------------------------------------------------------------------------
# Test Data Setup
# ------------------------------------------------------------------------------

test_header <- c("##fileformat=VCFv4.2", "##source=snplet-test")
test_samples <- c("sample1", "sample2")
test_variants <- data.frame(
    CHROM = c("chr1", "chr1"),
    POS = c(100, 200),
    ID = c(".", "."),
    REF = c("A", "G"),
    ALT = c("T", "C"),
    QUAL = c(".", "."),
    FILTER = c("PASS", "PASS"),
    INFO = c(".", "."),
    stringsAsFactors = FALSE
)

create_test_vcfdata <- function() {
    VCFData(header = test_header, samples = test_samples, variants = test_variants)
}

# Write a minimal VCF file to disk and return its path
write_test_vcf_file <- function(path, with_samples = TRUE) {
    header_lines <- c(
        "##fileformat=VCFv4.2",
        "##source=snplet-test"
    )
    if (with_samples) {
        col_line <- "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\tsample2"
        data_lines <- c(
            "chr1\t100\t.\tA\tT\t.\tPASS\t.\tGT\t0/1\t0/0",
            "chr1\t200\t.\tG\tC\t.\tPASS\t.\tGT\t1/1\t0/1"
        )
    } else {
        col_line <- "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
        data_lines <- c(
            "chr1\t100\t.\tA\tT\t.\tPASS\t.",
            "chr1\t200\t.\tG\tC\t.\tPASS\t."
        )
    }
    writeLines(c(header_lines, col_line, data_lines), path)
}

# ==============================================================================

test_that("VCFData() creates a valid S4 object", {
    vcf_data <- create_test_vcfdata()

    # Verify VCFData object is created successfully
    expect_s4_class(vcf_data, "VCFData")
})

test_that("VCFData() converts non-numeric POS to numeric", {
    variants_char_pos <- test_variants
    variants_char_pos$POS <- as.character(variants_char_pos$POS)

    vcf_data <- VCFData(header = test_header, samples = test_samples, variants = variants_char_pos)

    # Verify POS was coerced to numeric
    expect_true(is.numeric(get_variants(vcf_data)$POS))
})

test_that("VCFData() errors when variants is missing required VCF columns", {
    incomplete_variants <- data.frame(CHROM = "chr1", POS = 100)

    # Verify error lists the missing required columns
    expect_error(
        VCFData(header = test_header, samples = test_samples, variants = incomplete_variants),
        "Missing required VCF columns"
    )
})

test_that("get_header() returns the header lines", {
    vcf_data <- create_test_vcfdata()

    # Verify get_header returns the original header vector
    expect_equal(get_header(vcf_data), test_header)
})

test_that("get_samples() returns the sample names", {
    vcf_data <- create_test_vcfdata()

    # Verify get_samples returns the original sample vector
    expect_equal(get_samples(vcf_data), test_samples)
})

test_that("get_variants() returns the variants data frame", {
    vcf_data <- create_test_vcfdata()

    # Verify get_variants returns a data frame with expected columns
    expect_equal(colnames(get_variants(vcf_data)), colnames(test_variants))
})

test_that("nrow(), ncol(), and dim() report variant table dimensions", {
    vcf_data <- create_test_vcfdata()

    # Verify nrow matches number of variants
    expect_equal(nrow(vcf_data), 2)
    # Verify ncol matches number of variant columns
    expect_equal(ncol(vcf_data), ncol(test_variants))
    # Verify dim combines nrow and ncol
    expect_equal(dim(vcf_data), c(2, ncol(test_variants)))
})

test_that("rownames() and colnames() delegate to the variants data frame", {
    vcf_data <- create_test_vcfdata()

    # Verify colnames matches the variants data frame columns
    expect_equal(colnames(vcf_data), colnames(test_variants))
    # Verify rownames matches the variants data frame row names
    expect_equal(rownames(vcf_data), rownames(test_variants))
})

test_that("show() prints a summary without error", {
    vcf_data <- create_test_vcfdata()

    # Verify show() runs without error
    expect_no_error(capture.output(show(vcf_data)))
})

test_that("[ subsets variants by row and preserves header and samples", {
    vcf_data <- create_test_vcfdata()

    subset_data <- vcf_data[1, ]

    # Verify subsetting returns a VCFData object
    expect_s4_class(subset_data, "VCFData")
    # Verify only the selected row is retained
    expect_equal(nrow(subset_data), 1)
    # Verify header is preserved after subsetting
    expect_equal(get_header(subset_data), test_header)
    # Verify samples are preserved after subsetting
    expect_equal(get_samples(subset_data), test_samples)
})

test_that("[ with missing row index keeps all rows", {
    vcf_data <- create_test_vcfdata()

    subset_data <- vcf_data[, seq_len(ncol(vcf_data))]

    # Verify all rows are retained when i is missing
    expect_equal(nrow(subset_data), nrow(vcf_data))
})

test_that("[ with missing column index keeps all columns", {
    vcf_data <- create_test_vcfdata()

    subset_data <- vcf_data[1, ]

    # Verify all columns are retained when j is missing
    expect_equal(ncol(subset_data), ncol(vcf_data))
})

test_that("read_vcf() parses header, samples, and variants from a plain VCF file", {
    vcf_file <- withr::local_tempfile(fileext = ".vcf")
    write_test_vcf_file(vcf_file)

    vcf_data <- read_vcf(vcf_file)

    # Verify a VCFData object is returned
    expect_s4_class(vcf_data, "VCFData")
    # Verify only the "##" header lines are kept
    expect_equal(get_header(vcf_data), c("##fileformat=VCFv4.2", "##source=snplet-test"))
    # Verify sample names are extracted from columns after FORMAT
    expect_equal(get_samples(vcf_data), c("sample1", "sample2"))
    # Verify variant rows are parsed
    expect_equal(nrow(vcf_data), 2)
})

test_that("read_vcf() strips the leading '#' from the CHROM column name", {
    vcf_file <- withr::local_tempfile(fileext = ".vcf")
    write_test_vcf_file(vcf_file)

    vcf_data <- read_vcf(vcf_file)

    # Verify the CHROM column name has no leading '#'
    expect_equal(colnames(vcf_data)[1], "CHROM")
})

test_that("read_vcf() returns no sample names when there is no FORMAT column", {
    vcf_file <- withr::local_tempfile(fileext = ".vcf")
    write_test_vcf_file(vcf_file, with_samples = FALSE)

    vcf_data <- read_vcf(vcf_file)

    # Verify no samples are detected without a FORMAT column
    expect_equal(get_samples(vcf_data), character(0))
})

test_that("read_vcf() reads gzipped VCF files", {
    vcf_file <- withr::local_tempfile(fileext = ".vcf")
    write_test_vcf_file(vcf_file)
    vcf_gz_file <- paste0(vcf_file, ".gz")
    withr::defer(unlink(vcf_gz_file))
    R.utils_available <- requireNamespace("R.utils", quietly = TRUE)
    if (R.utils_available) {
        R.utils::gzip(vcf_file, destname = vcf_gz_file, remove = FALSE)
    } else {
        con <- gzfile(vcf_gz_file, "wb")
        writeLines(readLines(vcf_file), con)
        close(con)
    }

    vcf_data <- read_vcf(vcf_gz_file)

    # Verify variants are parsed correctly from a gzipped file
    expect_equal(nrow(vcf_data), 2)
    # Verify samples are still detected in a gzipped file
    expect_equal(get_samples(vcf_data), c("sample1", "sample2"))
})

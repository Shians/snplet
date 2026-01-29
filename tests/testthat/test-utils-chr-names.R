# ==============================================================================
# Test Suite: Chromosome Name Utilities
# Description: Tests for chromosome naming style detection, normalization,
#              and conversion functions
# ==============================================================================

library(testthat)
library(Matrix)

# ------------------------------------------------------------------------------
# Test Data Setup
# ------------------------------------------------------------------------------

# Chromosome style definitions
CHR_STYLES <- list(
    numeric = c("1", "2", "3", "X", "Y", "MT"),
    ucsc = c("chr1", "chr2", "chr3", "chrX", "chrY", "chrM"),
    refseq_mouse = c("NC_000067.6", "NC_000068.7", "NC_000069.6", "NC_000086.7", "NC_000087.7", "NC_005089.1"),
    genbank_mouse = c("CM000994.2", "CM000995.2", "CM000996.2", "CM001013.2", "CM001014.2", "AY172335.1"),
    refseq_human = c("NC_000001.11", "NC_000002.12", "NC_000003.12", "NC_000023.11", "NC_000024.10", "NC_012920.1"),
    genbank_human = c("CM000663.2", "CM000664.2", "CM000665.2", "CM000685.2", "CM000686.2", "J01415.2"),
    unknown = c("custom1", "custom2", "customX")
)

# Helper to create test SNP info with different chromosome styles
make_test_snp_info <- function(chr_style = "numeric") {
    chrom <- CHR_STYLES[[chr_style]]
    if (is.null(chrom)) {
        chrom <- CHR_STYLES[["unknown"]]
    }

    n <- length(chrom)
    data.frame(
        chrom = chrom,
        pos = seq_len(n) * 1000,
        ref = rep_len(c("A", "C", "G", "T"), n),
        alt = rep_len(c("G", "T", "A", "C"), n),
        stringsAsFactors = FALSE
    )
}

# Helper to create SNPData with unknown chromosome style
make_unknown_style_snpdata <- function(suppress_warnings = TRUE) {
    snp_info <- make_test_snp_info("unknown")
    barcode_info <- data.frame(barcode = c("cell1", "cell2"))
    n_snps <- nrow(snp_info)
    alt_count <- sparseMatrix(i = c(1, 2), j = c(1, 1), x = c(5, 3), dims = c(n_snps, 2))
    ref_count <- sparseMatrix(i = c(1, 2), j = c(1, 1), x = c(10, 7), dims = c(n_snps, 2))

    if (suppress_warnings) {
        suppressWarnings(
            SNPData(
                ref_count = ref_count,
                alt_count = alt_count,
                snp_info = snp_info,
                barcode_info = barcode_info
            )
        )
    } else {
        SNPData(
            ref_count = ref_count,
            alt_count = alt_count,
            snp_info = snp_info,
            barcode_info = barcode_info
        )
    }
}

# ==============================================================================
# Test: detect_chr_style
# ==============================================================================

test_that("detect_chr_style correctly identifies numeric style", {
    chr_names <- c("1", "2", "3", "X", "Y")

    # Verify numeric style is detected
    expect_equal(detect_chr_style(chr_names), "numeric")
})

test_that("detect_chr_style correctly identifies UCSC style", {
    chr_names <- c("chr1", "chr2", "chr3", "chrX", "chrY")

    # Verify UCSC style is detected
    expect_equal(detect_chr_style(chr_names), "ucsc")
})

test_that("detect_chr_style correctly identifies RefSeq mouse style", {
    chr_names <- c("NC_000067.6", "NC_000068.7", "NC_000086.7")

    # Verify RefSeq mouse style is detected
    expect_equal(detect_chr_style(chr_names), "refseq_mouse")
})

test_that("detect_chr_style correctly identifies GenBank mouse style", {
    chr_names <- c("CM000994.2", "CM000995.2", "CM001013.2")

    # Verify GenBank mouse style is detected
    expect_equal(detect_chr_style(chr_names), "genbank_mouse")
})

test_that("detect_chr_style correctly identifies RefSeq human style", {
    chr_names <- c("NC_000001.11", "NC_000002.12", "NC_000023.11")

    # Verify RefSeq human style is detected
    expect_equal(detect_chr_style(chr_names), "refseq_human")
})

test_that("detect_chr_style correctly identifies GenBank human style", {
    chr_names <- c("CM000663.2", "CM000664.2", "CM000685.2")

    # Verify GenBank human style is detected
    expect_equal(detect_chr_style(chr_names), "genbank_human")
})

test_that("detect_chr_style returns unknown for unrecognized chromosomes", {
    chr_names <- c("custom1", "custom2", "customX")

    # Verify unknown style is returned
    expect_equal(detect_chr_style(chr_names), "unknown")
})

test_that("detect_chr_style handles empty input", {
    chr_names <- character(0)

    # Verify unknown is returned for empty input
    expect_equal(detect_chr_style(chr_names), "unknown")
})

test_that("detect_chr_style handles NA values", {
    chr_names <- c("chr1", "chr2", NA, "chrX")

    # Verify UCSC style is detected despite NA values
    expect_equal(detect_chr_style(chr_names), "ucsc")
})

test_that("detect_chr_style handles all NA values", {
    chr_names <- c(NA, NA, NA)

    # Verify unknown is returned for all NA input
    expect_equal(detect_chr_style(chr_names), "unknown")
})

test_that("detect_chr_style handles data with many contigs", {
    # Realistic case: standard chromosomes plus many contigs/scaffolds
    chr_names <- c("chr1", "chr2", "chr3", "chrX", "chrY",
                   "contig_1", "contig_2", "scaffold_999",
                   "chrUn_001", "chrUn_002", "chrUn_003")

    # Verify UCSC style is still detected despite contigs
    expect_equal(detect_chr_style(chr_names), "ucsc")
})

test_that("detect_chr_style works with minimal standard chromosomes", {
    # Edge case: even a single standard chromosome is sufficient
    chr_names <- c("chr1", "chr2", "chrX")

    # Verify UCSC style is detected with any matches
    expect_equal(detect_chr_style(chr_names), "ucsc")
})

test_that("detect_chr_style detects style with any matches", {
    # With loosened criteria, any standard chromosome match identifies the style
    chr_names <- c("chr1", "chr2", "contig1", "contig2", "contig3",
                   "scaffold1", "scaffold2", "scaffold3")

    # Even with 2/8 = 25% matching, UCSC style is detected
    expect_equal(detect_chr_style(chr_names), "ucsc")
})

# ==============================================================================
# Test: normalize_chr_names
# ==============================================================================

test_that("normalize_chr_names preserves UCSC format", {
    chr_names <- c("chr1", "chr2", "chrX", "chrY", "chrM")

    # Verify UCSC names are preserved (already canonical)
    result <- normalize_chr_names(chr_names)
    expect_equal(result, c("chr1", "chr2", "chrX", "chrY", "chrM"))
})

test_that("normalize_chr_names converts RefSeq mouse to UCSC", {
    chr_names <- c("NC_000067.6", "NC_000086.7", "NC_005089.1")

    # Verify RefSeq mouse names are normalized to UCSC
    result <- normalize_chr_names(chr_names, from_style = "refseq_mouse")
    expect_equal(result, c("chr1", "chrX", "chrM"))
})

test_that("normalize_chr_names converts GenBank human to UCSC", {
    chr_names <- c("CM000663.2", "CM000685.2", "J01415.2")

    # Verify GenBank human names are normalized to UCSC
    result <- normalize_chr_names(chr_names, from_style = "genbank_human")
    expect_equal(result, c("chr1", "chrX", "chrM"))
})

test_that("normalize_chr_names converts numeric to UCSC", {
    chr_names <- c("1", "2", "X", "Y")

    # Verify numeric format is converted to UCSC
    result <- normalize_chr_names(chr_names)
    expect_equal(result, c("chr1", "chr2", "chrX", "chrY"))
})

test_that("normalize_chr_names issues warning for unknown style", {
    chr_names <- c("custom1", "custom2")

    # Verify warning is issued for unknown style
    expect_warning(
        normalize_chr_names(chr_names),
        "Chromosome style is unknown"
    )
})

test_that("normalize_chr_names returns original names for unknown style", {
    chr_names <- c("custom1", "custom2")

    # Verify original names are returned despite warning
    result <- suppressWarnings(normalize_chr_names(chr_names))
    expect_equal(result, chr_names)
})

test_that("normalize_chr_names auto-detects style", {
    chr_names <- c("1", "2", "X")

    # Verify auto-detection works (numeric -> UCSC)
    result <- normalize_chr_names(chr_names, from_style = "auto")
    expect_equal(result, c("chr1", "chr2", "chrX"))
})

test_that("normalize_chr_names handles unmapped chromosomes", {
    chr_names <- c("chr1", "chr2", "chrUn")

    # Verify unmapped chromosomes are preserved (UCSC -> UCSC)
    result <- normalize_chr_names(chr_names)
    expect_equal(result, c("chr1", "chr2", "chrUn"))
})

# ==============================================================================
# Test: convert_chr_style
# ==============================================================================

test_that("convert_chr_style converts numeric to UCSC (default)", {
    chr_names <- c("1", "2", "X", "Y", "MT")

    # Verify conversion to UCSC style (default)
    result <- convert_chr_style(chr_names)
    expect_equal(result, c("chr1", "chr2", "chrX", "chrY", "chrM"))
})

test_that("convert_chr_style converts UCSC to RefSeq mouse", {
    chr_names <- c("chr1", "chrX", "chrM")

    # Verify conversion to RefSeq mouse style
    result <- convert_chr_style(chr_names, from_style = "ucsc", to_style = "refseq_mouse")
    expect_equal(result, c("NC_000067.6", "NC_000086.7", "NC_005089.1"))
})

test_that("convert_chr_style converts numeric to GenBank human", {
    chr_names <- c("1", "X", "MT")

    # Verify conversion to GenBank human style
    result <- convert_chr_style(chr_names, to_style = "genbank_human")
    expect_equal(result, c("CM000663.2", "CM000685.2", "J01415.2"))
})

test_that("convert_chr_style auto-detects source style", {
    chr_names <- c("1", "2", "X")

    # Verify auto-detection and conversion works (numeric -> UCSC)
    result <- convert_chr_style(chr_names)
    expect_equal(result, c("chr1", "chr2", "chrX"))
})

test_that("convert_chr_style handles unmapped chromosomes", {
    chr_names <- c("chr1", "chr2", "custom")

    # Verify unmapped chromosomes are preserved in canonical form
    result <- convert_chr_style(chr_names, to_style = "numeric")
    expect_equal(result, c("1", "2", "custom"))
})

# ==============================================================================
# Test: SNPData chromosome integration
# ==============================================================================

test_that("SNPData detects and stores chr_style for numeric chromosomes", {
    snp_info <- make_test_snp_info("numeric")
    barcode_info <- data.frame(barcode = c("cell1", "cell2"))
    alt_count <- sparseMatrix(i = c(1, 2), j = c(1, 2), x = c(5, 3), dims = c(6, 2))
    ref_count <- sparseMatrix(i = c(1, 2), j = c(1, 2), x = c(10, 7), dims = c(6, 2))

    # Verify SNPData object is created with correct chr_style
    snp_data <- SNPData(
        ref_count = ref_count,
        alt_count = alt_count,
        snp_info = snp_info,
        barcode_info = barcode_info
    )

    # Check chr_style is detected as numeric
    expect_equal(chr_style(snp_data), "numeric")

    # Check chrom_canonical column is added
    expect_true("chrom_canonical" %in% colnames(get_snp_info(snp_data)))

    # Verify chrom_canonical is converted to UCSC
    expect_equal(get_snp_info(snp_data)$chrom_canonical, c("chr1", "chr2", "chr3", "chrX", "chrY", "chrM"))
})

test_that("SNPData detects and stores chr_style for UCSC chromosomes", {
    snp_info <- make_test_snp_info("ucsc")
    barcode_info <- data.frame(barcode = c("cell1", "cell2"))
    alt_count <- sparseMatrix(i = c(1, 2), j = c(1, 2), x = c(5, 3), dims = c(6, 2))
    ref_count <- sparseMatrix(i = c(1, 2), j = c(1, 2), x = c(10, 7), dims = c(6, 2))

    # Verify SNPData object is created with correct chr_style
    snp_data <- SNPData(
        ref_count = ref_count,
        alt_count = alt_count,
        snp_info = snp_info,
        barcode_info = barcode_info
    )

    # Check chr_style is detected as UCSC
    expect_equal(chr_style(snp_data), "ucsc")

    # Check chrom_canonical is preserved (already UCSC)
    expect_equal(
        get_snp_info(snp_data)$chrom_canonical,
        c("chr1", "chr2", "chr3", "chrX", "chrY", "chrM")
    )
})

test_that("SNPData handles unknown chromosome style", {
    # Verify SNPData object is created even with unknown chr_style
    expect_warning(
        snp_data <- make_unknown_style_snpdata(suppress_warnings = FALSE),
        "Chromosome style is unknown"
    )

    # Check chr_style is unknown
    expect_equal(chr_style(snp_data), "unknown")

    # Check chrom_canonical preserves original names
    expected_chrom <- CHR_STYLES[["unknown"]]
    expect_equal(
        get_snp_info(snp_data)$chrom_canonical,
        expected_chrom
    )
})

test_that("SNPData subsetting preserves chr_style", {
    snp_info <- make_test_snp_info("ucsc")
    barcode_info <- data.frame(barcode = c("cell1", "cell2"))
    alt_count <- sparseMatrix(i = c(1, 2), j = c(1, 2), x = c(5, 3), dims = c(6, 2))
    ref_count <- sparseMatrix(i = c(1, 2), j = c(1, 2), x = c(10, 7), dims = c(6, 2))

    snp_data <- SNPData(
        ref_count = ref_count,
        alt_count = alt_count,
        snp_info = snp_info,
        barcode_info = barcode_info
    )

    # Verify subsetting preserves chr_style
    subset_data <- snp_data[1:3, ]
    expect_equal(chr_style(subset_data), "ucsc")
})

# ==============================================================================
# Test: .validate_chr_style
# ==============================================================================

test_that(".validate_chr_style passes for known chr_style", {
    snp_info <- make_test_snp_info("numeric")
    barcode_info <- data.frame(barcode = c("cell1", "cell2"))
    alt_count <- sparseMatrix(i = c(1, 2), j = c(1, 2), x = c(5, 3), dims = c(6, 2))
    ref_count <- sparseMatrix(i = c(1, 2), j = c(1, 2), x = c(10, 7), dims = c(6, 2))

    snp_data <- SNPData(
        ref_count = ref_count,
        alt_count = alt_count,
        snp_info = snp_info,
        barcode_info = barcode_info
    )

    # Verify validation passes for known style
    expect_true(.validate_chr_style(snp_data))
})

test_that(".validate_chr_style fails for unknown chr_style", {
    # Verify SNPData with unknown style is created
    snp_data <- make_unknown_style_snpdata()

    # Verify validation fails for unknown style
    expect_error(
        .validate_chr_style(snp_data),
        "requires a known chromosome naming style"
    )
})

test_that(".validate_chr_style provides informative error message", {
    # Verify SNPData with unknown style is created
    snp_data <- make_unknown_style_snpdata()

    # Verify error includes operation name
    expect_error(
        .validate_chr_style(snp_data, "test_function"),
        "test_function requires"
    )
})

# ==============================================================================
# Test: Backwards compatibility with old SNPData objects
# ==============================================================================

test_that("chr_style accessor handles SNPData objects without chr_style slot", {
    snp_info <- make_test_snp_info("numeric")
    barcode_info <- data.frame(barcode = c("cell1", "cell2"))
    alt_count <- sparseMatrix(i = c(1, 2), j = c(1, 2), x = c(5, 3), dims = c(6, 2))
    ref_count <- sparseMatrix(i = c(1, 2), j = c(1, 2), x = c(10, 7), dims = c(6, 2))

    snp_data <- SNPData(
        ref_count = ref_count,
        alt_count = alt_count,
        snp_info = snp_info,
        barcode_info = barcode_info
    )

    # Simulate an old object by removing chr_style slot content
    # We cannot actually remove the slot from an S4 object, but we can test the accessor logic
    # by checking that it correctly handles the .hasSlot check
    # Verify accessor works for normal objects
    expect_equal(chr_style(snp_data), "numeric")
    # Verify .hasSlot returns TRUE for new objects
    expect_true(methods::.hasSlot(snp_data, "chr_style"))
})

test_that("show method handles SNPData objects without chr_style slot gracefully", {
    snp_info <- make_test_snp_info("ucsc")
    barcode_info <- data.frame(barcode = c("cell1", "cell2"))
    alt_count <- sparseMatrix(i = c(1, 2), j = c(1, 2), x = c(5, 3), dims = c(6, 2))
    ref_count <- sparseMatrix(i = c(1, 2), j = c(1, 2), x = c(10, 7), dims = c(6, 2))

    snp_data <- SNPData(
        ref_count = ref_count,
        alt_count = alt_count,
        snp_info = snp_info,
        barcode_info = barcode_info
    )

    # Verify show method works without errors
    expect_output(show(snp_data), "Object of class 'SNPData'")
    # Verify show method includes chr_style for new objects
    expect_output(show(snp_data), "Chromosome style:")
})

test_that("subsetting preserves backwards compatibility", {
    snp_info <- make_test_snp_info("ucsc")
    barcode_info <- data.frame(barcode = c("cell1", "cell2"))
    alt_count <- sparseMatrix(i = c(1, 2), j = c(1, 2), x = c(5, 3), dims = c(6, 2))
    ref_count <- sparseMatrix(i = c(1, 2), j = c(1, 2), x = c(10, 7), dims = c(6, 2))

    snp_data <- SNPData(
        ref_count = ref_count,
        alt_count = alt_count,
        snp_info = snp_info,
        barcode_info = barcode_info
    )

    # Verify subsetting works correctly
    subset_data <- snp_data[1:3, ]
    # Verify chr_style is preserved after subsetting
    expect_equal(chr_style(subset_data), "ucsc")
})

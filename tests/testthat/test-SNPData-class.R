# ==============================================================================
# Test Suite: SNPData S4 Class
# Description: Tests for SNPData class constructor, accessors, and core methods
# ==============================================================================

library(testthat)
library(Matrix)

# ------------------------------------------------------------------------------
# Helper Functions
# ------------------------------------------------------------------------------

# Helper function to compare objects without names
expect_equal_unnamed <- function(object, expected, ...) {
    expect_equal(unname(object), unname(expected), ...)
}

# ------------------------------------------------------------------------------
# Test Data Setup
# ------------------------------------------------------------------------------

# Create test matrices
test_alt_count <- Matrix::Matrix(matrix(c(1, 2, 3, 4), nrow = 2, ncol = 2))
test_ref_count <- Matrix::Matrix(matrix(c(5, 6, 7, 8), nrow = 2, ncol = 2))

# Create test metadata
test_snp_info <- data.frame(
    snp_id = c("snp_1", "snp_2"),
    pos = c(100, 200),
    stringsAsFactors = FALSE
)

test_barcode_info <- data.frame(
    cell_id = c("cell_1", "cell_2"),
    donor = c("donor_1", "donor_1"),
    clonotype = c("clonotype_1", "clonotype_2"),
    stringsAsFactors = FALSE
)

# ==============================================================================

test_that("SNPData constructor works correctly", {
    # Test basic construction
    snp_data <- SNPData(
        ref_count = test_ref_count,
        alt_count = test_alt_count,
        snp_info = test_snp_info,
        barcode_info = test_barcode_info
    )
    # Verify SNPData object is created successfully
    expect_s4_class(snp_data, "SNPData")

    # Test dimension compatibility
    # Check that dim() method returns correct dimensions
    expect_equal(dim(snp_data), c(2, 2))
    # Verify nrow() method works correctly
    expect_equal(nrow(snp_data), 2)
    # Verify ncol() method works correctly
    expect_equal(ncol(snp_data), 2)

    # Test rownames and colnames
    # Check that SNP IDs are used as row names
    expect_equal(rownames(snp_data), c("snp_1", "snp_2"))
    # Check that cell IDs are used as column names
    expect_equal(colnames(snp_data), c("cell_1", "cell_2"))
})

test_that("SNPData accessors work correctly", {
    snp_data <- SNPData(
        ref_count = test_ref_count,
        alt_count = test_alt_count,
        snp_info = test_snp_info,
        barcode_info = test_barcode_info
    )

    # Test ref_count accessor
    ref_count_matrix <- ref_count(snp_data)
    # Verify ref_count accessor returns correct matrix values
    expect_equal_unnamed(ref_count_matrix, test_ref_count)
    # Check that ref_count matrix has correct row names
    expect_equal(rownames(ref_count_matrix), c("snp_1", "snp_2"))
    # Check that ref_count matrix has correct column names
    expect_equal(colnames(ref_count_matrix), c("cell_1", "cell_2"))

    # Test alt_count accessor
    alt_count_matrix <- alt_count(snp_data)
    # Verify alt_count accessor returns correct matrix values
    expect_equal_unnamed(alt_count_matrix, test_alt_count)
    # Check that alt_count matrix has correct row names
    expect_equal(rownames(alt_count_matrix), c("snp_1", "snp_2"))
    # Check that alt_count matrix has correct column names
    expect_equal(colnames(alt_count_matrix), c("cell_1", "cell_2"))

    # Test get_snp_info accessor
    snp_info <- get_snp_info(snp_data)
    # Verify SNP IDs are preserved from input
    expect_equal(snp_info$snp_id, test_snp_info$snp_id)
    # Verify SNP positions are preserved from input
    expect_equal(snp_info$pos, test_snp_info$pos)
    # Check that coverage column is automatically added
    expect_true("coverage" %in% colnames(snp_info))
    # Check that non_zero_samples column is automatically added
    expect_true("non_zero_samples" %in% colnames(snp_info))
    # Verify coverage calculation: rowSums of alt_count + ref_count
    expect_equal(snp_info$coverage, c(16, 20))  # rowSums of alt_count + ref_count
    # Verify non_zero_samples count: all samples have counts
    expect_equal(snp_info$non_zero_samples, c(2, 2))  # all samples have counts

    # Test get_barcode_info accessor
    barcode_info <- get_barcode_info(snp_data)
    # Verify cell IDs are preserved from input
    expect_equal(barcode_info$cell_id, test_barcode_info$cell_id)
    # Verify donor information is preserved from input
    expect_equal(barcode_info$donor, test_barcode_info$donor)
    # Verify clonotype information is preserved from input
    expect_equal(barcode_info$clonotype, test_barcode_info$clonotype)
    # Check that library_size column is automatically added
    expect_true("library_size" %in% colnames(barcode_info))
    # Check that non_zero_snps column is automatically added
    expect_true("non_zero_snps" %in% colnames(barcode_info))
    # Verify library_size calculation: colSums of alt_count + ref_count
    expect_equal(barcode_info$library_size, c(14, 22))  # colSums of alt_count + ref_count
    # Verify non_zero_snps count: all SNPs have counts
    expect_equal(barcode_info$non_zero_snps, c(2, 2))  # all SNPs have counts
})

test_that("SNPData subsetting works correctly", {
    snp_data <- SNPData(
        alt_count = test_alt_count,
        ref_count = test_ref_count,
        snp_info = test_snp_info,
        barcode_info = test_barcode_info
    )

    # Test subsetting by index
    subset_data <- snp_data[1, 1]
    # Verify subsetting by index returns SNPData object
    expect_s4_class(subset_data, "SNPData")
    # Check that single element subset has correct dimensions
    expect_equal(dim(subset_data), c(1, 1))

    # Test subsetting by name
    subset_data <- snp_data["snp_1", "cell_1"]
    # Verify subsetting by name returns SNPData object
    expect_s4_class(subset_data, "SNPData")
    # Check that named subset has correct dimensions
    expect_equal(dim(subset_data), c(1, 1))

    # Test that drop parameter is ignored
    subset_data <- snp_data[1, 1, drop = TRUE]
    # Verify drop=TRUE is ignored and still returns SNPData object
    expect_s4_class(subset_data, "SNPData")
    # Check that drop=TRUE doesn't affect dimensions
    expect_equal(dim(subset_data), c(1, 1))

    # Test subsetting to single row
    subset_data <- snp_data[1, ]
    # Verify single row subset returns SNPData object
    expect_s4_class(subset_data, "SNPData")
    # Check that single row subset has correct dimensions
    expect_equal(dim(subset_data), c(1, 2))

    # Test subsetting to single column
    subset_data <- snp_data[, 1]
    # Verify single column subset returns SNPData object
    expect_s4_class(subset_data, "SNPData")
    # Check that single column subset has correct dimensions
    expect_equal(dim(subset_data), c(2, 1))
})

test_that("SNPData coverage calculations work correctly", {
    # Setup
    snp_data <- SNPData(
        alt_count = test_alt_count,
        ref_count = test_ref_count,
        snp_info = test_snp_info,
        barcode_info = test_barcode_info
    )

    # Test coverage method
    expected_coverage <- Matrix::Matrix(matrix(c(6, 8, 10, 12), nrow = 2, ncol = 2))
    # Verify coverage calculation: alt_count + ref_count
    expect_equal_unnamed(coverage(snp_data), expected_coverage)

    # Test ref_fraction method
    expected_ref_fraction <- Matrix::Matrix(matrix(c(5/6, 6/8, 7/10, 8/12), nrow = 2, ncol = 2))
    # Verify ref_fraction calculation: ref_count / (ref_count + alt_count)
    expect_equal_unnamed(ref_fraction(snp_data), expected_ref_fraction)

    # Test major_allele_frac method
    expected_major_allele_frac <- Matrix::Matrix(matrix(
        c(0.5 + abs(5/6 - 0.5), 0.5 + abs(6/8 - 0.5),
          0.5 + abs(7/10 - 0.5), 0.5 + abs(8/12 - 0.5)),
        nrow = 2, ncol = 2
    ))
    # Verify major_allele_frac calculation: 0.5 + abs(ref_fraction - 0.5)
    expect_equal_unnamed(major_allele_frac(snp_data), expected_major_allele_frac)
})
test_that("SNPData handles missing IDs correctly", {
    # Create data without snp_id and cell_id
    test_snp_info_no_id <- data.frame(pos = c(100, 200))
    test_barcode_info_no_id <- data.frame(
        donor = c("donor_1", "donor_1"),
        clonotype = c("clonotype_1", "clonotype_2")
    )

    snp_data <- SNPData(
        alt_count = test_alt_count,
        ref_count = test_ref_count,
        snp_info = test_snp_info_no_id,
        barcode_info = test_barcode_info_no_id
    )

    # Check that IDs were automatically assigned
    # Verify snp_id column was auto-generated when missing
    expect_true("snp_id" %in% colnames(get_snp_info(snp_data)))
    # Verify cell_id column was auto-generated when missing
    expect_true("cell_id" %in% colnames(get_barcode_info(snp_data)))
    # Check that auto-generated SNP IDs follow expected pattern
    expect_equal(get_snp_info(snp_data)$snp_id, c("snp_1", "snp_2"))
    # Check that auto-generated cell IDs follow expected pattern
    expect_equal(get_barcode_info(snp_data)$cell_id, c("cell_1", "cell_2"))
})

test_that("SNPData validates input dimensions", {
    # Test mismatched alt_count vs ref_count dimensions
    wrong_dim_alt_count <- Matrix::Matrix(matrix(c(1, 2, 3), nrow = 3, ncol = 1))
    # Verify error when alt_count and ref_count have different row counts
    expect_error(
        SNPData(
            ref_count = test_ref_count,
            alt_count = wrong_dim_alt_count,
            snp_info = test_snp_info,
            barcode_info = test_barcode_info
        ),
        "nrow\\(alt_count\\) == nrow\\(ref_count\\)"
    )

    # Test mismatched snp_info dimensions
    wrong_dim_snp_info <- data.frame(snp_id = c("snp_1", "snp_2", "snp_3"))
    # Verify error when snp_info rows don't match matrix rows
    expect_error(
        SNPData(
            ref_count = test_ref_count,
            alt_count = test_alt_count,
            snp_info = wrong_dim_snp_info,
            barcode_info = test_barcode_info
        ),
        "nrow\\(ref_count\\) == nrow\\(snp_info\\)"
    )

    # Test mismatched barcode_info dimensions
    wrong_dim_barcode_info <- data.frame(cell_id = c("cell_1", "cell_2", "cell_3"))
    # Verify error when barcode_info rows don't match matrix columns
    expect_error(
        SNPData(
            alt_count = test_alt_count,
            ref_count = test_ref_count,
            snp_info = test_snp_info,
            barcode_info = wrong_dim_barcode_info
        ),
        "ncol\\(alt_count\\) == nrow\\(barcode_info\\)"
    )
})

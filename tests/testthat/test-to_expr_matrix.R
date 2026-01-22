# ==============================================================================
# Test Suite: Expression Matrix Conversion
# Description: Tests for to_expr_matrix method in to_expr_matrix.R
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
test_alt_count <- Matrix::Matrix(matrix(c(1, 2, 3, 4, 5, 6), nrow = 3, ncol = 2))
test_ref_count <- Matrix::Matrix(matrix(c(5, 6, 7, 8, 9, 10), nrow = 3, ncol = 2))
rownames(test_alt_count) <- c("snp_1", "snp_2", "snp_3")
colnames(test_alt_count) <- c("cell_1", "cell_2")
rownames(test_ref_count) <- c("snp_1", "snp_2", "snp_3")
colnames(test_ref_count) <- c("cell_1", "cell_2")

# Create test metadata
test_snp_info <- data.frame(
    snp_id = c("snp_1", "snp_2", "snp_3"),
    pos = c(100, 200, 300),
    stringsAsFactors = FALSE
)

test_barcode_info <- data.frame(
    cell_id = c("cell_1", "cell_2"),
    donor = c("donor_1", "donor_2"),
    clonotype = c("clonotype_1", "clonotype_2"),
    stringsAsFactors = FALSE
)

test_snp_data <- SNPData(
    ref_count = test_ref_count,
    alt_count = test_alt_count,
    snp_info = test_snp_info,
    barcode_info = test_barcode_info
)

# ==============================================================================

test_that("to_expr_matrix works with barcode level", {
    # Test barcode level aggregation
    # Verify function returns matrix
    result <- to_expr_matrix(test_snp_data, level = "barcode")
    expect_s4_class(result, "Matrix")

    # Check dimensions
    # Verify matrix has correct dimensions
    expect_equal(dim(result), c(3, 2))

    # Check row and column names
    # Verify row names are set correctly
    expect_equal(rownames(result), rownames(test_ref_count))
    # Verify column names are set correctly
    expect_equal(colnames(result), colnames(test_ref_count))

    # Test calculation manually for first element
    # Calculate expected value: (ref - alt) / (depth + 1) * log1p(depth)
    # For snp_1, cell_1: ref=5, alt=1, depth=6
    expected_val <- (5 - 1) / (6 + 1) * log1p(6)
    # Verify calculation is correct
    expect_equal(result[1, 1], expected_val, tolerance = 1e-10)
})

test_that("to_expr_matrix works with clonotype level", {
    # Test clonotype level aggregation
    # Verify function returns matrix (may be regular matrix after grouping)
    result <- to_expr_matrix(test_snp_data, level = "clonotype")
    expect_true(is.matrix(result) || inherits(result, "Matrix"))

    # Check dimensions (should aggregate by clonotypes)
    # Verify matrix has correct dimensions (3 SNPs, 2 clonotypes)
    expect_equal(dim(result), c(3, 2))

    # Check that aggregation occurred
    # Verify column names reflect clonotype grouping
    expect_equal(colnames(result), c("clonotype_1", "clonotype_2"))
})

test_that("to_expr_matrix works with donor level", {
    # Test donor level aggregation
    # Verify function returns matrix (may be regular matrix after grouping)
    result <- to_expr_matrix(test_snp_data, level = "donor")
    expect_true(is.matrix(result) || inherits(result, "Matrix"))

    # Check dimensions (should aggregate by donors)
    # Verify matrix has correct dimensions (3 SNPs, 2 donors)
    expect_equal(dim(result), c(3, 2))

    # Check that aggregation occurred
    # Verify column names reflect donor grouping
    expect_equal(colnames(result), c("donor_1", "donor_2"))
})

test_that("to_expr_matrix handles missing clonotype column", {
    # Create SNPData without clonotype column
    barcode_no_clono <- data.frame(
        cell_id = c("cell_1", "cell_2"),
        donor = c("donor_1", "donor_2"),
        stringsAsFactors = FALSE
    )

    snp_data_no_clono <- SNPData(
        ref_count = test_ref_count,
        alt_count = test_alt_count,
        snp_info = test_snp_info,
        barcode_info = barcode_no_clono
    )

    # Verify error when clonotype column is missing
    expect_error(
        to_expr_matrix(snp_data_no_clono, level = "clonotype"),
        "Clonotype information not available.*add_barcode_metadata"
    )
})

test_that("to_expr_matrix handles missing donor column", {
    # Create SNPData without donor column
    barcode_no_donor <- data.frame(
        cell_id = c("cell_1", "cell_2"),
        clonotype = c("clonotype_1", "clonotype_2"),
        stringsAsFactors = FALSE
    )

    snp_data_no_donor <- SNPData(
        ref_count = test_ref_count,
        alt_count = test_alt_count,
        snp_info = test_snp_info,
        barcode_info = barcode_no_donor
    )

    # Verify error when donor column is missing
    expect_error(
        to_expr_matrix(snp_data_no_donor, level = "donor"),
        "No donor column in barcode_info"
    )
})

test_that("to_expr_matrix handles argument matching", {
    # Test default argument
    # Verify default level is "barcode"
    result_default <- to_expr_matrix(test_snp_data)
    result_barcode <- to_expr_matrix(test_snp_data, level = "barcode")
    expect_equal(result_default, result_barcode)

    # Test partial matching
    # Verify partial matching works for level argument
    result_partial <- to_expr_matrix(test_snp_data, level = "clo")
    result_full <- to_expr_matrix(test_snp_data, level = "clonotype")
    expect_equal(result_partial, result_full)

    # Test invalid level
    # Verify error with invalid level
    expect_error(
        to_expr_matrix(test_snp_data, level = "invalid"),
        "'arg' should be one of"
    )
})

# ==============================================================================
# Tests for Missing Clonotype Data
# ==============================================================================

test_that("to_expr_matrix errors when clonotype level requested but column missing", {
    # Setup - Create SNPData without clonotype column
    test_alt_count_mini <- Matrix::Matrix(matrix(c(1, 2), nrow = 1, ncol = 2))
    test_ref_count_mini <- Matrix::Matrix(matrix(c(5, 6), nrow = 1, ncol = 2))

    test_snp_info_mini <- data.frame(
        snp_id = c("snp_1"),
        stringsAsFactors = FALSE
    )

    barcode_info_no_clonotype <- data.frame(
        cell_id = c("cell_1", "cell_2"),
        donor = c("donor_1", "donor_1"),
        stringsAsFactors = FALSE
    )

    snp_data <- SNPData(
        alt_count = test_alt_count_mini,
        ref_count = test_ref_count_mini,
        snp_info = test_snp_info_mini,
        barcode_info = barcode_info_no_clonotype
    )

    # Test error when clonotype level requested but column is missing
    # Verify to_expr_matrix errors with appropriate message for clonotype
    expect_error(
        to_expr_matrix(snp_data, level = "clonotype"),
        "Clonotype information not available.*add_barcode_metadata"
    )
})

test_that("to_expr_matrix errors when clonotype level requested but all NA", {
    # Setup - Create SNPData with all NA clonotypes
    test_alt_count_mini <- Matrix::Matrix(matrix(c(1, 2), nrow = 1, ncol = 2))
    test_ref_count_mini <- Matrix::Matrix(matrix(c(5, 6), nrow = 1, ncol = 2))

    test_snp_info_mini <- data.frame(
        snp_id = c("snp_1"),
        stringsAsFactors = FALSE
    )

    barcode_info_na_clonotype <- data.frame(
        cell_id = c("cell_1", "cell_2"),
        donor = c("donor_1", "donor_1"),
        clonotype = c(NA_character_, NA_character_),
        stringsAsFactors = FALSE
    )

    snp_data <- SNPData(
        alt_count = test_alt_count_mini,
        ref_count = test_ref_count_mini,
        snp_info = test_snp_info_mini,
        barcode_info = barcode_info_na_clonotype
    )

    # Test error when clonotype level requested but all values are NA
    # Verify to_expr_matrix errors when all clonotypes are NA
    expect_error(
        to_expr_matrix(snp_data, level = "clonotype"),
        "All clonotype values are NA.*add_barcode_metadata"
    )
})

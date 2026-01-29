# ==============================================================================
# Test Suite: Expression Matrix Conversion
# Description: Tests for to_expr_matrix method in to_expr_matrix.R
# ==============================================================================

library(testthat)
library(Matrix)
library(snplet)

# ==============================================================================
# Helper Functions
# ==============================================================================

create_test_snp_data <- function(include_clonotype = TRUE, include_donor = TRUE,
                                 n_snp = 3, n_cell = 2) {
    alt_count <- Matrix::Matrix(matrix(seq_len(n_snp * n_cell), nrow = n_snp, ncol = n_cell))
    ref_count <- Matrix::Matrix(matrix(seq_len(n_snp * n_cell) + 100, nrow = n_snp, ncol = n_cell))
    rownames(alt_count) <- paste0("snp_", seq_len(n_snp))
    colnames(alt_count) <- paste0("cell_", seq_len(n_cell))
    rownames(ref_count) <- paste0("snp_", seq_len(n_snp))
    colnames(ref_count) <- paste0("cell_", seq_len(n_cell))

    snp_info <- data.frame(
        snp_id = paste0("snp_", seq_len(n_snp)),
        pos = seq(100, by = 100, length.out = n_snp),
        stringsAsFactors = FALSE
    )

    barcode_data <- data.frame(
        cell_id = paste0("cell_", seq_len(n_cell)),
        stringsAsFactors = FALSE
    )

    if (include_donor) {
        barcode_data$donor <- c("donor_1", "donor_2")[seq_len(n_cell)]
    }

    if (include_clonotype) {
        barcode_data$clonotype <- c("clonotype_1", "clonotype_2")[seq_len(n_cell)]
    }

    SNPData(
        ref_count = ref_count,
        alt_count = alt_count,
        snp_info = snp_info,
        barcode_info = barcode_data
    )
}

# ==============================================================================
# Test: Barcode Level Aggregation
# ==============================================================================

test_that("to_expr_matrix() returns Matrix class when using barcode level", {
    snp_data <- create_test_snp_data()

    result <- to_expr_matrix(snp_data, level = "barcode")

    expect_s4_class(result, "Matrix")
})

test_that("to_expr_matrix() returns correct dimensions at barcode level", {
    snp_data <- create_test_snp_data()

    result <- to_expr_matrix(snp_data, level = "barcode")

    expect_equal(dim(result), c(3, 2))
})

test_that("to_expr_matrix() preserves SNP names at barcode level", {
    snp_data <- create_test_snp_data()

    result <- to_expr_matrix(snp_data, level = "barcode")

    expect_equal(rownames(result), c("snp_1", "snp_2", "snp_3"))
})

test_that("to_expr_matrix() preserves cell names at barcode level", {
    snp_data <- create_test_snp_data()

    result <- to_expr_matrix(snp_data, level = "barcode")

    expect_equal(colnames(result), c("cell_1", "cell_2"))
})

test_that("to_expr_matrix() is consistent at barcode level", {
    # Create test data for consistency check
    alt_count <- Matrix::Matrix(matrix(c(1, 2), nrow = 1, ncol = 2))
    ref_count <- Matrix::Matrix(matrix(c(10, 20), nrow = 1, ncol = 2))
    rownames(alt_count) <- "snp_1"
    colnames(alt_count) <- c("cell_1", "cell_2")
    rownames(ref_count) <- "snp_1"
    colnames(ref_count) <- c("cell_1", "cell_2")

    snp_info <- data.frame(
        snp_id = "snp_1",
        pos = 100,
        stringsAsFactors = FALSE
    )

    barcode_info <- data.frame(
        cell_id = c("cell_1", "cell_2"),
        stringsAsFactors = FALSE
    )

    snp_data <- SNPData(
        ref_count = ref_count,
        alt_count = alt_count,
        snp_info = snp_info,
        barcode_info = barcode_info
    )

    # Call function twice with same input
    result1 <- to_expr_matrix(snp_data, level = "barcode")
    result2 <- to_expr_matrix(snp_data, level = "barcode")

    # Verify results are identical
    expect_equal(result1, result2)
})

test_that("to_expr_matrix() distinguishes REF and ALT counts at barcode level", {
    # Create two different scenarios: one with high REF, one with high ALT
    alt_high <- Matrix::Matrix(matrix(c(10), nrow = 1, ncol = 1))
    ref_low <- Matrix::Matrix(matrix(c(1), nrow = 1, ncol = 1))
    rownames(alt_high) <- "snp_1"
    colnames(alt_high) <- "cell_1"
    rownames(ref_low) <- "snp_1"
    colnames(ref_low) <- "cell_1"

    alt_low <- Matrix::Matrix(matrix(c(1), nrow = 1, ncol = 1))
    ref_high <- Matrix::Matrix(matrix(c(10), nrow = 1, ncol = 1))
    rownames(alt_low) <- "snp_1"
    colnames(alt_low) <- "cell_1"
    rownames(ref_high) <- "snp_1"
    colnames(ref_high) <- "cell_1"

    snp_info <- data.frame(snp_id = "snp_1", pos = 100, stringsAsFactors = FALSE)
    barcode_info <- data.frame(cell_id = "cell_1", stringsAsFactors = FALSE)

    snp_data_alt_bias <- SNPData(
        ref_count = ref_low,
        alt_count = alt_high,
        snp_info = snp_info,
        barcode_info = barcode_info
    )

    snp_data_ref_bias <- SNPData(
        ref_count = ref_high,
        alt_count = alt_low,
        snp_info = snp_info,
        barcode_info = barcode_info
    )

    result_alt_bias <- as.matrix(to_expr_matrix(snp_data_alt_bias, level = "barcode"))
    result_ref_bias <- as.matrix(to_expr_matrix(snp_data_ref_bias, level = "barcode"))

    # Results should be different (formula should distinguish ALT bias vs REF bias)
    expect_false(result_alt_bias[1, 1] == result_ref_bias[1, 1])
})

test_that("to_expr_matrix() aggregates clonotypes correctly", {
    # Create data where two cells belong to same clonotype
    alt_count <- Matrix::Matrix(matrix(c(1, 2, 3, 4), nrow = 2, ncol = 2))
    ref_count <- Matrix::Matrix(matrix(c(10, 20, 30, 40), nrow = 2, ncol = 2))
    rownames(alt_count) <- c("snp_1", "snp_2")
    colnames(alt_count) <- c("cell_1", "cell_2")
    rownames(ref_count) <- c("snp_1", "snp_2")
    colnames(ref_count) <- c("cell_1", "cell_2")

    snp_info <- data.frame(
        snp_id = c("snp_1", "snp_2"),
        pos = c(100, 200),
        stringsAsFactors = FALSE
    )

    barcode_info <- data.frame(
        cell_id = c("cell_1", "cell_2"),
        clonotype = c("clono_A", "clono_A"),
        stringsAsFactors = FALSE
    )

    snp_data <- SNPData(
        ref_count = ref_count,
        alt_count = alt_count,
        snp_info = snp_info,
        barcode_info = barcode_info
    )

    result_barcode <- to_expr_matrix(snp_data, level = "barcode")
    result_clono <- to_expr_matrix(snp_data, level = "clonotype")

    # Clonotype level should have single column (one clonotype)
    expect_equal(ncol(result_clono), 1)
    # Clonotype aggregation should produce different values than cell-level
    expect_false(result_clono[1, 1] == result_barcode[1, 1])
})

test_that("to_expr_matrix() aggregates donors correctly", {
    # Create data where two cells belong to same donor
    alt_count <- Matrix::Matrix(matrix(c(1, 2, 3, 4), nrow = 2, ncol = 2))
    ref_count <- Matrix::Matrix(matrix(c(10, 20, 30, 40), nrow = 2, ncol = 2))
    rownames(alt_count) <- c("snp_1", "snp_2")
    colnames(alt_count) <- c("cell_1", "cell_2")
    rownames(ref_count) <- c("snp_1", "snp_2")
    colnames(ref_count) <- c("cell_1", "cell_2")

    snp_info <- data.frame(
        snp_id = c("snp_1", "snp_2"),
        pos = c(100, 200),
        stringsAsFactors = FALSE
    )

    barcode_info <- data.frame(
        cell_id = c("cell_1", "cell_2"),
        donor = c("donor_A", "donor_A"),
        stringsAsFactors = FALSE
    )

    snp_data <- SNPData(
        ref_count = ref_count,
        alt_count = alt_count,
        snp_info = snp_info,
        barcode_info = barcode_info
    )

    result_barcode <- to_expr_matrix(snp_data, level = "barcode")
    result_donor <- to_expr_matrix(snp_data, level = "donor")

    # Donor level should have single column (one donor)
    expect_equal(ncol(result_donor), 1)
    # Donor aggregation should produce different values than cell-level
    expect_false(result_donor[1, 1] == result_barcode[1, 1])
})

# ==============================================================================
# Test: Clonotype Level Aggregation
# ==============================================================================

test_that("to_expr_matrix() returns matrix when using clonotype level", {
    snp_data <- create_test_snp_data()

    result <- to_expr_matrix(snp_data, level = "clonotype")

    expect_true(is.matrix(result) || inherits(result, "Matrix"))
})

test_that("to_expr_matrix() returns correct dimensions at clonotype level", {
    snp_data <- create_test_snp_data()

    result <- to_expr_matrix(snp_data, level = "clonotype")

    expect_equal(dim(result), c(3, 2))
})

test_that("to_expr_matrix() groups by clonotypes correctly", {
    snp_data <- create_test_snp_data()

    result <- to_expr_matrix(snp_data, level = "clonotype")

    expect_equal(colnames(result), c("clonotype_1", "clonotype_2"))
})

# ==============================================================================
# Test: Donor Level Aggregation
# ==============================================================================

test_that("to_expr_matrix() returns matrix when using donor level", {
    snp_data <- create_test_snp_data()

    result <- to_expr_matrix(snp_data, level = "donor")

    expect_true(is.matrix(result) || inherits(result, "Matrix"))
})

test_that("to_expr_matrix() returns correct dimensions at donor level", {
    snp_data <- create_test_snp_data()

    result <- to_expr_matrix(snp_data, level = "donor")

    expect_equal(dim(result), c(3, 2))
})

test_that("to_expr_matrix() groups by donors correctly", {
    snp_data <- create_test_snp_data()

    result <- to_expr_matrix(snp_data, level = "donor")

    expect_equal(colnames(result), c("donor_1", "donor_2"))
})

# ==============================================================================
# Test: Default Argument and Partial Matching
# ==============================================================================

test_that("to_expr_matrix() defaults to barcode level", {
    snp_data <- create_test_snp_data()

    result_default <- to_expr_matrix(snp_data)
    result_barcode <- to_expr_matrix(snp_data, level = "barcode")

    expect_equal(result_default, result_barcode)
})

test_that("to_expr_matrix() supports partial matching for level argument", {
    snp_data <- create_test_snp_data()

    result_partial <- to_expr_matrix(snp_data, level = "clo")
    result_full <- to_expr_matrix(snp_data, level = "clonotype")

    expect_equal(result_partial, result_full)
})

test_that("to_expr_matrix() throws error for invalid level", {
    snp_data <- create_test_snp_data()

    expect_error(
        to_expr_matrix(snp_data, level = "invalid"),
        "'arg' should be one of"
    )
})

# ==============================================================================
# Test: Error Handling for Missing Clonotype
# ==============================================================================

test_that("to_expr_matrix() errors when clonotype column missing but requested", {
    snp_data <- create_test_snp_data(include_clonotype = FALSE)

    expect_error(
        to_expr_matrix(snp_data, level = "clonotype"),
        "Clonotype information not available.*add_barcode_metadata"
    )
})

test_that("to_expr_matrix() errors when all clonotypes are NA", {
    snp_data <- create_test_snp_data()
    snp_data@barcode_info$clonotype <- NA_character_

    expect_error(
        to_expr_matrix(snp_data, level = "clonotype"),
        "All clonotype values are NA.*add_barcode_metadata"
    )
})

# ==============================================================================
# Test: Error Handling for Missing Donor
# ==============================================================================

test_that("to_expr_matrix() errors when donor column missing but requested", {
    snp_data <- create_test_snp_data(include_donor = FALSE)

    expect_error(
        to_expr_matrix(snp_data, level = "donor"),
        "No donor column in barcode_info"
    )
})

test_that("to_expr_matrix() handles all NA donors gracefully", {
    # Create data with all NA donors - should error during aggregation
    alt_count <- Matrix::Matrix(matrix(c(1, 2), nrow = 1, ncol = 2))
    ref_count <- Matrix::Matrix(matrix(c(5, 6), nrow = 1, ncol = 2))
    rownames(alt_count) <- "snp_1"
    colnames(alt_count) <- c("cell_1", "cell_2")
    rownames(ref_count) <- "snp_1"
    colnames(ref_count) <- c("cell_1", "cell_2")

    snp_info <- data.frame(
        snp_id = "snp_1",
        pos = 100,
        stringsAsFactors = FALSE
    )

    barcode_info <- data.frame(
        cell_id = c("cell_1", "cell_2"),
        donor = c(NA_character_, NA_character_),
        stringsAsFactors = FALSE
    )

    snp_data <- SNPData(
        ref_count = ref_count,
        alt_count = alt_count,
        snp_info = snp_info,
        barcode_info = barcode_info
    )

    # Should error when trying to aggregate with all NA donors
    expect_error(
        to_expr_matrix(snp_data, level = "donor"),
        class = "error"
    )
})

# ==============================================================================
# Test: Aggregation Produces Fewer Columns Than Input
# ==============================================================================

test_that("to_expr_matrix() reduces columns when aggregating clonotypes", {
    # Create data with 3 cells but only 2 clonotypes
    alt_count <- Matrix::Matrix(matrix(c(1, 2, 3, 4, 5, 6), nrow = 2, ncol = 3))
    ref_count <- Matrix::Matrix(matrix(c(5, 6, 7, 8, 9, 10), nrow = 2, ncol = 3))
    rownames(alt_count) <- c("snp_1", "snp_2")
    colnames(alt_count) <- c("cell_1", "cell_2", "cell_3")
    rownames(ref_count) <- c("snp_1", "snp_2")
    colnames(ref_count) <- c("cell_1", "cell_2", "cell_3")

    snp_info <- data.frame(
        snp_id = c("snp_1", "snp_2"),
        pos = c(100, 200),
        stringsAsFactors = FALSE
    )

    barcode_info <- data.frame(
        cell_id = c("cell_1", "cell_2", "cell_3"),
        clonotype = c("clono_A", "clono_A", "clono_B"),
        stringsAsFactors = FALSE
    )

    snp_data <- SNPData(
        ref_count = ref_count,
        alt_count = alt_count,
        snp_info = snp_info,
        barcode_info = barcode_info
    )

    result_barcode <- to_expr_matrix(snp_data, level = "barcode")
    result_clono <- to_expr_matrix(snp_data, level = "clonotype")

    # Barcode level has 3 columns (one per cell)
    expect_equal(ncol(result_barcode), 3)
    # Clonotype level has 2 columns (one per clonotype)
    expect_equal(ncol(result_clono), 2)
})

test_that("to_expr_matrix() reduces columns when aggregating donors", {
    # Create data with 3 cells but only 2 donors
    alt_count <- Matrix::Matrix(matrix(c(1, 2, 3, 4, 5, 6), nrow = 2, ncol = 3))
    ref_count <- Matrix::Matrix(matrix(c(5, 6, 7, 8, 9, 10), nrow = 2, ncol = 3))
    rownames(alt_count) <- c("snp_1", "snp_2")
    colnames(alt_count) <- c("cell_1", "cell_2", "cell_3")
    rownames(ref_count) <- c("snp_1", "snp_2")
    colnames(ref_count) <- c("cell_1", "cell_2", "cell_3")

    snp_info <- data.frame(
        snp_id = c("snp_1", "snp_2"),
        pos = c(100, 200),
        stringsAsFactors = FALSE
    )

    barcode_info <- data.frame(
        cell_id = c("cell_1", "cell_2", "cell_3"),
        donor = c("donor_A", "donor_A", "donor_B"),
        stringsAsFactors = FALSE
    )

    snp_data <- SNPData(
        ref_count = ref_count,
        alt_count = alt_count,
        snp_info = snp_info,
        barcode_info = barcode_info
    )

    result_barcode <- to_expr_matrix(snp_data, level = "barcode")
    result_donor <- to_expr_matrix(snp_data, level = "donor")

    # Barcode level has 3 columns (one per cell)
    expect_equal(ncol(result_barcode), 3)
    # Donor level has 2 columns (one per donor)
    expect_equal(ncol(result_donor), 2)
})

# ==============================================================================
# Test: Partial Matching for Different Levels
# ==============================================================================

test_that("to_expr_matrix() supports partial matching for donor level", {
    snp_data <- create_test_snp_data()

    result_partial <- to_expr_matrix(snp_data, level = "don")
    result_full <- to_expr_matrix(snp_data, level = "donor")

    expect_equal(result_partial, result_full)
})

# ==============================================================================
# Test: Edge Cases
# ==============================================================================

test_that("to_expr_matrix() handles single SNP", {
    snp_data <- create_test_snp_data(n_snp = 1)

    result <- to_expr_matrix(snp_data, level = "barcode")

    expect_equal(dim(result), c(1, 2))
    expect_equal(rownames(result), "snp_1")
})

test_that("to_expr_matrix() handles single cell", {
    snp_data <- create_test_snp_data(n_cell = 1)

    result <- to_expr_matrix(snp_data, level = "barcode")

    expect_equal(dim(result), c(3, 1))
    expect_equal(colnames(result), "cell_1")
})

test_that("to_expr_matrix() handles single SNP and single cell", {
    snp_data <- create_test_snp_data(n_snp = 1, n_cell = 1)

    result <- to_expr_matrix(snp_data, level = "barcode")

    expect_equal(dim(result), c(1, 1))
    expect_equal(rownames(result), "snp_1")
    expect_equal(colnames(result), "cell_1")
})

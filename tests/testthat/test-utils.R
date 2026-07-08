# ==============================================================================
# Test Suite: Utility Functions
# Description: Tests for utility functions including statistics and file operations
# ==============================================================================

library(testthat)
library(Matrix)

# ==============================================================================

test_that("percentile_summary works correctly", {
    # Setup - Test data
    x <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
    result <- percentile_summary(x)

    # Test default percentile names
    expected_names <- c("min", "p10", "p25", "median", "p75", "p90", "p95", "p99", "max")
    # Verify default percentiles produce expected names
    expect_equal(names(result), expected_names)

    # Test specific quantile values
    # Verify minimum value calculation
    expect_equal(result[["min"]], 1)
    # Verify maximum value calculation
    expect_equal(result[["max"]], 10)
    # Verify median calculation
    expect_equal(result[["median"]], 5.5)
    # Verify 25th percentile calculation
    expect_equal(result[["p25"]], 3.25)
    # Verify 75th percentile calculation
    expect_equal(result[["p75"]], 7.75)

    # Test custom percentiles
    custom_result <- percentile_summary(x, percentiles = c(0.2, 0.8))
    expected_custom_names <- c("min", "p20", "median", "p80", "max")
    # Verify custom percentiles produce expected names
    expect_equal(names(custom_result), expected_custom_names)

    # Test percentiles without median (all < 0.5)
    no_median_result <- percentile_summary(x, percentiles = c(0.1, 0.25))
    expected_no_median_names <- c("min", "p10", "p25", "max")
    # Verify median is not included when no percentiles > 0.5
    expect_equal(names(no_median_result), expected_no_median_names)

    # Test single value edge case
    single_result <- percentile_summary(5)
    # Verify min equals input value for single value
    expect_equal(single_result[["min"]], 5)
    # Verify max equals input value for single value
    expect_equal(single_result[["max"]], 5)
    # Verify median equals input value for single value
    expect_equal(single_result[["median"]], 5)
})

test_that("groupedRowSums works correctly", {
    # Setup - Create test matrix
    test_matrix <- Matrix::Matrix(matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12), nrow = 3, ncol = 4))
    rownames(test_matrix) <- c("row1", "row2", "row3")

    # Test basic grouping
    groups <- c("A", "A", "B", "B")
    result <- groupedRowSums(test_matrix, groups)

    # Test dimensions and names
    # Verify result has correct number of rows (same as input matrix)
    expect_equal(nrow(result), 3)
    # Verify result has correct number of columns (number of unique groups)
    expect_equal(ncol(result), 2)
    # Verify column names match group names in sorted order
    expect_equal(colnames(result), c("A", "B"))
    # Verify row names are preserved from input matrix
    expect_equal(rownames(result), c("row1", "row2", "row3"))

    # Test calculations: Group A = cols 1,2; Group B = cols 3,4
    # Row 1: A = 1+4 = 5, B = 7+10 = 17
    # Verify sum calculation for group A, row 1
    expect_equal(result[1, "A"], 5)
    # Verify sum calculation for group B, row 1
    expect_equal(result[1, "B"], 17)

    # Test uneven group sizes
    single_groups <- c("A", "B", "B", "B")
    single_result <- groupedRowSums(test_matrix, single_groups)
    # Verify result has correct number of groups
    expect_equal(ncol(single_result), 2)
    # Verify sum for single-member group A
    expect_equal(single_result[1, "A"], 1) # Only first column for group A

    # Test single group
    same_groups <- rep("A", 4)
    same_result <- groupedRowSums(test_matrix, same_groups)
    # Verify single group produces single column
    expect_equal(ncol(same_result), 1)
    # Verify column name matches group name
    expect_equal(colnames(same_result), "A")
})

test_that("groupedRowSums validates group length and missing values", {
    test_matrix <- Matrix::Matrix(matrix(1:12, nrow = 3, ncol = 4))

    # Verify error when group length does not match matrix columns
    expect_error(
        groupedRowSums(test_matrix, c("A", "B")),
        "Length of groups must match the number of columns in x"
    )

    # Verify error when groups contain NA values
    expect_error(
        groupedRowSums(test_matrix, c("A", NA, "B", "B")),
        "groups must not contain NA values"
    )
})

test_that("check_file works correctly", {
    # Setup - Create temporary file
    temp_file <- tempfile()
    writeLines("test content", temp_file)

    # Test existing file (should not error)
    # Verify no error when checking existing file
    expect_no_error(check_file(temp_file))

    # Cleanup
    unlink(temp_file)

    # Test non-existing file (should error)
    non_existing_file <- tempfile()
    # Verify error when checking non-existing file
    expect_error(
        check_file(non_existing_file),
        "Required file not found:"
    )
})

test_that("percentile_summary handles extreme percentiles", {
    # Setup - Test data
    x <- c(1, 2, 3, 4, 5)

    # Test with all percentiles below 0.5
    result <- percentile_summary(x, percentiles = c(0.1, 0.2, 0.3, 0.4))

    # Verify median is not included for all-low percentiles
    expect_false("median" %in% names(result))

    # Verify min and max are always included
    expect_true("min" %in% names(result))
    # Verify max is always included regardless of percentile selection
    expect_true("max" %in% names(result))
})

test_that("percentile_summary handles empty input", {
    # Test with empty vector
    expect_warning(
        result <- percentile_summary(numeric(0))
    )

    # Verify result contains expected structure
    expect_true(is.numeric(result))
    # Verify all core percentile names are present despite empty input
    expect_true(all(c("min", "median", "max") %in% names(result)))
})

test_that("groupedRowSums preserves matrix sparsity", {
    # Setup - Sparse matrix with values
    sparse_matrix <- Matrix::Matrix(
        matrix(c(1, 0, 0, 0, 2, 0, 3, 0, 0), nrow = 3, ncol = 3),
        sparse = TRUE
    )

    # Test with sparse matrix
    groups <- c("A", "B", "A")
    result <- groupedRowSums(sparse_matrix, groups)

    # Verify result structure
    expect_equal(nrow(result), 3)
    # Verify result has correct number of groups
    expect_equal(ncol(result), 2)

    # Verify calculations are correct with sparse matrices
    # Row 1: A columns (1,3) = 1+3, B column (2) = 0
    expect_equal(as.numeric(result[1, "A"]), 4)
    # Row 2: A columns (1,3) = 0+0, B column (2) = 2
    expect_equal(as.numeric(result[2, "B"]), 2)
})

# ==============================================================================

test_that("binom_test matches binom.test for alternative = 'greater'", {
    # Setup: a grid of (x, n) values and a fixed null probability
    n_vals <- c(1, 5, 10, 25, 100)
    p_null <- 0.3

    cases <- do.call(rbind, lapply(n_vals, function(n) {
        data.frame(x = 0:n, n = n)
    }))

    # Reference p-values from stats::binom.test
    expected <- mapply(
        function(x, n) stats::binom.test(x, n, p = p_null, alternative = "greater")$p.value,
        cases$x, cases$n
    )

    # Vectorised wrapper
    actual <- binom_test(cases$x, cases$n, p = p_null, alternative = "greater")

    # Verify all p-values match binom.test exactly across the grid
    expect_equal(actual, expected)
})

test_that("binom_test matches binom.test for alternative = 'less'", {
    # Setup: a grid of (x, n) values and a fixed null probability
    n_vals <- c(1, 5, 10, 25, 100)
    p_null <- 0.7

    cases <- do.call(rbind, lapply(n_vals, function(n) {
        data.frame(x = 0:n, n = n)
    }))

    # Reference p-values from stats::binom.test
    expected <- mapply(
        function(x, n) stats::binom.test(x, n, p = p_null, alternative = "less")$p.value,
        cases$x, cases$n
    )

    # Vectorised wrapper
    actual <- binom_test(cases$x, cases$n, p = p_null, alternative = "less")

    # Verify all p-values match binom.test exactly across the grid
    expect_equal(actual, expected)
})

test_that("binom_test handles boundary cases x = 0 and x = n", {
    # Verify p-value is 1 when x = 0 under 'greater' (P(X >= 0) = 1)
    expect_equal(binom_test(0, 10, 0.5, alternative = "greater"), 1)
    # Verify p-value is 1 when x = n under 'less' (P(X <= n) = 1)
    expect_equal(binom_test(10, 10, 0.5, alternative = "less"), 1)

    # Confirm boundary results agree with binom.test
    expect_equal(
        binom_test(0, 10, 0.5, alternative = "greater"),
        stats::binom.test(0, 10, 0.5, alternative = "greater")$p.value
    )
    expect_equal(
        binom_test(10, 10, 0.5, alternative = "less"),
        stats::binom.test(10, 10, 0.5, alternative = "less")$p.value
    )
})

test_that("binom_test recycles p across vectorised inputs", {
    # Setup: vary p alongside x and n
    x <- c(2, 5, 8)
    n <- c(10, 10, 10)
    p <- c(0.1, 0.5, 0.9)

    expected <- mapply(
        function(x, n, p) stats::binom.test(x, n, p, alternative = "greater")$p.value,
        x, n, p
    )

    actual <- binom_test(x, n, p, alternative = "greater")

    # Verify recycled p produces matching p-values element-wise
    expect_equal(actual, expected)
})

test_that("binom_test rejects invalid inputs", {
    # Ensure x outside [0, n] is rejected
    expect_error(binom_test(-1, 10, 0.5), "0 <= x <= n")
    # Ensure x greater than n is rejected
    expect_error(binom_test(11, 10, 0.5), "0 <= x <= n")
    # Ensure p outside [0, 1] is rejected
    expect_error(binom_test(2, 10, 1.5), "p must be in")
    # Confirm two-sided is not an accepted alternative
    expect_error(binom_test(2, 10, 0.5, alternative = "two.sided"))
})

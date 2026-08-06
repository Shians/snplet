# ==============================================================================
# Test Suite: Utility Functions
# Description: Tests for utility functions including statistics and file operations
# ==============================================================================

library(testthat)
library(Matrix)

# ==============================================================================

test_that("percentile_summary uses default percentile names", {
    x <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)

    result <- percentile_summary(x)

    expected_names <- c("min", "p10", "p25", "median", "p75", "p90", "p95", "p99", "max")
    # Verify default percentiles produce expected names
    expect_equal(names(result), expected_names)
})

test_that("percentile_summary calculates correct percentile values", {
    x <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)

    result <- percentile_summary(x)

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
})

test_that("percentile_summary honors custom percentiles", {
    x <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)

    custom_result <- percentile_summary(x, percentiles = c(0.2, 0.8))

    expected_custom_names <- c("min", "p20", "median", "p80", "max")
    # Verify custom percentiles produce expected names
    expect_equal(names(custom_result), expected_custom_names)
})

test_that("percentile_summary omits median when no percentile exceeds 0.5", {
    x <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)

    no_median_result <- percentile_summary(x, percentiles = c(0.1, 0.25))

    expected_no_median_names <- c("min", "p10", "p25", "max")
    # Verify median is not included when no percentiles > 0.5
    expect_equal(names(no_median_result), expected_no_median_names)
})

test_that("percentile_summary handles a single-value input", {
    single_result <- percentile_summary(5)

    # Verify min equals input value for single value
    expect_equal(single_result[["min"]], 5)
    # Verify max equals input value for single value
    expect_equal(single_result[["max"]], 5)
    # Verify median equals input value for single value
    expect_equal(single_result[["median"]], 5)
})

test_that("groupedRowSums returns correct dimensions and names for basic grouping", {
    test_matrix <- Matrix::Matrix(matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12), nrow = 3, ncol = 4))
    rownames(test_matrix) <- c("row1", "row2", "row3")
    groups <- c("A", "A", "B", "B")

    result <- groupedRowSums(test_matrix, groups)

    # Verify result has correct number of rows (same as input matrix)
    expect_equal(nrow(result), 3)
    # Verify result has correct number of columns (number of unique groups)
    expect_equal(ncol(result), 2)
    # Verify column names match group names in sorted order
    expect_equal(colnames(result), c("A", "B"))
    # Verify row names are preserved from input matrix
    expect_equal(rownames(result), c("row1", "row2", "row3"))
})

test_that("groupedRowSums sums the correct columns per group", {
    test_matrix <- Matrix::Matrix(matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12), nrow = 3, ncol = 4))
    rownames(test_matrix) <- c("row1", "row2", "row3")
    groups <- c("A", "A", "B", "B")

    result <- groupedRowSums(test_matrix, groups)

    # Group A = cols 1,2; Group B = cols 3,4. Row 1: A = 1+4 = 5, B = 7+10 = 17
    # Verify sum calculation for group A, row 1
    expect_equal(result[1, "A"], 5)
    # Verify sum calculation for group B, row 1
    expect_equal(result[1, "B"], 17)
})

test_that("groupedRowSums handles uneven group sizes", {
    test_matrix <- Matrix::Matrix(matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12), nrow = 3, ncol = 4))
    rownames(test_matrix) <- c("row1", "row2", "row3")
    single_groups <- c("A", "B", "B", "B")

    single_result <- groupedRowSums(test_matrix, single_groups)

    # Verify result has correct number of groups
    expect_equal(ncol(single_result), 2)
    # Verify sum for single-member group A (only first column for group A)
    expect_equal(single_result[1, "A"], 1)
})

test_that("groupedRowSums collapses to a single column when all groups match", {
    test_matrix <- Matrix::Matrix(matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12), nrow = 3, ncol = 4))
    rownames(test_matrix) <- c("row1", "row2", "row3")
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

test_that("check_file does not error for an existing file", {
    temp_file <- withr::local_tempfile()
    writeLines("test content", temp_file)

    # Verify no error when checking existing file
    expect_no_error(check_file(temp_file))
})

test_that("check_file errors for a non-existing file", {
    non_existing_file <- withr::local_tempfile()

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

    cases <- do.call(
        rbind,
        lapply(n_vals, function(n) {
            data.frame(x = 0:n, n = n)
        })
    )

    # Reference p-values from stats::binom.test
    expected <- mapply(
        function(x, n) stats::binom.test(x, n, p = p_null, alternative = "greater")$p.value,
        cases$x,
        cases$n
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

    cases <- do.call(
        rbind,
        lapply(n_vals, function(n) {
            data.frame(x = 0:n, n = n)
        })
    )

    # Reference p-values from stats::binom.test
    expected <- mapply(
        function(x, n) stats::binom.test(x, n, p = p_null, alternative = "less")$p.value,
        cases$x,
        cases$n
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
        x,
        n,
        p
    )

    actual <- binom_test(x, n, p, alternative = "greater")

    # Verify recycled p produces matching p-values element-wise
    expect_equal(actual, expected)
})

test_that("binom_test rejects x below 0", {
    # Ensure x outside [0, n] is rejected
    expect_error(binom_test(-1, 10, 0.5), "0 <= x <= n")
})

test_that("binom_test rejects x greater than n", {
    # Ensure x greater than n is rejected
    expect_error(binom_test(11, 10, 0.5), "0 <= x <= n")
})

test_that("binom_test rejects p outside [0, 1]", {
    # Ensure p outside [0, 1] is rejected
    expect_error(binom_test(2, 10, 1.5), "p must be in")
})

test_that("binom_test rejects unsupported alternative", {
    # Confirm two-sided is not an accepted alternative
    expect_error(binom_test(2, 10, 0.5, alternative = "two.sided"))
})

# ==============================================================================

test_that("betabinom_test reduces to binom_test when rho = 0", {
    # Setup: a grid of (x, n) values with rho = 0, where the beta-binomial
    # collapses to the plain binomial
    n_vals <- c(1, 5, 10, 25)
    p_null <- 0.3

    cases <- do.call(rbind, lapply(n_vals, function(n) data.frame(x = 0:n, n = n)))

    expected <- binom_test(cases$x, cases$n, p = p_null, alternative = "greater")
    actual <- betabinom_test(cases$x, cases$n, p = p_null, rho = 0, alternative = "greater")

    # Verify beta-binomial p-values match the binomial ones at rho = 0
    expect_equal(actual, expected, tolerance = 1e-6)
})

test_that("betabinom_test matches VGAM::pbetabinom directly for alternative = 'greater'", {
    # Setup: a grid of (x, n) values with non-trivial overdispersion
    n_vals <- c(5, 10, 25)
    p_null <- 0.1
    rho <- 0.3

    cases <- do.call(rbind, lapply(n_vals, function(n) data.frame(x = 0:n, n = n)))

    expected <- 1 - VGAM::pbetabinom(cases$x - 1, cases$n, p_null, rho)
    expected[cases$x == 0] <- 1
    actual <- betabinom_test(cases$x, cases$n, p = p_null, rho = rho, alternative = "greater")

    # Verify p-values match the direct beta-binomial CDF computation
    expect_equal(actual, expected)
})

test_that("betabinom_test matches VGAM::pbetabinom directly for alternative = 'less'", {
    # Setup: a grid of (x, n) values with non-trivial overdispersion
    n_vals <- c(5, 10, 25)
    p_null <- 0.7
    rho <- 0.3

    cases <- do.call(rbind, lapply(n_vals, function(n) data.frame(x = 0:n, n = n)))

    expected <- VGAM::pbetabinom(cases$x, cases$n, p_null, rho)
    expected[cases$x == cases$n] <- 1
    actual <- betabinom_test(cases$x, cases$n, p = p_null, rho = rho, alternative = "less")

    # Verify p-values match the direct beta-binomial CDF computation
    expect_equal(actual, expected)
})

test_that("betabinom_test handles boundary cases x = 0 and x = n", {
    # Verify p-value is 1 when x = 0 under 'greater' (P(X >= 0) = 1)
    expect_equal(betabinom_test(0, 10, 0.5, rho = 0.2, alternative = "greater"), 1)
    # Verify p-value is 1 when x = n under 'less' (P(X <= n) = 1)
    expect_equal(betabinom_test(10, 10, 0.5, rho = 0.2, alternative = "less"), 1)
})

test_that("betabinom_test produces larger p-values than binom_test for the same skewed observation", {
    # Higher overdispersion should widen the null distribution's tails, so an
    # observation that looks significant under a binomial null looks less so
    # under a beta-binomial null with the same mean.
    binom_p <- binom_test(8, 20, 0.1, alternative = "greater")
    betabinom_p <- betabinom_test(8, 20, 0.1, rho = 0.3, alternative = "greater")

    # Verify overdispersion inflates the p-value relative to the binomial null
    expect_true(betabinom_p > binom_p)
})

test_that("betabinom_test rejects x below 0", {
    # Ensure x outside [0, n] is rejected
    expect_error(betabinom_test(-1, 10, 0.5, rho = 0.2), "0 <= x <= n")
})

test_that("betabinom_test rejects x greater than n", {
    # Ensure x greater than n is rejected
    expect_error(betabinom_test(11, 10, 0.5, rho = 0.2), "0 <= x <= n")
})

test_that("betabinom_test rejects p outside [0, 1]", {
    # Ensure p outside [0, 1] is rejected
    expect_error(betabinom_test(2, 10, 1.5, rho = 0.2), "p must be in")
})

test_that("betabinom_test rejects rho outside [0, 1)", {
    # Ensure negative rho is rejected
    expect_error(betabinom_test(2, 10, 0.5, rho = -0.1), "rho must be in")
    # Ensure rho = 1 is rejected (degenerate beta-binomial)
    expect_error(betabinom_test(2, 10, 0.5, rho = 1), "rho must be in")
})

test_that("betabinom_test rejects unsupported alternative", {
    # Confirm two-sided is not an accepted alternative
    expect_error(betabinom_test(2, 10, 0.5, rho = 0.2, alternative = "two.sided"))
})

test_that("betabinom_test clamps p-values that would go slightly negative from VGAM floating-point overshoot", {
    # VGAM::pbetabinom(x - 1, n, p, rho = 0) is observed to overshoot 1 by a
    # few multiples of machine epsilon for some inputs (e.g. x = 337, n = 875,
    # p = 0.0336), which would otherwise make 1 - pbetabinom(...) negative.
    cases <- data.frame(
        x = c(337, 2524, 80),
        n = c(875, 10315, 200),
        p = c(0.0336, 0.0315, 0.0236)
    )

    actual <- betabinom_test(cases$x, cases$n, cases$p, rho = 0, alternative = "greater")

    # Verify no p-value is negative despite the underlying CDF overshooting 1
    expect_true(all(actual >= 0))
    # Verify no p-value exceeds 1
    expect_true(all(actual <= 1))
})

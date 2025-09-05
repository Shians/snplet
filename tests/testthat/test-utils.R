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

test_that("groupedRowMeans works correctly", {
    # Setup - Create test matrix
    test_matrix <- Matrix::Matrix(matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12), nrow = 3, ncol = 4))
    rownames(test_matrix) <- c("row1", "row2", "row3")

    # Test basic grouping
    groups <- c("A", "A", "B", "B")
    result <- groupedRowMeans(test_matrix, groups)

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
    # Row 1: A = (1+4)/2 = 2.5, B = (7+10)/2 = 8.5
    # Verify mean calculation for group A, row 1
    expect_equal(result[1, "A"], 2.5)
    # Verify mean calculation for group B, row 1
    expect_equal(result[1, "B"], 8.5)

    # Test uneven group sizes
    single_groups <- c("A", "B", "B", "B")
    single_result <- groupedRowMeans(test_matrix, single_groups)
    # Verify result has correct number of groups
    expect_equal(ncol(single_result), 2)
    # Verify mean for single-member group A
    expect_equal(single_result[1, "A"], 1) # Only first column for group A

    # Test single group
    same_groups <- rep("A", 4)
    same_result <- groupedRowMeans(test_matrix, same_groups)
    # Verify single group produces single column
    expect_equal(ncol(same_result), 1)
    # Verify column name matches group name
    expect_equal(colnames(same_result), "A")
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

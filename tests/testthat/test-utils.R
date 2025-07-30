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
    expect_equal(names(result), expected_names)
    
    # Test specific quantile values
    expect_equal(result[["min"]], 1)
    expect_equal(result[["max"]], 10)
    expect_equal(result[["median"]], 5.5)
    expect_equal(result[["p25"]], 3.25)
    expect_equal(result[["p75"]], 7.75)
    
    # Test custom percentiles
    custom_result <- percentile_summary(x, percentiles = c(0.2, 0.8))
    expected_custom_names <- c("min", "p20", "median", "p80", "max")
    expect_equal(names(custom_result), expected_custom_names)
    
    # Test percentiles without median (all < 0.5)
    no_median_result <- percentile_summary(x, percentiles = c(0.1, 0.25))
    expected_no_median_names <- c("min", "p10", "p25", "max")
    expect_equal(names(no_median_result), expected_no_median_names)
    
    # Test single value edge case
    single_result <- percentile_summary(5)
    expect_equal(single_result[["min"]], 5)
    expect_equal(single_result[["max"]], 5)
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
    expect_equal(nrow(result), 3)
    expect_equal(ncol(result), 2)
    expect_equal(colnames(result), c("A", "B"))
    expect_equal(rownames(result), c("row1", "row2", "row3"))
    
    # Test calculations: Group A = cols 1,2; Group B = cols 3,4
    # Row 1: A = (1+4)/2 = 2.5, B = (7+10)/2 = 8.5
    expect_equal(result[1, "A"], 2.5)
    expect_equal(result[1, "B"], 8.5)
    
    # Test uneven group sizes
    single_groups <- c("A", "B", "B", "B")
    single_result <- groupedRowMeans(test_matrix, single_groups)
    expect_equal(ncol(single_result), 2)
    expect_equal(single_result[1, "A"], 1)  # Only first column for group A
    
    # Test single group
    same_groups <- rep("A", 4)
    same_result <- groupedRowMeans(test_matrix, same_groups)
    expect_equal(ncol(same_result), 1)
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
    expect_equal(nrow(result), 3)
    expect_equal(ncol(result), 2)
    expect_equal(colnames(result), c("A", "B"))
    expect_equal(rownames(result), c("row1", "row2", "row3"))
    
    # Test calculations: Group A = cols 1,2; Group B = cols 3,4
    # Row 1: A = 1+4 = 5, B = 7+10 = 17
    expect_equal(result[1, "A"], 5)
    expect_equal(result[1, "B"], 17)
    
    # Test uneven group sizes
    single_groups <- c("A", "B", "B", "B")
    single_result <- groupedRowSums(test_matrix, single_groups)
    expect_equal(ncol(single_result), 2)
    expect_equal(single_result[1, "A"], 1)  # Only first column for group A
    
    # Test single group
    same_groups <- rep("A", 4)
    same_result <- groupedRowSums(test_matrix, same_groups)
    expect_equal(ncol(same_result), 1)
    expect_equal(colnames(same_result), "A")
})

test_that("check_file works correctly", {
    # Setup - Create temporary file
    temp_file <- tempfile()
    writeLines("test content", temp_file)
    
    # Test existing file (should not error)
    expect_no_error(check_file(temp_file))
    
    # Cleanup
    unlink(temp_file)
    
    # Test non-existing file (should error)
    non_existing_file <- tempfile()
    expect_error(
        check_file(non_existing_file),
        "Required file not found:"
    )
})
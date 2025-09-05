# ==============================================================================
# Test Suite: Minor Allele Frequency Testing
# Description: Tests for statistical significance testing of minor allele frequencies
# ==============================================================================

library(testthat)
library(dplyr)

# ==============================================================================

test_that("test_maf works correctly with basic input", {
    # Setup - Create test data frame with known values
    test_df <- data.frame(
        ref_count = c(10, 5, 20),
        alt_count = c(2, 8, 1),
        total_count = c(12, 13, 21)
    )

    # Execute test_maf function
    result <- test_maf(test_df)

    # Verify result is a data frame
    expect_s3_class(result, "data.frame")

    # Check that all original columns are preserved
    expect_true("ref_count" %in% colnames(result))
    expect_true("alt_count" %in% colnames(result))
    expect_true("total_count" %in% colnames(result))

    # Check that new columns are added
    expect_true("minor_allele_count" %in% colnames(result))
    expect_true("p_val" %in% colnames(result))
    expect_true("adj_p_val" %in% colnames(result))

    # Verify minor_allele_count calculation (minimum of ref_count and alt_count)
    expected_minor_allele <- pmin(test_df$ref_count, test_df$alt_count)
    expect_equal(result$minor_allele_count, expected_minor_allele)

    # Check that result has same number of rows as input
    expect_equal(nrow(result), nrow(test_df))

    # Verify p-values are numeric and within valid range [0,1]
    expect_type(result$p_val, "double")
    expect_true(all(result$p_val >= 0 & result$p_val <= 1))

    # Verify adjusted p-values are numeric and within valid range [0,1]
    expect_type(result$adj_p_val, "double")
    expect_true(all(result$adj_p_val >= 0 & result$adj_p_val <= 1))

    # Check that adjusted p-values are >= original p-values (BH correction property)
    expect_true(all(result$adj_p_val >= result$p_val))
})

test_that("test_maf handles custom null hypothesis probability", {
    # Setup - Create simple test data
    test_df <- data.frame(
        ref_count = c(8, 15),
        alt_count = c(2, 5),
        total_count = c(10, 20)
    )

    # Test with different null hypothesis probabilities
    result_p05 <- test_maf(test_df, p = 0.05)
    result_p20 <- test_maf(test_df, p = 0.20)

    # Verify both results have required columns
    expect_true(all(c("minor_allele_count", "p_val", "adj_p_val") %in% colnames(result_p05)))
    expect_true(all(c("minor_allele_count", "p_val", "adj_p_val") %in% colnames(result_p20)))

    # Check that different null hypotheses produce different p-values
    expect_false(identical(result_p05$p_val, result_p20$p_val))

    # Verify p-values are in valid range for both tests
    expect_true(all(result_p05$p_val >= 0 & result_p05$p_val <= 1))
    expect_true(all(result_p20$p_val >= 0 & result_p20$p_val <= 1))
})

test_that("test_maf handles edge cases correctly", {
    # Setup - Test with edge case values
    edge_cases_df <- data.frame(
        ref_count = c(0, 10, 5),
        alt_count = c(10, 0, 5),
        total_count = c(10, 10, 10)
    )

    # Execute test
    result <- test_maf(edge_cases_df)

    # Verify minor allele count calculation for edge cases
    # Case 1: ref=0, alt=10 -> minor_allele = 0
    expect_equal(result$minor_allele_count[1], 0)
    # Case 2: ref=10, alt=0 -> minor_allele = 0
    expect_equal(result$minor_allele_count[2], 0)
    # Case 3: ref=5, alt=5 -> minor_allele = 5
    expect_equal(result$minor_allele_count[3], 5)

    # Check that all results are valid
    expect_true(all(result$p_val >= 0 & result$p_val <= 1))
    expect_true(all(result$adj_p_val >= 0 & result$adj_p_val <= 1))
})

test_that("test_maf handles fractional total_count values", {
    # Setup - Data with fractional total counts (should be ceiling'd)
    fractional_df <- data.frame(
        ref_count = c(8, 12),
        alt_count = c(2, 3),
        total_count = c(10.7, 15.3)
    )

    # Execute test
    result <- test_maf(fractional_df)

    # Verify function completes without error
    expect_s3_class(result, "data.frame")

    # Check that results are valid (function should ceiling the total_count internally)
    expect_true(all(result$p_val >= 0 & result$p_val <= 1))
    expect_true(all(result$adj_p_val >= 0 & result$adj_p_val <= 1))
})

test_that("test_maf validates input data frame", {
    # Test with non-data.frame input
    # Verify error when input is not a data frame
    expect_error(
        test_maf("not_a_dataframe"),
        "is\\(x, \"data\\.frame\"\\) is not TRUE"
    )

    # Verify error with NULL input
    expect_error(
        test_maf(NULL),
        "is\\(x, \"data\\.frame\"\\) is not TRUE"
    )

    # Verify error with list input
    expect_error(
        test_maf(list(a = 1, b = 2)),
        "is\\(x, \"data\\.frame\"\\) is not TRUE"
    )
})

test_that("test_maf validates required columns", {
    # Test with missing ref_count column
    missing_ref <- data.frame(
        alt_count = c(2, 5),
        total_count = c(10, 15)
    )

    # Verify error when ref_count column is missing
    expect_error(
        test_maf(missing_ref),
        "Missing required columns: ref_count"
    )

    # Test with missing alt_count column
    missing_alt <- data.frame(
        ref_count = c(8, 10),
        total_count = c(10, 15)
    )

    # Verify error when alt_count column is missing
    expect_error(
        test_maf(missing_alt),
        "Missing required columns: alt_count"
    )

    # Test with missing total_count column
    missing_total <- data.frame(
        ref_count = c(8, 10),
        alt_count = c(2, 5)
    )

    # Verify error when total_count column is missing
    expect_error(
        test_maf(missing_total),
        "Missing required columns: total_count"
    )

    # Test with multiple missing columns
    missing_multiple <- data.frame(
        other_col = c(1, 2, 3)
    )

    # Verify error message includes all missing columns
    expect_error(
        test_maf(missing_multiple),
        "Missing required columns: ref_count, alt_count, total_count"
    )
})

test_that("test_maf handles empty data frame", {
    # Setup - Empty data frame with correct column structure
    empty_df <- data.frame(
        ref_count = numeric(0),
        alt_count = numeric(0),
        total_count = numeric(0)
    )

    # Execute test
    result <- test_maf(empty_df)

    # Verify result structure
    expect_s3_class(result, "data.frame")

    # Check that result has zero rows but correct columns
    expect_equal(nrow(result), 0)
    expect_true(all(
        c("ref_count", "alt_count", "total_count", "minor_allele_count", "p_val", "adj_p_val") %in% colnames(result)
    ))
})

test_that("test_maf handles single row input", {
    # Setup - Single row data frame
    single_row_df <- data.frame(
        ref_count = 7,
        alt_count = 3,
        total_count = 10
    )

    # Execute test
    result <- test_maf(single_row_df)

    # Verify result structure
    expect_s3_class(result, "data.frame")
    expect_equal(nrow(result), 1)

    # Check calculations for single row
    expect_equal(result$minor_allele_count, 3)
    expect_true(result$p_val >= 0 & result$p_val <= 1)
    # For single row, adjusted p-value should equal original p-value
    expect_equal(result$adj_p_val, result$p_val)
})

test_that("test_maf preserves additional columns", {
    # Setup - Data frame with additional columns beyond required ones
    extended_df <- data.frame(
        snp_id = c("snp_1", "snp_2", "snp_3"),
        ref_count = c(15, 8, 12),
        alt_count = c(5, 7, 3),
        total_count = c(20, 15, 15),
        gene_name = c("GENE1", "GENE2", "GENE3"),
        stringsAsFactors = FALSE
    )

    # Execute test
    result <- test_maf(extended_df)

    # Verify all original columns are preserved
    expect_true("snp_id" %in% colnames(result))
    expect_true("gene_name" %in% colnames(result))
    expect_equal(result$snp_id, extended_df$snp_id)
    expect_equal(result$gene_name, extended_df$gene_name)

    # Check that new columns are added
    expect_true("minor_allele_count" %in% colnames(result))
    expect_true("p_val" %in% colnames(result))
    expect_true("adj_p_val" %in% colnames(result))
})

test_that("test_maf produces consistent results", {
    # Setup - Test data with reproducible results
    test_df <- data.frame(
        ref_count = c(18, 12, 25),
        alt_count = c(2, 8, 5),
        total_count = c(20, 20, 30)
    )

    # Execute test multiple times
    result1 <- test_maf(test_df, p = 0.1)
    result2 <- test_maf(test_df, p = 0.1)

    # Verify results are identical (function should be deterministic)
    expect_equal(result1$minor_allele_count, result2$minor_allele_count)
    expect_equal(result1$p_val, result2$p_val)
    expect_equal(result1$adj_p_val, result2$adj_p_val)
})

test_that("test_maf handles large datasets efficiently", {
    # Setup - Create larger dataset to test performance
    large_df <- data.frame(
        ref_count = sample(1:50, 100, replace = TRUE),
        alt_count = sample(1:20, 100, replace = TRUE),
        total_count = sample(20:70, 100, replace = TRUE)
    )

    # Ensure total_count >= ref_count + alt_count for biological validity
    large_df$total_count <- pmax(large_df$total_count, large_df$ref_count + large_df$alt_count)

    # Execute test and measure that it completes
    # Verify function handles larger datasets without error
    expect_no_error(result <- test_maf(large_df))

    # Check result structure
    expect_s3_class(result, "data.frame")
    expect_equal(nrow(result), 100)

    # Verify all calculations are valid
    expect_true(all(result$p_val >= 0 & result$p_val <= 1))
    expect_true(all(result$adj_p_val >= 0 & result$adj_p_val <= 1))
})

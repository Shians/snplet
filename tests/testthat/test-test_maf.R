# ==============================================================================
# Test Suite: Minor Allele Frequency Testing
# Description: Tests for statistical significance testing of minor allele frequencies
# ==============================================================================

library(testthat)
library(dplyr)

# ------------------------------------------------------------------------------
# Test Helpers
# ------------------------------------------------------------------------------

expect_valid_maf_result <- function(result) {
    expect_true(all(c("minor_allele_count", "p_val", "adj_p_val") %in% colnames(result)))
    expect_true(all(result$p_val >= 0 & result$p_val <= 1))
    expect_true(all(result$adj_p_val >= 0 & result$adj_p_val <= 1))
}

# ==============================================================================

test_that("test_maf returns correct structure and calculations for basic input", {
    test_df <- data.frame(
        ref_count = c(10, 5, 20),
        alt_count = c(2, 8, 1),
        total_count = c(12, 13, 21)
    )

    result <- test_maf(test_df)

    # Verify result is a data frame
    expect_s3_class(result, "data.frame")
    # Check that result has same number of rows as input
    expect_equal(nrow(result), nrow(test_df))
    # Verify all original columns are preserved
    expect_true("ref_count" %in% colnames(result))
    expect_true("alt_count" %in% colnames(result))
    expect_true("total_count" %in% colnames(result))
    # Verify new columns are added
    expect_true("minor_allele_count" %in% colnames(result))
    expect_true("p_val" %in% colnames(result))
    expect_true("adj_p_val" %in% colnames(result))
    # Verify minor_allele_count calculation (minimum of ref_count and alt_count)
    expected_minor_allele <- pmin(test_df$ref_count, test_df$alt_count)
    expect_equal(result$minor_allele_count, expected_minor_allele)
    # Verify p-values are numeric and within valid range [0,1]
    expect_type(result$p_val, "double")
    expect_true(all(result$p_val >= 0 & result$p_val <= 1))
    # Verify adjusted p-values are numeric and within valid range [0,1]
    expect_type(result$adj_p_val, "double")
    expect_true(all(result$adj_p_val >= 0 & result$adj_p_val <= 1))
    # Verify adjusted p-values satisfy BH correction property
    expect_true(all(result$adj_p_val >= result$p_val))
})

test_that("test_maf accepts custom null hypothesis probability", {
    test_df <- data.frame(
        ref_count = c(8, 15),
        alt_count = c(2, 5),
        total_count = c(10, 20)
    )

    result <- test_maf(test_df, p = 0.05)

    # Verify result has required columns
    expect_valid_maf_result(result)
})

test_that("test_maf produces different p-values for different null hypotheses", {
    test_df <- data.frame(
        ref_count = c(8, 15),
        alt_count = c(2, 5),
        total_count = c(10, 20)
    )

    result_p05 <- test_maf(test_df, p = 0.05)
    result_p20 <- test_maf(test_df, p = 0.20)

    # Verify different null hypotheses produce different p-values
    expect_false(identical(result_p05$p_val, result_p20$p_val))
})

test_that("test_maf calculates minor allele counts correctly for edge cases", {
    edge_cases_df <- data.frame(
        ref_count = c(0, 10, 5),
        alt_count = c(10, 0, 5),
        total_count = c(10, 10, 10)
    )

    result <- test_maf(edge_cases_df)

    # Verify minor allele count for zero ref, zero alt, and equal counts
    expect_equal(result$minor_allele_count[1], 0)
    expect_equal(result$minor_allele_count[2], 0)
    expect_equal(result$minor_allele_count[3], 5)
    # Verify p-values are valid
    expect_true(all(result$p_val >= 0 & result$p_val <= 1))
    expect_true(all(result$adj_p_val >= 0 & result$adj_p_val <= 1))
})

test_that("test_maf handles fractional total_count values", {
    fractional_df <- data.frame(
        ref_count = c(8, 12),
        alt_count = c(2, 3),
        total_count = c(10.7, 15.3)
    )

    result <- test_maf(fractional_df)

    # Verify function completes successfully with valid results
    expect_s3_class(result, "data.frame")
    expect_true(all(result$p_val >= 0 & result$p_val <= 1))
    expect_true(all(result$adj_p_val >= 0 & result$adj_p_val <= 1))
})

test_that("test_maf throws error when input is not a data frame", {
    # Verify error for string input
    expect_error(
        test_maf("not_a_dataframe"),
        "is\\(x, \"data\\.frame\"\\) is not TRUE"
    )
    # Verify error for NULL input
    expect_error(
        test_maf(NULL),
        "is\\(x, \"data\\.frame\"\\) is not TRUE"
    )
    # Verify error for list input
    expect_error(
        test_maf(list(a = 1, b = 2)),
        "is\\(x, \"data\\.frame\"\\) is not TRUE"
    )
})

test_that("test_maf throws error when required columns are missing", {
    missing_ref <- data.frame(
        alt_count = c(2, 5),
        total_count = c(10, 15)
    )

    # Verify error when ref_count is missing
    expect_error(
        test_maf(missing_ref),
        "Missing required columns: ref_count"
    )
})

test_that("test_maf throws error when alt_count column is missing", {
    missing_alt <- data.frame(
        ref_count = c(8, 10),
        total_count = c(10, 15)
    )

    expect_error(
        test_maf(missing_alt),
        "Missing required columns: alt_count"
    )
})

test_that("test_maf throws error when total_count column is missing", {
    missing_total <- data.frame(
        ref_count = c(8, 10),
        alt_count = c(2, 5)
    )

    expect_error(
        test_maf(missing_total),
        "Missing required columns: total_count"
    )
})

test_that("test_maf throws error listing all missing columns", {
    missing_multiple <- data.frame(
        other_col = c(1, 2, 3)
    )

    expect_error(
        test_maf(missing_multiple),
        "Missing required columns: ref_count, alt_count, total_count"
    )
})

test_that("test_maf returns empty data frame with correct structure for empty input", {
    empty_df <- data.frame(
        ref_count = numeric(0),
        alt_count = numeric(0),
        total_count = numeric(0)
    )

    result <- test_maf(empty_df)

    # Verify result has zero rows but correct columns
    expect_s3_class(result, "data.frame")
    expect_equal(nrow(result), 0)
    expect_true(all(
        c("ref_count", "alt_count", "total_count", "minor_allele_count", "p_val", "adj_p_val") %in% colnames(result)
    ))
})

test_that("test_maf calculates correct values for single row input", {
    single_row_df <- data.frame(
        ref_count = 7,
        alt_count = 3,
        total_count = 10
    )

    result <- test_maf(single_row_df)

    # Verify result structure
    expect_s3_class(result, "data.frame")
    expect_equal(nrow(result), 1)
    # Verify minor allele count calculation
    expect_equal(result$minor_allele_count, 3)
    # Verify p-value is valid
    expect_true(result$p_val >= 0 & result$p_val <= 1)
    # Verify adjusted p-value equals original p-value for single row
    expect_equal(result$adj_p_val, result$p_val)
})

test_that("test_maf preserves additional columns beyond required ones", {
    extended_df <- data.frame(
        snp_id = c("snp_1", "snp_2", "snp_3"),
        ref_count = c(15, 8, 12),
        alt_count = c(5, 7, 3),
        total_count = c(20, 15, 15),
        gene_name = c("GENE1", "GENE2", "GENE3"),
        stringsAsFactors = FALSE
    )

    result <- test_maf(extended_df)

    # Verify all original columns are preserved with correct values
    expect_true("snp_id" %in% colnames(result))
    expect_true("gene_name" %in% colnames(result))
    expect_equal(result$snp_id, extended_df$snp_id)
    expect_equal(result$gene_name, extended_df$gene_name)
    # Verify new columns are added
    expect_true("minor_allele_count" %in% colnames(result))
    expect_true("p_val" %in% colnames(result))
    expect_true("adj_p_val" %in% colnames(result))
})

test_that("test_maf processes larger datasets without error", {
    n_rows <- 100
    max_ref_count <- 30
    max_alt_count <- 20
    min_total_count <- 60

    large_df <- data.frame(
        ref_count = sample(1:max_ref_count, n_rows, replace = TRUE),
        alt_count = sample(1:max_alt_count, n_rows, replace = TRUE),
        total_count = sample(min_total_count:100, n_rows, replace = TRUE)
    )

    expect_no_error(result <- test_maf(large_df))

    # Verify result structure and validity
    expect_s3_class(result, "data.frame")
    expect_equal(nrow(result), n_rows)
    expect_true(all(result$p_val >= 0 & result$p_val <= 1))
    expect_true(all(result$adj_p_val >= 0 & result$adj_p_val <= 1))
})

test_that("test_maf throws error when ref_count plus alt_count exceeds total_count", {
    inconsistent_df <- data.frame(
        ref_count = c(10, 8),
        alt_count = c(5, 7),
        total_count = c(20, 10)
    )

    # Verify error is thrown for invalid row (row 2: 8 + 7 = 15 > 10)
    expect_error(
        test_maf(inconsistent_df),
        "Invalid data: ref_count \\+ alt_count > total_count in row\\(s\\): 2"
    )
})

test_that("test_maf throws error for single invalid row", {
    single_invalid <- data.frame(
        ref_count = c(10, 15, 8),
        alt_count = c(5, 10, 2),
        total_count = c(15, 20, 12)
    )

    # Verify error is thrown for row 2 (15 + 10 = 25 > 20)
    expect_error(
        test_maf(single_invalid),
        "Invalid data: ref_count \\+ alt_count > total_count in row\\(s\\): 2"
    )
})

test_that("test_maf throws error listing multiple invalid rows", {
    multiple_invalid <- data.frame(
        ref_count = c(10, 15, 8, 12),
        alt_count = c(5, 10, 7, 8),
        total_count = c(15, 20, 10, 25)
    )

    # Verify error lists all invalid rows (2: 15+10=25>20, 3: 8+7=15>10)
    expect_error(
        test_maf(multiple_invalid),
        "Invalid data: ref_count \\+ alt_count > total_count in row\\(s\\): 2, 3"
    )
})

test_that("test_maf throws error with truncated row list for many invalid rows", {
    many_invalid <- data.frame(
        ref_count = rep(10, 10),
        alt_count = rep(10, 10),
        total_count = rep(15, 10)
    )

    # Verify error truncates long list of invalid rows
    expect_error(
        test_maf(many_invalid),
        "Invalid data: ref_count \\+ alt_count > total_count in row\\(s\\): 1, 2, 3, 4, 5, \\.\\.\\."
    )
})

test_that("test_maf accepts data where ref_count plus alt_count equals total_count", {
    exact_match_df <- data.frame(
        ref_count = c(10, 15, 8),
        alt_count = c(5, 5, 2),
        total_count = c(15, 20, 10)
    )

    result <- test_maf(exact_match_df)

    # Verify function accepts exact matches
    expect_s3_class(result, "data.frame")
    expect_valid_maf_result(result)
})

test_that("test_maf accepts data where ref_count plus alt_count is less than total_count", {
    valid_df <- data.frame(
        ref_count = c(10, 15, 8),
        alt_count = c(5, 5, 2),
        total_count = c(20, 25, 15)
    )

    result <- test_maf(valid_df)

    # Verify function accepts valid data
    expect_s3_class(result, "data.frame")
    expect_valid_maf_result(result)
})

test_that("test_maf handles all zero counts in a row", {
    zero_row_df <- data.frame(
        ref_count = c(0, 10, 0),
        alt_count = c(0, 5, 0),
        total_count = c(0, 15, 0)
    )

    result <- test_maf(zero_row_df)

    # Verify function handles all-zero row
    expect_s3_class(result, "data.frame")
    expect_equal(nrow(result), 3)
    expect_equal(result$minor_allele_count[1], 0)
    expect_equal(result$minor_allele_count[3], 0)
})

test_that("test_maf handles NA values in count columns", {
    na_df <- data.frame(
        ref_count = c(10, NA, 15),
        alt_count = c(2, 5, NA),
        total_count = c(12, 13, 18)
    )

    result <- test_maf(na_df)

    # Verify function handles NA values
    expect_s3_class(result, "data.frame")
    expect_true(is.na(result$minor_allele_count[2]))
    expect_true(is.na(result$minor_allele_count[3]))
})

test_that("test_maf handles extreme null hypothesis probability values", {
    test_df <- data.frame(
        ref_count = c(10, 15),
        alt_count = c(2, 5),
        total_count = c(12, 20)
    )

    result_p_low <- test_maf(test_df, p = 0.01)
    result_p_high <- test_maf(test_df, p = 0.99)

    # Verify function handles extreme p values
    expect_valid_maf_result(result_p_low)
    expect_valid_maf_result(result_p_high)
})

test_that("test_maf handles boundary p parameter values", {
    test_df <- data.frame(
        ref_count = c(10, 15),
        alt_count = c(2, 5),
        total_count = c(12, 20)
    )

    result_p_zero <- test_maf(test_df, p = 0)
    result_p_one <- test_maf(test_df, p = 1)

    # Verify function handles boundary p values
    expect_s3_class(result_p_zero, "data.frame")
    expect_s3_class(result_p_one, "data.frame")
    expect_true(all(result_p_zero$p_val >= 0 & result_p_zero$p_val <= 1))
    expect_true(all(result_p_one$p_val >= 0 & result_p_one$p_val <= 1))
})

test_that("test_maf handles negative count values", {
    negative_df <- data.frame(
        ref_count = c(-5, 10),
        alt_count = c(2, -3),
        total_count = c(12, 20)
    )

    result <- test_maf(negative_df)

    # Verify function processes negative values (behavior may vary)
    expect_s3_class(result, "data.frame")
})

test_that("test_maf handles non-numeric values in count columns", {
    non_numeric_df <- data.frame(
        ref_count = c("10", "15"),
        alt_count = c("2", "5"),
        total_count = c("12", "20"),
        stringsAsFactors = FALSE
    )

    # Verify function behavior with character columns
    expect_error(test_maf(non_numeric_df))
})

test_that("test_maf handles very large count values", {
    large_count_df <- data.frame(
        ref_count = c(1e6, 5e6),
        alt_count = c(1e5, 2e5),
        total_count = c(1.1e6, 5.2e6)
    )

    result <- test_maf(large_count_df)

    # Verify function handles large numeric values
    expect_s3_class(result, "data.frame")
    expect_valid_maf_result(result)
})

test_that("test_maf handles mixed valid and edge case rows", {
    mixed_df <- data.frame(
        ref_count = c(10, 0, 100, 5),
        alt_count = c(5, 20, 0, 5),
        total_count = c(15, 20, 100, 10)
    )

    result <- test_maf(mixed_df)

    # Verify function processes mixed data correctly
    expect_s3_class(result, "data.frame")
    expect_equal(nrow(result), 4)
    expect_equal(result$minor_allele_count[1], 5)
    expect_equal(result$minor_allele_count[2], 0)
    expect_equal(result$minor_allele_count[3], 0)
    expect_equal(result$minor_allele_count[4], 5)
    expect_valid_maf_result(result)
})

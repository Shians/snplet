# ==============================================================================
# Test Suite: XCI Escape Testing
# Description: Tests for statistical significance testing of XCI escape fraction
# ==============================================================================

library(testthat)
library(dplyr)

# ------------------------------------------------------------------------------
# Test Helpers
# ------------------------------------------------------------------------------

expect_valid_escape_result <- function(result) {
    # Verify required output columns are present
    expect_true(all(c("coverage", "p_val", "adj_p_val") %in% colnames(result)))
    # Verify p-values fall within the valid [0, 1] range
    expect_true(all(result$p_val >= 0 & result$p_val <= 1))
    # Verify adjusted p-values fall within the valid [0, 1] range
    expect_true(all(result$adj_p_val >= 0 & result$adj_p_val <= 1))
}

# ==============================================================================

test_that("test_escape preserves original columns and adds expected columns", {
    test_df <- data.frame(
        active_count = c(10, 5, 20),
        inactive_count = c(1, 6, 2)
    )

    result <- test_escape(test_df, rho = 0.3)

    # Verify result is a data frame
    expect_s3_class(result, "data.frame")
    # Check that result has same number of rows as input
    expect_equal(nrow(result), nrow(test_df))
    # Verify all original columns are preserved
    expect_true("active_count" %in% colnames(result))
    # Verify inactive_count column is preserved
    expect_true("inactive_count" %in% colnames(result))
    # Verify new columns are added
    expect_true("coverage" %in% colnames(result))
    # Verify p_val column is added
    expect_true("p_val" %in% colnames(result))
    # Verify adj_p_val column is added
    expect_true("adj_p_val" %in% colnames(result))
})

test_that("test_escape calculates coverage as the sum of active and inactive counts", {
    test_df <- data.frame(
        active_count = c(10, 5, 20),
        inactive_count = c(1, 6, 2)
    )

    result <- test_escape(test_df, rho = 0.3)

    # Verify coverage calculation (active_count + inactive_count)
    expect_equal(result$coverage, test_df$active_count + test_df$inactive_count)
})

test_that("test_escape produces valid p-values and adjusted p-values for basic input", {
    test_df <- data.frame(
        active_count = c(10, 5, 20),
        inactive_count = c(1, 6, 2)
    )

    result <- test_escape(test_df, rho = 0.3)

    # Verify p-values are numeric
    expect_type(result$p_val, "double")
    # Check that p-values fall within [0, 1]
    expect_true(all(result$p_val >= 0 & result$p_val <= 1))
    # Verify adjusted p-values are numeric
    expect_type(result$adj_p_val, "double")
    # Check that adjusted p-values fall within [0, 1]
    expect_true(all(result$adj_p_val >= 0 & result$adj_p_val <= 1))
})

test_that("test_escape adjusted p-values satisfy the BH correction property", {
    test_df <- data.frame(
        active_count = c(10, 5, 20),
        inactive_count = c(1, 6, 2)
    )

    result <- test_escape(test_df, rho = 0.3)

    # Verify adjusted p-values satisfy BH correction property
    expect_true(all(result$adj_p_val >= result$p_val))
})

test_that("test_escape matches test_maf/binom_test when rho = 0", {
    # At rho = 0 the beta-binomial null collapses to the binomial one, so
    # test_escape on active/inactive counts should match test_maf on the same
    # counts recast as ref/alt.
    test_df <- data.frame(
        active_count = c(10, 5, 20),
        inactive_count = c(1, 6, 2)
    )

    escape_result <- test_escape(test_df, p = 0.1, rho = 0)
    maf_df <- data.frame(
        ref_count = test_df$active_count,
        alt_count = test_df$inactive_count,
        total_count = test_df$active_count + test_df$inactive_count
    )
    expected <- binom_test(maf_df$alt_count, maf_df$total_count, 0.1, alternative = "greater")

    # Verify beta-binomial p-values collapse to the binomial ones at rho = 0
    expect_equal(escape_result$p_val, expected, tolerance = 1e-6)
})

test_that("test_escape produces different p-values for different rho", {
    test_df <- data.frame(
        active_count = c(10, 15),
        inactive_count = c(4, 5)
    )

    result_low_rho <- test_escape(test_df, rho = 0.01)
    result_high_rho <- test_escape(test_df, rho = 0.5)

    # Verify higher overdispersion produces different (larger) p-values
    expect_false(identical(result_low_rho$p_val, result_high_rho$p_val))
    # Verify overdispersion inflates p-values relative to a near-binomial null
    expect_true(all(result_high_rho$p_val >= result_low_rho$p_val))
})

test_that("test_escape accepts a rho vector recycled across rows", {
    test_df <- data.frame(
        active_count = c(10, 10),
        inactive_count = c(4, 4)
    )

    result <- test_escape(test_df, rho = c(0.05, 0.4))
    result_low <- test_escape(test_df[1, , drop = FALSE], rho = 0.05)
    result_high <- test_escape(test_df[2, , drop = FALSE], rho = 0.4)

    # Verify per-row rho values are applied rather than a single scalar
    expect_equal(result$p_val[1], result_low$p_val)
    # Verify the second row used its own rho value
    expect_equal(result$p_val[2], result_high$p_val)
})

test_that("test_escape accepts a p vector recycled across rows, e.g. a per-donor null", {
    # Mirrors the donor-specific-null usage pattern: p and rho both vary by
    # row (as they would after joining donor_info's xci_median_pi_g/xci_rho).
    test_df <- data.frame(
        donor = c("donor0", "donor1"),
        active_count = c(10, 10),
        inactive_count = c(4, 4)
    )
    donor_null <- data.frame(donor = c("donor0", "donor1"), p = c(0.02, 0.15), rho = c(0.1, 0.4))
    joined <- dplyr::left_join(test_df, donor_null, by = "donor")

    result <- test_escape(joined, p = joined$p, rho = joined$rho)
    result_donor0 <- test_escape(joined[1, , drop = FALSE], p = 0.02, rho = 0.1)
    result_donor1 <- test_escape(joined[2, , drop = FALSE], p = 0.15, rho = 0.4)

    # Verify each row used its own donor-specific null rather than a shared scalar
    expect_equal(result$p_val[1], result_donor0$p_val)
    expect_equal(result$p_val[2], result_donor1$p_val)
    # Verify the two donors' different nulls actually produced different p-values
    expect_false(identical(result$p_val[1], result$p_val[2]))
})

test_that("test_escape accepts custom null hypothesis probability", {
    test_df <- data.frame(
        active_count = c(8, 15),
        inactive_count = c(2, 5)
    )

    result <- test_escape(test_df, p = 0.05, rho = 0.2)

    # Verify result has required columns
    expect_valid_escape_result(result)
})

test_that("test_escape produces different p-values for different null hypotheses", {
    test_df <- data.frame(
        active_count = c(8, 15),
        inactive_count = c(2, 5)
    )

    result_p05 <- test_escape(test_df, p = 0.05, rho = 0.2)
    result_p20 <- test_escape(test_df, p = 0.20, rho = 0.2)

    # Verify different null hypotheses produce different p-values
    expect_false(identical(result_p05$p_val, result_p20$p_val))
})

test_that("test_escape throws error when input is a string", {
    # Verify error for string input
    expect_error(
        test_escape("not_a_dataframe"),
        "is\\(x, \"data\\.frame\"\\) is not TRUE"
    )
})

test_that("test_escape throws error when input is NULL", {
    # Verify error for NULL input
    expect_error(
        test_escape(NULL),
        "is\\(x, \"data\\.frame\"\\) is not TRUE"
    )
})

test_that("test_escape throws error when required columns are missing", {
    missing_active <- data.frame(
        inactive_count = c(2, 5)
    )

    # Verify error when active_count is missing
    expect_error(
        test_escape(missing_active),
        "Missing required columns: active_count"
    )
})

test_that("test_escape throws error listing all missing columns", {
    missing_multiple <- data.frame(
        other_col = c(1, 2, 3)
    )

    # Verify error lists both missing required columns
    expect_error(
        test_escape(missing_multiple),
        "Missing required columns: active_count, inactive_count"
    )
})

test_that("test_escape rejects negative active_count values", {
    negative_active_df <- data.frame(
        active_count = c(-5, 10),
        inactive_count = c(2, 3)
    )

    # Verify error when active_count is negative
    expect_error(
        test_escape(negative_active_df),
        "active_count < 0"
    )
})

test_that("test_escape rejects negative inactive_count values", {
    negative_inactive_df <- data.frame(
        active_count = c(5, 10),
        inactive_count = c(2, -3)
    )

    # Verify error when inactive_count is negative
    expect_error(
        test_escape(negative_inactive_df),
        "inactive_count < 0"
    )
})

test_that("test_escape handles all-zero counts in a row", {
    zero_row_df <- data.frame(
        active_count = c(0, 10, 0),
        inactive_count = c(0, 5, 0)
    )

    result <- test_escape(zero_row_df, rho = 0.2)

    # Verify function handles all-zero row without error
    expect_s3_class(result, "data.frame")
    # Verify all rows are retained despite an all-zero row
    expect_equal(nrow(result), 3)
    # Verify coverage is zero for the all-zero rows
    expect_equal(result$coverage[1], 0)
    expect_equal(result$coverage[3], 0)
    # Verify p-value is 1 for a row with zero inactive reads (no evidence of escape)
    expect_equal(result$p_val[1], 1)
})

test_that("test_escape returns empty data frame with correct structure for empty input", {
    empty_df <- data.frame(
        active_count = numeric(0),
        inactive_count = numeric(0)
    )

    result <- test_escape(empty_df, rho = 0.2)

    # Verify result has zero rows but correct columns
    expect_s3_class(result, "data.frame")
    # Verify result has zero rows for empty input
    expect_equal(nrow(result), 0)
    # Verify all expected columns are present despite empty input
    expect_true(all(
        c("active_count", "inactive_count", "coverage", "p_val", "adj_p_val") %in% colnames(result)
    ))
})

test_that("test_escape preserves additional columns beyond required ones", {
    extended_df <- data.frame(
        gene_name = c("PRKX", "SMC1A", "DDX3X"),
        active_count = c(15, 8, 12),
        inactive_count = c(5, 7, 3),
        stringsAsFactors = FALSE
    )

    result <- test_escape(extended_df, rho = 0.3)

    # Verify additional columns are preserved with correct values
    expect_true("gene_name" %in% colnames(result))
    # Verify gene_name values are unchanged
    expect_equal(result$gene_name, extended_df$gene_name)
})

test_that("test_escape flags a gene with strong biallelic signal as more significant than a tightly silenced one", {
    # A gene close to 50/50 active/inactive looks like escape; one that is
    # almost entirely active-allele reads does not.
    escaping_df <- data.frame(active_count = 12, inactive_count = 10)
    silenced_df <- data.frame(active_count = 20, inactive_count = 1)

    escaping_result <- test_escape(escaping_df, p = 0.1, rho = 0.3)
    silenced_result <- test_escape(silenced_df, p = 0.1, rho = 0.3)

    # Verify the near-biallelic gene gets a smaller (more significant) p-value
    expect_true(escaping_result$p_val < silenced_result$p_val)
})

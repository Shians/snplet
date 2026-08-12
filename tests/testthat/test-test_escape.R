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
    # Verify no method dispatches for a type that is neither counts nor an object
    expect_error(
        test_escape("not_a_dataframe"),
        "unable to find an inherited method"
    )
})

test_that("test_escape throws error when input is NULL", {
    # Verify error for NULL input
    expect_error(
        test_escape(NULL),
        "unable to find an inherited method"
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
    # An uncovered row carries no evidence either way, so it is untested rather
    # than assigned a p-value of 1, which would otherwise pad the BH denominator
    # and cost the testable rows power
    expect_true(is.na(result$p_val[1]))
    expect_true(is.na(result$adj_p_val[1]))
    # Confirm the covered row is corrected as the only test, not one of three
    expect_equal(result$adj_p_val[2], result$p_val[2])
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

# ------------------------------------------------------------------------------
# SNPData method: counts, null and correction all taken from the object
# ------------------------------------------------------------------------------

# Two donors x two genes, one SNP each, with the XCI diagnostics written by hand
# rather than fitted: the point here is what test_escape() does with a stored
# fit, not whether the EM recovers one. X1 carries REF everywhere; cells 1-2 are
# X1-active and cells 3-4 X2-active in both donors.
#
# geneA leaks 2/10 reads from the silenced haplotype in each group, geneB none.
make_escape_fixture <- function(rho = 0.05, median_pi_g = 0.05) {
    ref <- rbind(
        c(8L, 8L, 2L, 2L, 8L, 8L, 2L, 2L),
        c(10L, 10L, 0L, 0L, 10L, 10L, 0L, 0L)
    )
    alt <- rbind(
        c(2L, 2L, 8L, 8L, 2L, 2L, 8L, 8L),
        c(0L, 0L, 10L, 10L, 0L, 0L, 10L, 10L)
    )
    snp_info <- data.frame(
        chrom = "X",
        pos = c(1000L, 2000L),
        ref = "A",
        alt = "G",
        gene_name = c("geneA", "geneB"),
        stringsAsFactors = FALSE
    )
    barcode_info <- data.frame(
        barcode = paste0("cell", 1:8),
        donor = rep(c("donor0", "donor1"), each = 4),
        stringsAsFactors = FALSE
    )

    obj <- SNPData(
        ref_count = Matrix(ref, sparse = TRUE),
        alt_count = Matrix(alt, sparse = TRUE),
        snp_info = snp_info,
        barcode_info = barcode_info
    )

    obj <- add_donor_snp_metadata(
        obj,
        expand.grid(
            snp_id = snp_info(obj)$snp_id,
            donor = c("donor0", "donor1"),
            stringsAsFactors = FALSE
        ) %>%
            dplyr::mutate(xci_informative = TRUE, allele_on_x1 = "REF"),
        join_by = c("snp_id", "donor"),
        overwrite = TRUE
    )
    obj <- add_barcode_metadata(
        obj,
        data.frame(
            cell_id = barcode_info(obj)$cell_id,
            active_x = rep(c("X1", "X1", "X2", "X2"), 2),
            stringsAsFactors = FALSE
        ),
        join_by = "cell_id",
        overwrite = TRUE
    )
    add_donor_metadata(
        obj,
        data.frame(
            donor = c("donor0", "donor1"),
            xci_median_pi_g = median_pi_g,
            xci_rho = rho,
            stringsAsFactors = FALSE
        ),
        join_by = "donor",
        overwrite = TRUE
    )
}

test_that("test_escape() on a SNPData tests one row per donor and gene", {
    obj <- make_escape_fixture()

    result <- test_escape(obj)

    # Confirm the grain is one row per (donor, gene): testing the two active-X
    # groups separately would double the tests and halve each row's coverage
    expect_equal(nrow(result), 4L)
    expect_equal(nrow(dplyr::distinct(result, donor, gene_name)), nrow(result))
    # Verify the counts are pooled over both groups (geneA: 4 cells x 10 reads
    # per donor, which is the two 20-read groups summed)
    expect_equal(dplyr::filter(result, gene_name == "geneA")$coverage, c(40, 40))
    expect_valid_escape_result(result)
})

test_that("test_escape() on a SNPData takes the null from each donor's own fit", {
    lenient <- test_escape(make_escape_fixture(median_pi_g = 0.4))
    strict <- test_escape(make_escape_fixture(median_pi_g = 0.01))

    # Verify the stored per-donor null is used rather than a fixed default: the
    # same counts tested against a higher null escape fraction are less significant
    expect_true(all(lenient$p_val >= strict$p_val))
    expect_false(identical(lenient$p_val, strict$p_val))
})

test_that("test_escape() on a SNPData corrects within each donor by default", {
    obj <- make_escape_fixture()

    within <- test_escape(obj)
    across <- test_escape(obj, by_donor = FALSE)

    # Each donor contributes 2 genes, so correcting within a donor divides by 2
    # while correcting across the table divides by 4; the raw p-values are
    # untouched either way
    expect_equal(within$p_val, across$p_val)
    expect_true(all(across$adj_p_val >= within$adj_p_val))
})

test_that("test_escape() on a SNPData reports where the counts came from", {
    result <- test_escape(make_escape_fixture())

    # Without stored molecule phase the counts come from per-SNP reads, and the
    # column says so: the same call on the same data after add_molecule_phase()
    # counts differently, so the source cannot be left implicit
    expect_equal(unique(result$count_source), "snp")
})

test_that("test_escape() on a SNPData returns a donor with no stored fit untested", {
    obj <- make_escape_fixture()
    obj <- add_donor_metadata(
        obj,
        data.frame(donor = "donor1", xci_rho = NA_real_, stringsAsFactors = FALSE),
        join_by = "donor",
        overwrite = TRUE
    )

    # Verify the donor is named rather than silently producing bad p-values
    expect_warning(result <- test_escape(obj), NA)

    # donor1 has no overdispersion to test against, so its genes come back NA
    # rather than being dropped from the output or tested against a guess
    expect_true(all(is.na(dplyr::filter(result, donor == "donor1")$p_val)))
    # Confirm donor0 is unaffected
    expect_false(any(is.na(dplyr::filter(result, donor == "donor0")$p_val)))
})

test_that("test_escape() on a SNPData requires stored XCI diagnostics", {
    ref <- Matrix(matrix(1L, 2, 2), sparse = TRUE)
    obj <- SNPData(
        ref_count = ref,
        alt_count = ref,
        snp_info = data.frame(chrom = "X", pos = c(1L, 2L), ref = "A", alt = "G"),
        barcode_info = data.frame(barcode = c("c1", "c2"))
    )

    # Verify the error names the step that produces what is missing
    expect_error(test_escape(obj), "Run assign_xci")
})

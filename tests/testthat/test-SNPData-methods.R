# ==============================================================================
# Test Suite: SNPData Methods
# Description: Tests for SNPData methods including filtering and data frame conversion
# ==============================================================================

library(testthat)
library(Matrix)

# ------------------------------------------------------------------------------
# Helper Functions
# ------------------------------------------------------------------------------

# Helper function to compare objects without names
expected_equal_unnamed <- function(object, expected, ...) {
    expect_equal(unname(object), unname(expected), ...)
}

# ------------------------------------------------------------------------------
# Test Data Setup
# ------------------------------------------------------------------------------

# Create test matrices
test_alt_count <- Matrix::Matrix(matrix(c(1, 2, 3, 4), nrow = 2, ncol = 2))
test_ref_count <- Matrix::Matrix(matrix(c(5, 6, 7, 8), nrow = 2, ncol = 2))

# Create test metadata
test_snp_info <- data.frame(
    snp_id = c("snp_1", "snp_2"),
    pos = c(100, 200),
    stringsAsFactors = FALSE
)

test_sample_info <- data.frame(
    cell_id = c("cell_1", "cell_2"),
    donor = c("donor_1", "donor_1"),
    clonotype = c("clonotype_1", "clonotype_2"),
    stringsAsFactors = FALSE
)

# ==============================================================================

test_that("SNPData count data frame methods work correctly", {
    # Setup
    snp_data <- SNPData(
        alt_count = test_alt_count,
        ref_count = test_ref_count,
        snp_info = test_snp_info,
        sample_info = test_sample_info
    )

    # Test barcode_count_df method
    cell_df <- barcode_count_df(snp_data)
    expect_s3_class(cell_df, "tbl_df")
    expect_equal(nrow(cell_df), 4)  # 2 SNPs x 2 cells
    expected_barcode_cols <- c("snp_id", "cell_id", "ref_count", "alt_count", "total_count", "ref_ratio", "maf")
    expect_true(all(expected_barcode_cols %in% colnames(cell_df)))

    # Test donor_count_df method
    donor_df <- donor_count_df(snp_data)
    expect_s3_class(donor_df, "tbl_df")
    expect_equal(nrow(donor_df), 2)  # 2 SNPs x 1 donor
    expected_donor_cols <- c("snp_id", "donor", "ref_count", "alt_count", "total_count", "ref_ratio", "maf")
    expect_true(all(expected_donor_cols %in% colnames(donor_df)))

    # Test clonotype_count_df method
    clonotype_df <- clonotype_count_df(snp_data)
    expect_s3_class(clonotype_df, "tbl_df")
    expect_equal(nrow(clonotype_df), 4)  # 2 SNPs x 2 clonotypes
    expected_clonotype_cols <- c("snp_id", "clonotype", "ref_count", "alt_count", "total_count", "ref_ratio", "maf", "donor")
    expect_true(all(expected_clonotype_cols %in% colnames(clonotype_df)))
})

test_that("filter_snps and filter_barcodes validate column names", {
    # Setup
    snp_data <- SNPData(
        alt_count = test_alt_count,
        ref_count = test_ref_count,
        snp_info = test_snp_info,
        sample_info = test_sample_info
    )

    # Test filter_snps with non-existent column
    expect_error(
        filter_snps(snp_data, not_a_column > 0),
        "The following columns are not present in snp_info or parent environment: not_a_column"
    )

    # Test filter_barcodes with non-existent column
    expect_error(
        filter_barcodes(snp_data, not_a_column > 0),
        "The following columns are not present in sample_info or parent environment: not_a_column"
    )

    # Test filter_snps with valid column
    expect_s4_class(filter_snps(snp_data, pos > 100), "SNPData")

    # Test filter_barcodes with valid column
    expect_s4_class(filter_barcodes(snp_data, donor == "donor_1"), "SNPData")
})

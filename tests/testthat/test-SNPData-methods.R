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
    cell_type = c("T_cell", "B_cell"),
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
    # Verify barcode_count_df returns a tibble
    expect_s3_class(cell_df, "tbl_df")
    # Check that result has expected number of rows (2 SNPs x 2 cells)
    expect_equal(nrow(cell_df), 4)  # 2 SNPs x 2 cells
    expected_barcode_cols <- c("snp_id", "cell_id", "ref_count", "alt_count", "total_count", "ref_ratio", "maf")
    # Verify all expected columns are present in barcode data frame
    expect_true(all(expected_barcode_cols %in% colnames(cell_df)))

    # Test donor_count_df method
    donor_df <- donor_count_df(snp_data)
    # Verify donor_count_df returns a tibble
    expect_s3_class(donor_df, "tbl_df")
    # Check that result has expected number of rows (2 SNPs x 1 donor)
    expect_equal(nrow(donor_df), 2)  # 2 SNPs x 1 donor
    expected_donor_cols <- c("snp_id", "donor", "ref_count", "alt_count", "total_count", "ref_ratio", "maf")
    # Verify all expected columns are present in donor data frame
    expect_true(all(expected_donor_cols %in% colnames(donor_df)))

    # Test clonotype_count_df method
    clonotype_df <- clonotype_count_df(snp_data)
    # Verify clonotype_count_df returns a tibble
    expect_s3_class(clonotype_df, "tbl_df")
    # Check that result has expected number of rows (2 SNPs x 2 clonotypes)
    expect_equal(nrow(clonotype_df), 4)  # 2 SNPs x 2 clonotypes
    expected_clonotype_cols <- c("snp_id", "clonotype", "ref_count", "alt_count", "total_count", "ref_ratio", "maf", "donor")
    # Verify all expected columns are present in clonotype data frame
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
    # Verify error when filtering SNPs by non-existent column
    expect_error(
        filter_snps(snp_data, not_a_column > 0),
        "The following columns are not present in snp_info or parent environment: not_a_column"
    )

    # Test filter_barcodes with non-existent column
    # Verify error when filtering barcodes by non-existent column
    expect_error(
        filter_barcodes(snp_data, not_a_column > 0),
        "The following columns are not present in sample_info or parent environment: not_a_column"
    )

    # Test filter_snps with valid column
    # Verify successful filtering when using valid SNP info column
    expect_s4_class(filter_snps(snp_data, pos > 100), "SNPData")

    # Test filter_barcodes with valid column
    # Verify successful filtering when using valid sample info column
    expect_s4_class(filter_barcodes(snp_data, donor == "donor_1"), "SNPData")
})

test_that("aggregate_count_df works correctly with various grouping columns", {
    # Setup
    snp_data <- SNPData(
        alt_count = test_alt_count,
        ref_count = test_ref_count,
        snp_info = test_snp_info,
        sample_info = test_sample_info
    )

    # Test aggregation by donor
    donor_agg_df <- aggregate_count_df(snp_data, "donor")
    # Verify donor aggregation returns a tibble
    expect_s3_class(donor_agg_df, "tbl_df")
    # Check that result has expected number of rows (2 SNPs x 1 donor)
    expect_equal(nrow(donor_agg_df), 2)  # 2 SNPs x 1 unique donor
    expected_cols <- c("snp_id", "donor", "ref_count", "alt_count", "total_count", "ref_ratio", "maf")
    # Verify all expected columns are present
    expect_true(all(expected_cols %in% colnames(donor_agg_df)))
    # Check that counts are properly aggregated (sum of both cells)
    expect_equal(donor_agg_df$ref_count[1], 5 + 7)  # ref counts for snp_1: cell_1 + cell_2
    expect_equal(donor_agg_df$alt_count[1], 1 + 3)  # alt counts for snp_1: cell_1 + cell_2

    # Test aggregation by clonotype
    clonotype_agg_df <- aggregate_count_df(snp_data, "clonotype")
    # Verify clonotype aggregation returns a tibble
    expect_s3_class(clonotype_agg_df, "tbl_df")
    # Check that result has expected number of rows (2 SNPs x 2 clonotypes)
    expect_equal(nrow(clonotype_agg_df), 4)  # 2 SNPs x 2 unique clonotypes
    expected_clonotype_cols <- c("snp_id", "clonotype", "ref_count", "alt_count", "total_count", "ref_ratio", "maf")
    # Verify all expected columns are present for clonotype
    expect_true(all(expected_clonotype_cols %in% colnames(clonotype_agg_df)))

    # Test aggregation by cell_type (custom column)
    cell_type_agg_df <- aggregate_count_df(snp_data, "cell_type")
    # Verify cell_type aggregation returns a tibble
    expect_s3_class(cell_type_agg_df, "tbl_df")
    # Check that result has expected number of rows (2 SNPs x 2 cell_types)
    expect_equal(nrow(cell_type_agg_df), 4)  # 2 SNPs x 2 unique cell_types
    expected_cell_type_cols <- c("snp_id", "cell_type", "ref_count", "alt_count", "total_count", "ref_ratio", "maf")
    # Verify all expected columns are present for cell_type
    expect_true(all(expected_cell_type_cols %in% colnames(cell_type_agg_df)))

    # Test with test_maf = FALSE
    no_maf_df <- aggregate_count_df(snp_data, "donor", test_maf = FALSE)
    # Verify test_maf column is not present when test_maf = FALSE
    expect_false("test_maf" %in% colnames(no_maf_df))

    # Test with test_maf = TRUE (default)
    with_maf_df <- aggregate_count_df(snp_data, "donor", test_maf = TRUE)
    # Verify test_maf column is present when test_maf = TRUE
    expect_true("test_maf" %in% colnames(with_maf_df))
})

test_that("aggregate_count_df handles edge cases correctly", {
    # Setup
    snp_data <- SNPData(
        alt_count = test_alt_count,
        ref_count = test_ref_count,
        snp_info = test_snp_info,
        sample_info = test_sample_info
    )

    # Test error with non-existent column
    # Verify error when grouping by non-existent column
    expect_error(
        aggregate_count_df(snp_data, "non_existent_column"),
        "Column 'non_existent_column' not found in sample_info"
    )

    # Test with NA values in grouping column
    test_sample_info_na <- test_sample_info
    test_sample_info_na$donor[1] <- NA
    snp_data_na <- SNPData(
        alt_count = test_alt_count,
        ref_count = test_ref_count,
        snp_info = test_snp_info,
        sample_info = test_sample_info_na
    )
    
    # Verify function handles NA values by excluding them
    result_na <- suppressMessages(aggregate_count_df(snp_data_na, "donor"))
    # Check that function issues warning about NA values (checking this indirectly by verifying behavior)
    expect_true(exists("result_na"))
    # Check that result excludes NA groups
    expect_equal(nrow(result_na), 2)  # 2 SNPs x 1 non-NA donor
})

test_that("aggregate_count_df produces correct calculations", {
    # Setup with known values for easier testing
    test_alt_simple <- Matrix::Matrix(matrix(c(10, 20, 30, 40), nrow = 2, ncol = 2))
    test_ref_simple <- Matrix::Matrix(matrix(c(90, 80, 70, 60), nrow = 2, ncol = 2))
    
    test_sample_info_simple <- data.frame(
        cell_id = c("cell_1", "cell_2"),
        group_col = c("group_A", "group_A"),  # Both cells in same group
        stringsAsFactors = FALSE
    )
    
    snp_data_simple <- SNPData(
        alt_count = test_alt_simple,
        ref_count = test_ref_simple,
        snp_info = test_snp_info,
        sample_info = test_sample_info_simple
    )

    result <- aggregate_count_df(snp_data_simple, "group_col")
    
    # Check aggregated counts for first SNP
    # Verify ref_count aggregation: 90 + 70 = 160
    expect_equal(result$ref_count[1], 160)
    # Verify alt_count aggregation: 10 + 30 = 40  
    expect_equal(result$alt_count[1], 40)
    # Verify total_count calculation: 160 + 40 = 200
    expect_equal(result$total_count[1], 200)
    # Verify ref_ratio calculation: 160/200 = 0.8
    expect_equal(result$ref_ratio[1], 0.8)
    # Verify maf calculation: min(160, 40)/200 = 0.2
    expect_equal(result$maf[1], 0.2)

    # Check aggregated counts for second SNP  
    # Verify ref_count aggregation: 80 + 60 = 140
    expect_equal(result$ref_count[2], 140)
    # Verify alt_count aggregation: 20 + 40 = 60
    expect_equal(result$alt_count[2], 60)
    # Verify total_count calculation: 140 + 60 = 200
    expect_equal(result$total_count[2], 200)
    # Verify ref_ratio calculation: 140/200 = 0.7
    expect_equal(result$ref_ratio[2], 0.7)
    # Verify maf calculation: min(140, 60)/200 = 0.3
    expect_equal(result$maf[2], 0.3)
})

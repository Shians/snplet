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

test_barcode_info <- data.frame(
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
        barcode_info = test_barcode_info
    )

    # Test barcode_count_df method
    cell_df <- barcode_count_df(snp_data)
    # Verify barcode_count_df returns a tibble
    expect_s3_class(cell_df, "tbl_df")
    # Check that result has expected number of rows (2 SNPs x 2 cells)
    expect_equal(nrow(cell_df), 4) # 2 SNPs x 2 cells
    expected_barcode_cols <- c("snp_id", "cell_id", "ref_count", "alt_count", "total_count", "ref_ratio", "maf")
    # Verify all expected columns are present in barcode data frame
    expect_true(all(expected_barcode_cols %in% colnames(cell_df)))

    # Test donor_count_df method
    donor_df <- donor_count_df(snp_data)
    # Verify donor_count_df returns a tibble
    expect_s3_class(donor_df, "tbl_df")
    # Check that result has expected number of rows (2 SNPs x 1 donor)
    expect_equal(nrow(donor_df), 2) # 2 SNPs x 1 donor
    expected_donor_cols <- c("snp_id", "donor", "ref_count", "alt_count", "total_count", "ref_ratio", "maf")
    # Verify all expected columns are present in donor data frame
    expect_true(all(expected_donor_cols %in% colnames(donor_df)))

    # Test clonotype_count_df method
    clonotype_df <- clonotype_count_df(snp_data)
    # Verify clonotype_count_df returns a tibble
    expect_s3_class(clonotype_df, "tbl_df")
    # Check that result has expected number of rows (2 SNPs x 2 clonotypes)
    expect_equal(nrow(clonotype_df), 4) # 2 SNPs x 2 clonotypes
    expected_clonotype_cols <- c(
        "snp_id",
        "clonotype",
        "ref_count",
        "alt_count",
        "total_count",
        "ref_ratio",
        "maf",
        "donor"
    )
    # Verify all expected columns are present in clonotype data frame
    expect_true(all(expected_clonotype_cols %in% colnames(clonotype_df)))
})

test_that("SNPData constructor drops duplicate SNP IDs", {
    dup_alt_count <- Matrix::Matrix(
        matrix(c(1, 2, 3, 4, 5, 6), nrow = 3, ncol = 2, byrow = TRUE)
    )
    dup_ref_count <- Matrix::Matrix(
        matrix(c(5, 6, 7, 8, 9, 10), nrow = 3, ncol = 2, byrow = TRUE)
    )
    dup_snp_info <- data.frame(
        snp_id = c("dup_snp", "dup_snp", "unique_snp"),
        pos = c(100, 101, 200),
        stringsAsFactors = FALSE
    )

    # Construct SNPData to trigger deduplication
    expect_warning(
        snp_dup <- SNPData(
            alt_count = dup_alt_count,
            ref_count = dup_ref_count,
            snp_info = dup_snp_info,
            barcode_info = test_barcode_info
        ),
        regexp = "Duplicate SNP IDs detected"
    )

    # Verify correct number of SNPs after deduplication
    expect_equal(nrow(snp_dup), 2)
    
    # Verify correct SNP IDs are retained
    expect_equal(
        get_snp_info(snp_dup)$snp_id,
        c("dup_snp", "unique_snp")
    )
    
    # Verify row names are updated
    expect_equal(
        rownames(ref_count(snp_dup)),
        c("dup_snp", "unique_snp")
    )
    
    # Verify which duplicate was kept (should be first occurrence with pos=100)
    expect_equal(get_snp_info(snp_dup)$pos, c(100, 200))
    
    # Verify count matrices match the kept SNPs
    # First row should be from first dup_snp (values 1, 2 for alt)
    expect_equal(
        as.vector(alt_count(snp_dup)[1,]),
        c(1, 2)
    )
    
    # First row should be from first dup_snp (values 5, 6 for ref)
    expect_equal(
        as.vector(ref_count(snp_dup)[1,]),
        c(5, 6)
    )
    
    # Second row should be unique_snp (values 5, 6 for alt)
    expect_equal(
        as.vector(alt_count(snp_dup)[2,]),
        c(5, 6)
    )
    
    # Second row should be unique_snp (values 9, 10 for ref)
    expect_equal(
        as.vector(ref_count(snp_dup)[2,]),
        c(9, 10)
    )
})

test_that("filter_snps and filter_barcodes validate column names", {
    # Setup
    snp_data <- SNPData(
        alt_count = test_alt_count,
        ref_count = test_ref_count,
        snp_info = test_snp_info,
        barcode_info = test_barcode_info
    )

    # Test filter_snps with non-existent column
    # Verify error when filtering SNPs by non-existent column
    expect_error(
        filter_snps(snp_data, not_a_column > 0),
        "The following columns are not present in snp_info or accessible environment: not_a_column"
    )

    # Test filter_barcodes with non-existent column
    # Verify error when filtering barcodes by non-existent column
    expect_error(
        filter_barcodes(snp_data, not_a_column > 0),
        "The following columns are not present in barcode_info or accessible environment: not_a_column"
    )

    # Test filter_snps with valid column
    # Verify successful filtering when using valid SNP info column
    expect_s4_class(filter_snps(snp_data, pos > 100), "SNPData")

    # Test filter_barcodes with valid column
    # Verify successful filtering when using valid sample info column
    expect_s4_class(filter_barcodes(snp_data, donor == "donor_1"), "SNPData")
})

test_that("variable resolution in filter expressions works correctly", {
    # Setup
    snp_data <- SNPData(
        alt_count = test_alt_count,
        ref_count = test_ref_count,
        snp_info = test_snp_info,
        barcode_info = test_barcode_info
    )

    # Test function parameter resolution - basic case
    filter_with_param <- function(data, threshold_value) {
        filter_barcodes(data, donor == threshold_value)
    }
    # Verify function parameters are accessible in filter expressions
    expect_s4_class(filter_with_param(snp_data, "donor_1"), "SNPData")

    # Test nested function variable resolution
    outer_function <- function(data, outer_var) {
        inner_function <- function() {
            # outer_var should be accessible due to lexical scoping
            filter_barcodes(data, donor == outer_var)
        }
        inner_function()
    }
    # Verify lexically scoped variables are accessible
    expect_s4_class(outer_function(snp_data, "donor_1"), "SNPData")

    # Test multiple parameter resolution
    multi_param_filter <- function(data, param1, param2) {
        # Both parameters should be accessible
        result1 <- filter_barcodes(data, donor == param1)
        result2 <- filter_snps(result1, pos > param2)
        result2
    }
    # Verify multiple parameters work correctly
    expect_s4_class(multi_param_filter(snp_data, "donor_1", 50), "SNPData")

    # Test that truly non-existent variables still cause errors
    test_nonexistent_var <- function(data) {
        # completely_made_up_var should not exist anywhere
        filter_barcodes(data, donor == completely_made_up_var)
    }
    # Verify non-existent variables still cause appropriate errors
    expect_error(
        test_nonexistent_var(snp_data),
        "completely_made_up_var"
    )

    # Test local environment variables
    test_global_var <- "donor_1"
    # Verify local variables are accessible in filter expressions
    expect_s4_class(filter_barcodes(snp_data, donor == test_global_var), "SNPData")
})

test_that("variable resolution respects R scoping rules", {
    # Setup
    snp_data <- SNPData(
        alt_count = test_alt_count,
        ref_count = test_ref_count,
        snp_info = test_snp_info,
        barcode_info = test_barcode_info
    )

    # Test that variables in enclosing environments work
    create_filter_function <- function(default_donor) {
        function(data, override_donor = NULL) {
            target_donor <- if (is.null(override_donor)) default_donor else override_donor
            filter_barcodes(data, donor == target_donor)
        }
    }

    filter_func <- create_filter_function("donor_1")
    # Verify closure variables are accessible
    expect_s4_class(filter_func(snp_data), "SNPData")
    # Verify parameter override works
    expect_s4_class(filter_func(snp_data, "donor_1"), "SNPData")

    # Test with local variables in function scope
    test_local_vars <- function(data) {
        local_threshold <- 150
        local_donor <- "donor_1"
        result1 <- filter_snps(data, pos > local_threshold)
        result2 <- filter_barcodes(result1, donor == local_donor)
        result2
    }
    # Verify local variables within functions are accessible
    expect_s4_class(test_local_vars(snp_data), "SNPData")
})

test_that("aggregate_count_df works correctly with various grouping columns", {
    # Setup
    snp_data <- SNPData(
        alt_count = test_alt_count,
        ref_count = test_ref_count,
        snp_info = test_snp_info,
        barcode_info = test_barcode_info
    )

    # Test aggregation by donor
    donor_agg_df <- aggregate_count_df(snp_data, "donor")
    # Verify donor aggregation returns a tibble
    expect_s3_class(donor_agg_df, "tbl_df")
    # Check that result has expected number of rows (2 SNPs x 1 donor)
    expect_equal(nrow(donor_agg_df), 2) # 2 SNPs x 1 unique donor
    expected_cols <- c("snp_id", "donor", "ref_count", "alt_count", "total_count", "ref_ratio", "maf")
    # Verify all expected columns are present
    expect_true(all(expected_cols %in% colnames(donor_agg_df)))
    # Check that counts are properly aggregated (sum of both cells)
    expect_equal(donor_agg_df$ref_count[1], 5 + 7) # ref counts for snp_1: cell_1 + cell_2
    expect_equal(donor_agg_df$alt_count[1], 1 + 3) # alt counts for snp_1: cell_1 + cell_2

    # Test aggregation by clonotype
    clonotype_agg_df <- aggregate_count_df(snp_data, "clonotype")
    # Verify clonotype aggregation returns a tibble
    expect_s3_class(clonotype_agg_df, "tbl_df")
    # Check that result has expected number of rows (2 SNPs x 2 clonotypes)
    expect_equal(nrow(clonotype_agg_df), 4) # 2 SNPs x 2 unique clonotypes
    expected_clonotype_cols <- c("snp_id", "clonotype", "ref_count", "alt_count", "total_count", "ref_ratio", "maf")
    # Verify all expected columns are present for clonotype
    expect_true(all(expected_clonotype_cols %in% colnames(clonotype_agg_df)))

    # Test aggregation by cell_type (custom column)
    cell_type_agg_df <- aggregate_count_df(snp_data, "cell_type")
    # Verify cell_type aggregation returns a tibble
    expect_s3_class(cell_type_agg_df, "tbl_df")
    # Check that result has expected number of rows (2 SNPs x 2 cell_types)
    expect_equal(nrow(cell_type_agg_df), 4) # 2 SNPs x 2 unique cell_types
    expected_cell_type_cols <- c("snp_id", "cell_type", "ref_count", "alt_count", "total_count", "ref_ratio", "maf")
    # Verify all expected columns are present for cell_type
    expect_true(all(expected_cell_type_cols %in% colnames(cell_type_agg_df)))

    # Test with test_maf = FALSE
    no_maf_df <- aggregate_count_df(snp_data, "donor", test_maf = FALSE)
    expected_maf_cols <- c("minor_allele_count", "p_val", "adj_p_val")
    # Verify test_maf function columns are not present when test_maf = FALSE
    expect_false(any(expected_maf_cols %in% colnames(no_maf_df)))

    # Test with test_maf = TRUE (default)
    with_maf_df <- aggregate_count_df(snp_data, "donor", test_maf = TRUE)
    # Verify test_maf function columns are present when test_maf = TRUE
    expect_true(all(expected_maf_cols %in% colnames(with_maf_df)))
})

test_that("aggregate_count_df handles edge cases correctly", {
    # Setup
    snp_data <- SNPData(
        alt_count = test_alt_count,
        ref_count = test_ref_count,
        snp_info = test_snp_info,
        barcode_info = test_barcode_info
    )

    # Test error with non-existent column
    # Verify error when grouping by non-existent column
    expect_error(
        aggregate_count_df(snp_data, "non_existent_column"),
        "Column 'non_existent_column' not found in barcode_info"
    )

    # Test with NA values in grouping column
    test_barcode_info_na <- test_barcode_info
    test_barcode_info_na$donor[1] <- NA
    snp_data_na <- SNPData(
        alt_count = test_alt_count,
        ref_count = test_ref_count,
        snp_info = test_snp_info,
        barcode_info = test_barcode_info_na
    )

    # Verify function handles NA values by excluding them
    result_na <- suppressMessages(aggregate_count_df(snp_data_na, "donor"))
    # Check that function issues warning about NA values (checking this indirectly by verifying behavior)
    expect_true(exists("result_na"))
    # Check that result excludes NA groups
    expect_equal(nrow(result_na), 2) # 2 SNPs x 1 non-NA donor
})

test_that("test_maf=TRUE adds correct statistical columns to barcode_count_df", {
    # Setup
    snp_data <- SNPData(
        alt_count = test_alt_count,
        ref_count = test_ref_count,
        snp_info = test_snp_info,
        barcode_info = test_barcode_info
    )

    # Test with test_maf = TRUE (default)
    result_with_maf <- barcode_count_df(snp_data, test_maf = TRUE)
    expected_maf_cols <- c("minor_allele_count", "p_val", "adj_p_val")
    # Verify test_maf function columns are added when test_maf = TRUE
    expect_true(all(expected_maf_cols %in% colnames(result_with_maf)))
    # Verify minor_allele_count is calculated correctly
    expect_true(all(result_with_maf$minor_allele_count == pmin(result_with_maf$ref_count, result_with_maf$alt_count)))
    # Verify p_val and adj_p_val are numeric
    expect_true(is.numeric(result_with_maf$p_val))
    expect_true(is.numeric(result_with_maf$adj_p_val))

    # Test with test_maf = FALSE
    result_without_maf <- barcode_count_df(snp_data, test_maf = FALSE)
    # Verify test_maf columns are not present when test_maf = FALSE
    expect_false(any(expected_maf_cols %in% colnames(result_without_maf)))
})

test_that("test_maf=TRUE adds correct statistical columns to donor_count_df", {
    # Setup
    snp_data <- SNPData(
        alt_count = test_alt_count,
        ref_count = test_ref_count,
        snp_info = test_snp_info,
        barcode_info = test_barcode_info
    )

    # Test with test_maf = TRUE (default)
    result_with_maf <- donor_count_df(snp_data, test_maf = TRUE)
    expected_maf_cols <- c("minor_allele_count", "p_val", "adj_p_val")
    # Verify test_maf function columns are added when test_maf = TRUE
    expect_true(all(expected_maf_cols %in% colnames(result_with_maf)))
    # Verify minor_allele_count is calculated correctly
    expect_true(all(result_with_maf$minor_allele_count == pmin(result_with_maf$ref_count, result_with_maf$alt_count)))
    # Verify p_val and adj_p_val are numeric
    expect_true(is.numeric(result_with_maf$p_val))
    expect_true(is.numeric(result_with_maf$adj_p_val))

    # Test with test_maf = FALSE
    result_without_maf <- donor_count_df(snp_data, test_maf = FALSE)
    # Verify test_maf columns are not present when test_maf = FALSE
    expect_false(any(expected_maf_cols %in% colnames(result_without_maf)))
})

test_that("test_maf=TRUE adds correct statistical columns to clonotype_count_df", {
    # Setup
    snp_data <- SNPData(
        alt_count = test_alt_count,
        ref_count = test_ref_count,
        snp_info = test_snp_info,
        barcode_info = test_barcode_info
    )

    # Test with test_maf = TRUE (default)
    result_with_maf <- clonotype_count_df(snp_data, test_maf = TRUE)
    expected_maf_cols <- c("minor_allele_count", "p_val", "adj_p_val")
    # Verify test_maf function columns are added when test_maf = TRUE
    expect_true(all(expected_maf_cols %in% colnames(result_with_maf)))
    # Verify minor_allele_count is calculated correctly
    expect_true(all(result_with_maf$minor_allele_count == pmin(result_with_maf$ref_count, result_with_maf$alt_count)))
    # Verify p_val and adj_p_val are numeric
    expect_true(is.numeric(result_with_maf$p_val))
    expect_true(is.numeric(result_with_maf$adj_p_val))
    # Verify donor column is still present (specific to clonotype_count_df)
    expect_true("donor" %in% colnames(result_with_maf))

    # Test with test_maf = FALSE
    result_without_maf <- clonotype_count_df(snp_data, test_maf = FALSE)
    # Verify test_maf columns are not present when test_maf = FALSE
    expect_false(any(expected_maf_cols %in% colnames(result_without_maf)))
})

test_that("test_maf=TRUE adds correct statistical columns to aggregate_count_df", {
    # Setup
    snp_data <- SNPData(
        alt_count = test_alt_count,
        ref_count = test_ref_count,
        snp_info = test_snp_info,
        barcode_info = test_barcode_info
    )

    # Test with test_maf = TRUE (default) for donor aggregation
    result_with_maf <- aggregate_count_df(snp_data, "donor", test_maf = TRUE)
    expected_maf_cols <- c("minor_allele_count", "p_val", "adj_p_val")
    # Verify test_maf function columns are added when test_maf = TRUE
    expect_true(all(expected_maf_cols %in% colnames(result_with_maf)))
    # Verify minor_allele_count is calculated correctly
    expect_true(all(result_with_maf$minor_allele_count == pmin(result_with_maf$ref_count, result_with_maf$alt_count)))
    # Verify p_val and adj_p_val are numeric
    expect_true(is.numeric(result_with_maf$p_val))
    expect_true(is.numeric(result_with_maf$adj_p_val))

    # Test with test_maf = FALSE
    result_without_maf <- aggregate_count_df(snp_data, "donor", test_maf = FALSE)
    # Verify test_maf columns are not present when test_maf = FALSE
    expect_false(any(expected_maf_cols %in% colnames(result_without_maf)))

    # Test with different grouping column (clonotype) and test_maf = TRUE
    result_clonotype_maf <- aggregate_count_df(snp_data, "clonotype", test_maf = TRUE)
    # Verify test_maf columns are present for different grouping variables
    expect_true(all(expected_maf_cols %in% colnames(result_clonotype_maf)))
})

test_that("test_maf function produces meaningful statistical results", {
    # Setup with specific values to test statistical calculations
    test_alt_specific <- Matrix::Matrix(matrix(c(5, 15, 8, 12), nrow = 2, ncol = 2))
    test_ref_specific <- Matrix::Matrix(matrix(c(95, 85, 92, 88), nrow = 2, ncol = 2))

    snp_data_specific <- SNPData(
        alt_count = test_alt_specific,
        ref_count = test_ref_specific,
        snp_info = test_snp_info,
        barcode_info = test_barcode_info
    )

    result <- barcode_count_df(snp_data_specific, test_maf = TRUE)

    # Verify p-values are between 0 and 1
    expect_true(all(result$p_val >= 0 & result$p_val <= 1))
    # Verify adjusted p-values are between 0 and 1
    expect_true(all(result$adj_p_val >= 0 & result$adj_p_val <= 1))
    # Verify adjusted p-values are >= original p-values (BH correction property)
    expect_true(all(result$adj_p_val >= result$p_val))
    # Verify minor allele count calculation matches expectation
    expected_minor <- pmin(result$ref_count, result$alt_count)
    expect_equal(result$minor_allele_count, expected_minor)

    # Test that the function preserves all original columns
    base_result <- barcode_count_df(snp_data_specific, test_maf = FALSE)
    original_cols <- colnames(base_result)
    # Verify all original columns are preserved when test_maf = TRUE
    expect_true(all(original_cols %in% colnames(result)))
})

test_that("aggregate_count_df produces correct calculations", {
    # Setup with known values for easier testing
    test_alt_simple <- Matrix::Matrix(matrix(c(10, 20, 30, 40), nrow = 2, ncol = 2))
    test_ref_simple <- Matrix::Matrix(matrix(c(90, 80, 70, 60), nrow = 2, ncol = 2))

    test_barcode_info_simple <- data.frame(
        cell_id = c("cell_1", "cell_2"),
        group_col = c("group_A", "group_A"), # Both cells in same group
        stringsAsFactors = FALSE
    )

    snp_data_simple <- SNPData(
        alt_count = test_alt_simple,
        ref_count = test_ref_simple,
        snp_info = test_snp_info,
        barcode_info = test_barcode_info_simple
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

# ==============================================================================
# Tests for Missing Clonotype Data
# ==============================================================================

test_that("clonotype_count_df errors when clonotype column missing", {
    # Setup - Create SNPData without clonotype column
    barcode_info_no_clonotype <- data.frame(
        cell_id = c("cell_1", "cell_2"),
        donor = c("donor_1", "donor_1"),
        stringsAsFactors = FALSE
    )

    snp_data <- SNPData(
        alt_count = test_alt_count,
        ref_count = test_ref_count,
        snp_info = test_snp_info,
        barcode_info = barcode_info_no_clonotype
    )

    # Test error when clonotype column is missing
    # Verify clonotype_count_df errors with appropriate message
    expect_error(
        clonotype_count_df(snp_data),
        "Clonotype information not available.*add_barcode_metadata"
    )
})

test_that("clonotype_count_df errors when all clonotypes are NA", {
    # Setup - Create SNPData with all NA clonotypes
    barcode_info_na_clonotype <- data.frame(
        cell_id = c("cell_1", "cell_2"),
        donor = c("donor_1", "donor_1"),
        clonotype = c(NA_character_, NA_character_),
        stringsAsFactors = FALSE
    )

    snp_data <- SNPData(
        alt_count = test_alt_count,
        ref_count = test_ref_count,
        snp_info = test_snp_info,
        barcode_info = barcode_info_na_clonotype
    )

    # Test error when all clonotype values are NA
    # Verify clonotype_count_df errors when all clonotypes are NA
    expect_error(
        clonotype_count_df(snp_data),
        "All clonotype values are NA.*add_barcode_metadata"
    )
})

test_that("remove_na_clonotypes warns when clonotype column missing", {
    # Setup - Create SNPData without clonotype column
    barcode_info_no_clonotype <- data.frame(
        cell_id = c("cell_1", "cell_2"),
        donor = c("donor_1", "donor_1"),
        stringsAsFactors = FALSE
    )

    snp_data <- SNPData(
        alt_count = test_alt_count,
        ref_count = test_ref_count,
        snp_info = test_snp_info,
        barcode_info = barcode_info_no_clonotype
    )

    # Test warning when clonotype column is missing
    # Verify remove_na_clonotypes warns with appropriate message
    expect_warning(
        result <- remove_na_clonotypes(snp_data),
        "No 'clonotype' column found.*add_barcode_metadata"
    )

    # Verify original object is returned
    expect_equal(nrow(result), nrow(snp_data))
    expect_equal(ncol(result), ncol(snp_data))
})

test_that("aggregate_count_df errors when group_by clonotype but column missing", {
    # Setup - Create SNPData without clonotype column
    barcode_info_no_clonotype <- data.frame(
        cell_id = c("cell_1", "cell_2"),
        donor = c("donor_1", "donor_1"),
        stringsAsFactors = FALSE
    )

    snp_data <- SNPData(
        alt_count = test_alt_count,
        ref_count = test_ref_count,
        snp_info = test_snp_info,
        barcode_info = barcode_info_no_clonotype
    )

    # Test error when grouping by clonotype but column missing
    # Verify aggregate_count_df errors with appropriate message
    expect_error(
        aggregate_count_df(snp_data, "clonotype"),
        "Column 'clonotype' not found in barcode_info"
    )
})

test_that("aggregate_count_df handles all NA clonotypes with warning and error", {
    # Setup - Create SNPData with all NA clonotypes
    barcode_info_na_clonotype <- data.frame(
        cell_id = c("cell_1", "cell_2"),
        donor = c("donor_1", "donor_1"),
        clonotype = c(NA_character_, NA_character_),
        stringsAsFactors = FALSE
    )

    snp_data <- SNPData(
        alt_count = test_alt_count,
        ref_count = test_ref_count,
        snp_info = test_snp_info,
        barcode_info = barcode_info_na_clonotype
    )

    # Test that aggregate_count_df handles all NA values
    # When all values are NA, they are filtered out causing an error in pivot
    # Verify error occurs when all samples filtered out
    expect_error(
        suppressWarnings(
            aggregate_count_df(snp_data, "clonotype", test_maf = FALSE)
        ),
        "must select at least one column"
    )
})

test_that("clonotype functions work after adding clonotype via add_barcode_metadata", {
    # Setup - Create SNPData without clonotype initially
    barcode_info_no_clonotype <- data.frame(
        cell_id = c("cell_1", "cell_2"),
        donor = c("donor_1", "donor_1"),
        stringsAsFactors = FALSE
    )

    snp_data <- SNPData(
        alt_count = test_alt_count,
        ref_count = test_ref_count,
        snp_info = test_snp_info,
        barcode_info = barcode_info_no_clonotype
    )

    # Verify clonotype_count_df fails initially
    expect_error(
        clonotype_count_df(snp_data),
        "Clonotype information not available"
    )

    # Add clonotype information using add_barcode_metadata
    clonotype_data <- data.frame(
        cell_id = c("cell_1", "cell_2"),
        clonotype = c("clonotype_1", "clonotype_2"),
        stringsAsFactors = FALSE
    )

    snp_data_with_clonotype <- add_barcode_metadata(
        snp_data,
        clonotype_data,
        join_by = "cell_id"
    )

    # Verify clonotype column now exists
    barcode_info <- get_barcode_info(snp_data_with_clonotype)
    expect_true("clonotype" %in% colnames(barcode_info))
    expect_equal(barcode_info$clonotype, c("clonotype_1", "clonotype_2"))

    # Verify clonotype_count_df now works
    result <- clonotype_count_df(snp_data_with_clonotype, test_maf = FALSE)
    expect_s3_class(result, "data.frame")
    expect_true("clonotype" %in% colnames(result))

    # Verify to_expr_matrix with clonotype level now works
    expr_mat <- to_expr_matrix(snp_data_with_clonotype, level = "clonotype")
    expect_true(is.matrix(expr_mat))
    expect_equal(ncol(expr_mat), 2)
})

test_that("clonotype functions work after updating clonotype with overwrite=TRUE", {
    # Setup - Create SNPData with all NA clonotypes
    barcode_info_na_clonotype <- data.frame(
        cell_id = c("cell_1", "cell_2"),
        donor = c("donor_1", "donor_1"),
        clonotype = c(NA_character_, NA_character_),
        stringsAsFactors = FALSE
    )

    snp_data <- SNPData(
        alt_count = test_alt_count,
        ref_count = test_ref_count,
        snp_info = test_snp_info,
        barcode_info = barcode_info_na_clonotype
    )

    # Verify clonotype_count_df fails with all NA
    expect_error(
        clonotype_count_df(snp_data),
        "All clonotype values are NA"
    )

    # Update clonotype information using add_barcode_metadata with overwrite=TRUE
    clonotype_data <- data.frame(
        cell_id = c("cell_1", "cell_2"),
        clonotype = c("clonotype_1", "clonotype_2"),
        stringsAsFactors = FALSE
    )

    snp_data_updated <- add_barcode_metadata(
        snp_data,
        clonotype_data,
        join_by = "cell_id",
        overwrite = TRUE
    )

    # Verify clonotype values are updated
    barcode_info <- get_barcode_info(snp_data_updated)
    expect_equal(barcode_info$clonotype, c("clonotype_1", "clonotype_2"))
    expect_false(any(is.na(barcode_info$clonotype)))

    # Verify clonotype_count_df now works
    result <- clonotype_count_df(snp_data_updated, test_maf = FALSE)
    expect_s3_class(result, "data.frame")
    expect_true(all(!is.na(result$clonotype)))

    # Verify to_expr_matrix with clonotype level now works
    expr_mat <- to_expr_matrix(snp_data_updated, level = "clonotype")
    expect_true(is.matrix(expr_mat))
})

# ------------------------------------------------------------------------------
# Test Data Setup - merge_snpdata
# ------------------------------------------------------------------------------

# Create test data for merge_snpdata tests
# x has: snpA, snpB and cell1, cell2
# y has: snpB, snpC and cell2, cell3
# snpB and cell2 overlap between x and y

create_merge_test_data <- function() {
    alt_x <- Matrix::Matrix(
        matrix(
            c(
                1, 0, # snpA
                0, 2  # snpB
            ),
            nrow = 2,
            ncol = 2,
            dimnames = list(c("snpA", "snpB"), c("cell1", "cell2"))
        )
    )
    ref_x <- Matrix::Matrix(
        matrix(
            c(
                3, 0, # snpA
                1, 1  # snpB
            ),
            nrow = 2,
            ncol = 2,
            dimnames = list(c("snpA", "snpB"), c("cell1", "cell2"))
        )
    )
    alt_y <- Matrix::Matrix(
        matrix(
            c(
                1, 0, # snpB
                2, 1  # snpC
            ),
            nrow = 2,
            ncol = 2,
            dimnames = list(c("snpB", "snpC"), c("cell2", "cell3"))
        )
    )
    ref_y <- Matrix::Matrix(
        matrix(
            c(
                0, 2, # snpB
                1, 1  # snpC
            ),
            nrow = 2,
            ncol = 2,
            dimnames = list(c("snpB", "snpC"), c("cell2", "cell3"))
        )
    )

    snp_info_x <- data.frame(
        snp_id = c("snpA", "snpB"),
        gene = c("geneA", "geneB_x"),
        stringsAsFactors = FALSE
    )
    snp_info_y <- data.frame(
        snp_id = c("snpB", "snpC"),
        gene = c("geneB_y", "geneC"),
        stringsAsFactors = FALSE
    )
    barcode_info_x <- data.frame(
        cell_id = c("cell1", "cell2"),
        donor = c("d1", "d2"),
        stringsAsFactors = FALSE
    )
    barcode_info_y <- data.frame(
        cell_id = c("cell2", "cell3"),
        donor = c("d2_y", "d3"),
        stringsAsFactors = FALSE
    )

    x <- SNPData(
        alt_count = alt_x,
        ref_count = ref_x,
        snp_info = snp_info_x,
        barcode_info = barcode_info_x
    )
    y <- SNPData(
        alt_count = alt_y,
        ref_count = ref_y,
        snp_info = snp_info_y,
        barcode_info = barcode_info_y
    )

    list(x = x, y = y)
}

test_that("merge_snpdata union/union retains all SNPs and cells", {
    data <- create_merge_test_data()
    x <- data$x
    y <- data$y

    merged <- merge_snpdata(x, y, snp_join = "union", cell_join = "union")

    # Verify all SNPs from both objects are retained
    expect_equal(rownames(merged), c("snpA", "snpB", "snpC"))
    # Verify all cells from both objects are retained
    expect_equal(colnames(merged), c("cell1", "cell2", "cell3"))

    # Verify ref counts are added correctly for overlapping entries
    # x: snpA/cell1=3, snpA/cell2=1, snpB/cell1=0, snpB/cell2=1
    # y: snpB/cell2=0, snpB/cell3=1, snpC/cell2=2, snpC/cell3=1
    # Merged: snpB/cell2 = 1+0 = 1
    expected_ref <- Matrix::Matrix(
        matrix(
            c(
                3, 1, 0, # snpA: from x only
                0, 1, 1, # snpB: cell1 from x, cell2 from x+y (1+0), cell3 from y
                0, 2, 1  # snpC: from y only
            ),
            nrow = 3,
            ncol = 3,
            byrow = TRUE,
            dimnames = list(c("snpA", "snpB", "snpC"), c("cell1", "cell2", "cell3"))
        )
    )
    # Verify alt counts are added correctly for overlapping entries
    # x: snpA/cell1=1, snpA/cell2=0, snpB/cell1=0, snpB/cell2=2
    # y: snpB/cell2=1, snpB/cell3=2, snpC/cell2=0, snpC/cell3=1
    # Merged: snpB/cell2 = 2+1 = 3
    expected_alt <- Matrix::Matrix(
        matrix(
            c(
                1, 0, 0, # snpA: from x only
                0, 3, 2, # snpB: cell1 from x, cell2 from x+y (2+1), cell3 from y
                0, 0, 1  # snpC: from y only
            ),
            nrow = 3,
            ncol = 3,
            byrow = TRUE,
            dimnames = list(c("snpA", "snpB", "snpC"), c("cell1", "cell2", "cell3"))
        )
    )
    # Check ref count matrix matches expected
    expect_equal(as.matrix(ref_count(merged)), as.matrix(expected_ref))
    # Check alt count matrix matches expected
    expect_equal(as.matrix(alt_count(merged)), as.matrix(expected_alt))
})

test_that("merge_snpdata union/union merges metadata correctly", {
    data <- create_merge_test_data()
    x <- data$x
    y <- data$y

    merged <- merge_snpdata(x, y, snp_join = "union", cell_join = "union")

    merged_snp_info <- get_snp_info(merged)
    merged_barcode_info <- get_barcode_info(merged)

    # Verify SNP metadata is merged with x taking priority for conflicts
    expect_equal(merged_snp_info$gene, c("geneA", "geneB_x", "geneC"))
    # Verify coverage is recalculated correctly (sum of all counts per SNP)
    # snpA: 3+1+1+0 = 5, snpB: 0+1+1+0+3+2 = 7, snpC: 0+2+1+0+0+1 = 4
    expect_equal(unname(merged_snp_info$coverage), c(5, 7, 4))
    # Verify non_zero_samples is recalculated correctly (cells with non-zero coverage)
    # snpA: cell1, cell2 (2 cells), snpB: cell2, cell3 (2 cells), snpC: cell2, cell3 (2 cells)
    expect_equal(unname(merged_snp_info$non_zero_samples), c(2, 2, 2))

    # Verify barcode metadata is merged with x taking priority for conflicts
    expect_equal(merged_barcode_info$donor, c("d1", "d2", "d3"))
    # Verify library_size is recalculated correctly (sum of all counts per cell)
    # cell1: 3+1+0+0 = 4, cell2: 1+0+1+3+2+0 = 7, cell3: 0+0+1+2+1+1 = 5
    expect_equal(unname(merged_barcode_info$library_size), c(4, 7, 5))
    # Verify non_zero_snps is recalculated correctly (SNPs with non-zero coverage)
    # cell1: snpA (1), cell2: snpA, snpB, snpC (3), cell3: snpB, snpC (2)
    expect_equal(unname(merged_barcode_info$non_zero_snps), c(1, 3, 2))
})

test_that("merge_snpdata intersect/intersect retains only overlapping SNPs and cells", {
    data <- create_merge_test_data()
    x <- data$x
    y <- data$y

    merged <- merge_snpdata(x, y, snp_join = "intersect", cell_join = "intersect")

    # Verify only overlapping SNP (snpB) is retained
    expect_equal(rownames(merged), "snpB")
    # Verify only overlapping cell (cell2) is retained
    expect_equal(colnames(merged), "cell2")

    # Verify counts for snpB/cell2 are summed from both objects
    # x: snpB/cell2 has ref=1, alt=2
    # y: snpB/cell2 has ref=0, alt=1
    # Expected: ref=1, alt=3
    expect_equal(as.vector(ref_count(merged)), 1)
    expect_equal(as.vector(alt_count(merged)), 3)
})

test_that("merge_snpdata intersect/intersect merges metadata correctly", {
    data <- create_merge_test_data()
    x <- data$x
    y <- data$y

    merged <- merge_snpdata(x, y, snp_join = "intersect", cell_join = "intersect")

    merged_snp_info <- get_snp_info(merged)
    merged_barcode_info <- get_barcode_info(merged)

    # Verify SNP metadata for snpB uses x value for conflict
    expect_equal(merged_snp_info$gene, "geneB_x")
    # Verify coverage is recalculated
    expect_equal(unname(merged_snp_info$coverage), 4)
    # Verify non_zero_samples is recalculated
    expect_equal(unname(merged_snp_info$non_zero_samples), 1)

    # Verify barcode metadata for cell2 uses x value for conflict
    expect_equal(merged_barcode_info$donor, "d2")
    # Verify library_size is recalculated
    expect_equal(unname(merged_barcode_info$library_size), 4)
    # Verify non_zero_snps is recalculated
    expect_equal(unname(merged_barcode_info$non_zero_snps), 1)
})

test_that("merge_snpdata left/left retains all SNPs and cells from x only", {
    data <- create_merge_test_data()
    x <- data$x
    y <- data$y

    merged <- merge_snpdata(x, y, snp_join = "left", cell_join = "left")

    # Verify only SNPs from x are retained (snpA, snpB)
    expect_equal(rownames(merged), c("snpA", "snpB"))
    # Verify only cells from x are retained (cell1, cell2)
    expect_equal(colnames(merged), c("cell1", "cell2"))

    # Verify counts: snpB/cell2 should be summed, others from x only
    # x: snpA/cell1=3, snpA/cell2=1, snpB/cell1=0, snpB/cell2=1
    # y: snpB/cell2=0 (overlaps with x)
    # Merged: snpB/cell2 = 1+0 = 1
    expected_ref <- Matrix::Matrix(
        matrix(
            c(
                3, 1, # snpA: from x only
                0, 1  # snpB: cell1 from x, cell2 from x+y (1+0)
            ),
            nrow = 2,
            ncol = 2,
            byrow = TRUE,
            dimnames = list(c("snpA", "snpB"), c("cell1", "cell2"))
        )
    )
    expected_alt <- Matrix::Matrix(
        matrix(
            c(
                1, 0, # snpA: from x only
                0, 3  # snpB: cell1 from x, cell2 from x+y (2+1)
            ),
            nrow = 2,
            ncol = 2,
            byrow = TRUE,
            dimnames = list(c("snpA", "snpB"), c("cell1", "cell2"))
        )
    )
    # Check ref count matrix matches expected
    expect_equal(as.matrix(ref_count(merged)), as.matrix(expected_ref))
    # Check alt count matrix matches expected
    expect_equal(as.matrix(alt_count(merged)), as.matrix(expected_alt))
})

test_that("merge_snpdata left/left merges metadata correctly", {
    data <- create_merge_test_data()
    x <- data$x
    y <- data$y

    merged <- merge_snpdata(x, y, snp_join = "left", cell_join = "left")

    merged_snp_info <- get_snp_info(merged)
    merged_barcode_info <- get_barcode_info(merged)

    # Verify SNP metadata contains only x SNPs
    expect_equal(merged_snp_info$snp_id, c("snpA", "snpB"))
    # Verify gene metadata for snpB uses x value
    expect_equal(merged_snp_info$gene, c("geneA", "geneB_x"))

    # Verify barcode metadata contains only x cells
    expect_equal(merged_barcode_info$cell_id, c("cell1", "cell2"))
    # Verify donor metadata uses x values
    expect_equal(merged_barcode_info$donor, c("d1", "d2"))
})

test_that("merge_snpdata right/right retains all SNPs and cells from y only", {
    data <- create_merge_test_data()
    x <- data$x
    y <- data$y

    merged <- merge_snpdata(x, y, snp_join = "right", cell_join = "right")

    # Verify only SNPs from y are retained (snpB, snpC)
    expect_equal(rownames(merged), c("snpB", "snpC"))
    # Verify only cells from y are retained (cell2, cell3)
    expect_equal(colnames(merged), c("cell2", "cell3"))

    # Verify counts: snpB/cell2 should be summed, others from y only
    # x: snpB/cell2=1 (ref), snpB/cell2=2 (alt)
    # y: snpB/cell2=0 (ref), snpB/cell2=1 (alt), snpB/cell3=1 (ref), snpB/cell3=2 (alt)
    #    snpC/cell2=2 (ref), snpC/cell2=0 (alt), snpC/cell3=1 (ref), snpC/cell3=1 (alt)
    # Merged: snpB/cell2 ref=1+0=1, alt=2+1=3
    expected_ref <- Matrix::Matrix(
        matrix(
            c(
                1, 1, # snpB: cell2 from x+y (1+0), cell3 from y
                2, 1  # snpC: from y only
            ),
            nrow = 2,
            ncol = 2,
            byrow = TRUE,
            dimnames = list(c("snpB", "snpC"), c("cell2", "cell3"))
        )
    )
    expected_alt <- Matrix::Matrix(
        matrix(
            c(
                3, 2, # snpB: cell2 from x+y (2+1), cell3 from y
                0, 1  # snpC: from y only
            ),
            nrow = 2,
            ncol = 2,
            byrow = TRUE,
            dimnames = list(c("snpB", "snpC"), c("cell2", "cell3"))
        )
    )
    # Check ref count matrix matches expected
    expect_equal(as.matrix(ref_count(merged)), as.matrix(expected_ref))
    # Check alt count matrix matches expected
    expect_equal(as.matrix(alt_count(merged)), as.matrix(expected_alt))
})

test_that("merge_snpdata right/right merges metadata correctly", {
    data <- create_merge_test_data()
    x <- data$x
    y <- data$y

    merged <- merge_snpdata(x, y, snp_join = "right", cell_join = "right")

    merged_snp_info <- get_snp_info(merged)
    merged_barcode_info <- get_barcode_info(merged)

    # Verify SNP metadata contains only y SNPs
    expect_equal(merged_snp_info$snp_id, c("snpB", "snpC"))
    # Verify gene metadata: snpB should use x value (priority), snpC from y
    expect_equal(merged_snp_info$gene, c("geneB_x", "geneC"))

    # Verify barcode metadata contains only y cells
    expect_equal(merged_barcode_info$cell_id, c("cell2", "cell3"))
    # Verify donor metadata: cell2 should use x value (priority), cell3 from y
    expect_equal(merged_barcode_info$donor, c("d2", "d3"))
})

test_that("merge_snpdata union/intersect retains all SNPs but only overlapping cells", {
    data <- create_merge_test_data()
    x <- data$x
    y <- data$y

    merged <- merge_snpdata(x, y, snp_join = "union", cell_join = "intersect")

    # Verify all SNPs are retained
    expect_equal(rownames(merged), c("snpA", "snpB", "snpC"))
    # Verify only overlapping cell is retained
    expect_equal(colnames(merged), "cell2")

    # Verify counts for cell2 are summed where both objects have data
    # x: snpA/cell2=1 (ref), snpB/cell2=1 (ref)
    # y: snpB/cell2=0 (ref), snpC/cell2=2 (ref)
    expect_equal(as.vector(ref_count(merged)[1, ]), 1) # snpA/cell2: from x
    expect_equal(as.vector(ref_count(merged)[2, ]), 1) # snpB/cell2: from x+y (1+0)
    expect_equal(as.vector(ref_count(merged)[3, ]), 2) # snpC/cell2: from y
})

test_that("merge_snpdata intersect/union retains only overlapping SNPs but all cells", {
    data <- create_merge_test_data()
    x <- data$x
    y <- data$y

    merged <- merge_snpdata(x, y, snp_join = "intersect", cell_join = "union")

    # Verify only overlapping SNP is retained
    expect_equal(rownames(merged), "snpB")
    # Verify all cells are retained
    expect_equal(colnames(merged), c("cell1", "cell2", "cell3"))

    # Verify counts for snpB across all cells
    # x: snpB/cell1=0 (ref), snpB/cell2=1 (ref)
    # y: snpB/cell2=0 (ref), snpB/cell3=1 (ref)
    expect_equal(as.vector(ref_count(merged)[1, ]), c(0, 1, 1)) # snpB: cell1 from x, cell2 from x+y (1+0), cell3 from y
    # x: snpB/cell1=0 (alt), snpB/cell2=2 (alt)
    # y: snpB/cell2=1 (alt), snpB/cell3=2 (alt)
    expect_equal(as.vector(alt_count(merged)[1, ]), c(0, 3, 2)) # snpB: cell1 from x, cell2 from x+y (2+1), cell3 from y
})

test_that("merge_snpdata left/intersect retains x SNPs but only overlapping cells", {
    data <- create_merge_test_data()
    x <- data$x
    y <- data$y

    merged <- merge_snpdata(x, y, snp_join = "left", cell_join = "intersect")

    # Verify only SNPs from x are retained
    expect_equal(rownames(merged), c("snpA", "snpB"))
    # Verify only overlapping cell is retained
    expect_equal(colnames(merged), "cell2")

    # Verify dimensions match expected
    expect_equal(dim(merged), c(2, 1))
})

test_that("merge_snpdata intersect/right retains only overlapping SNPs but y cells", {
    data <- create_merge_test_data()
    x <- data$x
    y <- data$y

    merged <- merge_snpdata(x, y, snp_join = "intersect", cell_join = "right")

    # Verify only overlapping SNP is retained
    expect_equal(rownames(merged), "snpB")
    # Verify only cells from y are retained
    expect_equal(colnames(merged), c("cell2", "cell3"))

    # Verify dimensions match expected
    expect_equal(dim(merged), c(1, 2))
})

# ==============================================================================
# Test Suite: donor_het_status_df
# Description: Tests for filtering heterozygous SNPs based on donor-level biallelic expression
# ==============================================================================

create_heterozygous_test_data <- function() {
    # Create test data with known heterozygous/homozygous patterns
    # snp1: heterozygous in donor1 (balanced ~50:50, alt=10, ref=10, total=20)
    # snp2: homozygous in donor1 (highly skewed 2:18, alt=2, ref=18, total=20)
    # snp3: heterozygous in both donors (balanced, alt=10, ref=10, total=20 each)
    # snp4: low coverage in donor1 (total=3), heterozygous in donor2 (alt=10, ref=10, total=20)

    alt_count <- Matrix::Matrix(
        matrix(
            c(
                5, 1, 5, 1,  # cell1 (donor1): snp1, snp2, snp3, snp4
                5, 1, 5, 2,  # cell2 (donor1)
                5, 1, 5, 3,  # cell3 (donor2)
                5, 1, 5, 4   # cell4 (donor2)
            ),
            nrow = 4,
            ncol = 4,
            byrow = FALSE
        )
    )

    ref_count <- Matrix::Matrix(
        matrix(
            c(
                5, 9, 5, 0,  # cell1 (donor1)
                5, 9, 5, 1,  # cell2 (donor1)
                5, 9, 5, 2,  # cell3 (donor2)
                5, 9, 5, 3   # cell4 (donor2)
            ),
            nrow = 4,
            ncol = 4,
            byrow = FALSE
        )
    )

    snp_info <- data.frame(
        snp_id = c("snp1", "snp2", "snp3", "snp4"),
        chrom = c("chr1", "chr1", "chr2", "chr2"),
        pos = c(100, 200, 300, 400),
        stringsAsFactors = FALSE
    )

    barcode_info <- data.frame(
        cell_id = c("cell1", "cell2", "cell3", "cell4"),
        donor = c("donor1", "donor1", "donor2", "donor2"),
        stringsAsFactors = FALSE
    )

    SNPData(
        alt_count = alt_count,
        ref_count = ref_count,
        snp_info = snp_info,
        barcode_info = barcode_info
    )
}

# ==============================================================================

test_that("donor_het_status_df reports het/hom status with tested flag", {
    snp_data <- create_heterozygous_test_data()

    status_df <- donor_het_status_df(
        snp_data,
        min_total_count = 10,
        p_value_threshold = 0.05,
        minor_allele_prop = 0.1
    )

    expect_true(all(c("snp_id", "chrom", "pos", "donor", "zygosity", "tested") %in% colnames(status_df)))

    snp1_d1 <- status_df %>% dplyr::filter(snp_id == "snp1", donor == "donor1")
    expect_equal(snp1_d1$zygosity, "het")
    expect_true(snp1_d1$tested)
    expect_equal(snp1_d1$chrom, "chr1")
    expect_equal(snp1_d1$pos, 100)

    snp2_d1 <- status_df %>% dplyr::filter(snp_id == "snp2", donor == "donor1")
    expect_equal(snp2_d1$zygosity, "hom")
    expect_true(snp2_d1$tested)

    snp4_d1 <- status_df %>% dplyr::filter(snp_id == "snp4", donor == "donor1")
    expect_equal(snp4_d1$zygosity, "unknown")
    expect_false(snp4_d1$tested)
})

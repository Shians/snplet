# ==============================================================================
# Test Suite: Data Cleaning Functions
# Description: Tests for functions that remove doublets, NA genes, and NA clonotypes
# ==============================================================================

library(testthat)
library(Matrix)

# ------------------------------------------------------------------------------
# Test Data Setup
# ------------------------------------------------------------------------------

# Create test matrices
test_alt_count <- Matrix::Matrix(matrix(c(1, 2, 3, 4, 5, 6), nrow = 3, ncol = 2))
test_ref_count <- Matrix::Matrix(matrix(c(7, 8, 9, 10, 11, 12), nrow = 3, ncol = 2))

# Create test SNP info with some NA gene names
test_snp_info <- data.frame(
    snp_id = c("snp_1", "snp_2", "snp_3"),
    pos = c(100, 200, 300),
    gene_name = c("GENE1", NA, "GENE3"),
    stringsAsFactors = FALSE
)

# Create test sample info with doublets and NA clonotypes
test_sample_info <- data.frame(
    cell_id = c("cell_1", "cell_2"),
    donor = c("donor_1", "doublet"),
    clonotype = c("clonotype_1", NA),
    stringsAsFactors = FALSE
)

# Create test SNPData object
create_test_snpdata <- function() {
    SNPData(
        alt_count = test_alt_count,
        ref_count = test_ref_count,
        snp_info = test_snp_info,
        sample_info = test_sample_info
    )
}

# ==============================================================================

test_that("remove_doublets works correctly", {
    # Setup - Create test data with doublets
    snp_data <- create_test_snpdata()

    # Test basic doublet removal
    filtered_data <- remove_doublets(snp_data)

    # Verify filtered data is still SNPData object
    expect_s4_class(filtered_data, "SNPData")

    # Check that doublet cell was removed (only 1 cell should remain)
    expect_equal(ncol(filtered_data), 1)

    # Verify remaining cell is not the doublet
    remaining_sample_info <- get_barcode_info(filtered_data)
    expect_equal(remaining_sample_info$donor, "donor_1")
    expect_false("doublet" %in% remaining_sample_info$donor)

    # Check that matrices were properly subsetted
    expect_equal(dim(alt_count(filtered_data)), c(3, 1))
    expect_equal(dim(ref_count(filtered_data)), c(3, 1))

    # Verify SNP info unchanged (row subsetting not affected)
    expect_equal(nrow(get_snp_info(filtered_data)), 3)
})

test_that("remove_doublets handles missing donor column", {
    # Setup - Create SNPData without donor column
    sample_info_no_donor <- data.frame(
        cell_id = c("cell_1", "cell_2"),
        clonotype = c("clonotype_1", "clonotype_2"),
        stringsAsFactors = FALSE
    )

    snp_data_no_donor <- SNPData(
        alt_count = test_alt_count,
        ref_count = test_ref_count,
        snp_info = test_snp_info,
        sample_info = sample_info_no_donor
    )

    # Test with missing donor column
    # Verify warning is issued when donor column is missing
    expect_warning(
        result <- remove_doublets(snp_data_no_donor),
        "No 'donor' column found in sample_info, returning original object"
    )

    # Check that original object is returned unchanged
    expect_equal(ncol(result), ncol(snp_data_no_donor))
    expect_equal(nrow(result), nrow(snp_data_no_donor))
})

test_that("remove_doublets handles NA values correctly", {
    # Setup - Create data with NA donor values
    sample_info_with_na <- data.frame(
        cell_id = c("cell_1", "cell_2", "cell_3"),
        donor = c("donor_1", NA, "donor_2"),
        clonotype = c("clonotype_1", "clonotype_2", "clonotype_3"),
        stringsAsFactors = FALSE
    )

    alt_count_3col <- Matrix::Matrix(matrix(1:9, nrow = 3, ncol = 3))
    ref_count_3col <- Matrix::Matrix(matrix(10:18, nrow = 3, ncol = 3))

    snp_data_with_na <- SNPData(
        alt_count = alt_count_3col,
        ref_count = ref_count_3col,
        snp_info = test_snp_info,
        sample_info = sample_info_with_na
    )

    # Test with drop_na = TRUE (default)
    filtered_drop_na <- remove_doublets(snp_data_with_na, drop_na = TRUE)
    # Verify NA donor cells are removed when drop_na = TRUE
    expect_equal(ncol(filtered_drop_na), 2)  # Should remove cell with NA donor

    # Test with drop_na = FALSE
    filtered_keep_na <- remove_doublets(snp_data_with_na, drop_na = FALSE)
    # Verify NA donor cells are kept when drop_na = FALSE
    expect_equal(ncol(filtered_keep_na), 3)  # Should keep all cells
})

test_that("remove_doublets validates input", {
    # Test with non-SNPData object
    # Verify error when input is not SNPData object
    expect_error(
        remove_doublets("not_snpdata"),
        "Input must be a SNPData object"
    )

    # Verify error with NULL input
    expect_error(
        remove_doublets(NULL),
        "Input must be a SNPData object"
    )
})

test_that("remove_na_genes works correctly", {
    # Setup - Use test data with NA gene names
    snp_data <- create_test_snpdata()

    # Test basic NA gene removal
    filtered_data <- remove_na_genes(snp_data)

    # Verify filtered data is still SNPData object
    expect_s4_class(filtered_data, "SNPData")

    # Check that SNP with NA gene_name was removed (2 SNPs should remain)
    expect_equal(nrow(filtered_data), 2)

    # Verify remaining SNPs don't have NA gene names
    remaining_snp_info <- get_snp_info(filtered_data)
    expect_false(any(is.na(remaining_snp_info$gene_name)))
    expect_equal(remaining_snp_info$gene_name, c("GENE1", "GENE3"))

    # Check that matrices were properly subsetted
    expect_equal(dim(alt_count(filtered_data)), c(2, 2))
    expect_equal(dim(ref_count(filtered_data)), c(2, 2))

    # Verify sample info unchanged (column subsetting not affected)
    expect_equal(ncol(get_barcode_info(filtered_data)), ncol(get_barcode_info(snp_data)))
})

test_that("remove_na_genes handles missing gene column", {
    # Setup - Create SNPData without gene_name column
    snp_info_no_gene <- data.frame(
        snp_id = c("snp_1", "snp_2", "snp_3"),
        pos = c(100, 200, 300),
        stringsAsFactors = FALSE
    )

    snp_data_no_gene <- SNPData(
        alt_count = test_alt_count,
        ref_count = test_ref_count,
        snp_info = snp_info_no_gene,
        sample_info = test_sample_info
    )

    # Test with missing gene column
    # Verify warning is issued when gene column is missing
    expect_warning(
        result <- remove_na_genes(snp_data_no_gene),
        "No 'gene_name' column found in snp_info, returning original object"
    )

    # Check that original object is returned unchanged
    expect_equal(nrow(result), nrow(snp_data_no_gene))
    expect_equal(ncol(result), ncol(snp_data_no_gene))
})

test_that("remove_na_genes handles custom gene column", {
    # Setup - Create data with custom gene column name
    snp_info_custom <- data.frame(
        snp_id = c("snp_1", "snp_2", "snp_3"),
        pos = c(100, 200, 300),
        custom_gene = c("GENE1", NA, "GENE3"),
        stringsAsFactors = FALSE
    )

    snp_data_custom <- SNPData(
        alt_count = test_alt_count,
        ref_count = test_ref_count,
        snp_info = snp_info_custom,
        sample_info = test_sample_info
    )

    # Test with custom gene column name
    filtered_data <- remove_na_genes(snp_data_custom, gene_col = "custom_gene")

    # Verify filtering worked with custom column
    expect_equal(nrow(filtered_data), 2)
    remaining_snp_info <- get_snp_info(filtered_data)
    expect_false(any(is.na(remaining_snp_info$custom_gene)))
})

test_that("remove_na_genes validates input", {
    # Test with non-SNPData object
    # Verify error when input is not SNPData object
    expect_error(
        remove_na_genes("not_snpdata"),
        "Input must be a SNPData object"
    )
})

test_that("remove_na_clonotypes works correctly", {
    # Setup - Use test data with NA clonotypes
    snp_data <- create_test_snpdata()

    # Test basic NA clonotype removal
    filtered_data <- remove_na_clonotypes(snp_data)

    # Verify filtered data is still SNPData object
    expect_s4_class(filtered_data, "SNPData")

    # Check that cell with NA clonotype was removed (1 cell should remain)
    expect_equal(ncol(filtered_data), 1)

    # Verify remaining cell doesn't have NA clonotype
    remaining_sample_info <- get_barcode_info(filtered_data)
    expect_false(any(is.na(remaining_sample_info$clonotype)))
    expect_equal(remaining_sample_info$clonotype, "clonotype_1")

    # Check that matrices were properly subsetted
    expect_equal(dim(alt_count(filtered_data)), c(3, 1))
    expect_equal(dim(ref_count(filtered_data)), c(3, 1))

    # Verify SNP info unchanged (row subsetting not affected)
    expect_equal(nrow(get_snp_info(filtered_data)), 3)
})

test_that("remove_na_clonotypes handles missing clonotype column", {
    # Setup - Create SNPData without clonotype column
    sample_info_no_clonotype <- data.frame(
        cell_id = c("cell_1", "cell_2"),
        donor = c("donor_1", "donor_2"),
        stringsAsFactors = FALSE
    )

    snp_data_no_clonotype <- SNPData(
        alt_count = test_alt_count,
        ref_count = test_ref_count,
        snp_info = test_snp_info,
        sample_info = sample_info_no_clonotype
    )

    # Test with missing clonotype column
    # Verify warning is issued when clonotype column is missing
    expect_warning(
        result <- remove_na_clonotypes(snp_data_no_clonotype),
        "No 'clonotype' column found in sample_info, returning original object"
    )

    # Check that original object is returned unchanged
    expect_equal(ncol(result), ncol(snp_data_no_clonotype))
    expect_equal(nrow(result), nrow(snp_data_no_clonotype))
})

test_that("remove_na_clonotypes handles custom clonotype column", {
    # Setup - Create data with custom clonotype column name
    sample_info_custom <- data.frame(
        cell_id = c("cell_1", "cell_2"),
        donor = c("donor_1", "donor_2"),
        custom_clonotype = c("clonotype_1", NA),
        stringsAsFactors = FALSE
    )

    snp_data_custom <- SNPData(
        alt_count = test_alt_count,
        ref_count = test_ref_count,
        snp_info = test_snp_info,
        sample_info = sample_info_custom
    )

    # Test with custom clonotype column name
    filtered_data <- remove_na_clonotypes(snp_data_custom, clonotype_col = "custom_clonotype")

    # Verify filtering worked with custom column
    expect_equal(ncol(filtered_data), 1)
    remaining_sample_info <- get_barcode_info(filtered_data)
    expect_false(any(is.na(remaining_sample_info$custom_clonotype)))
})

test_that("remove_na_clonotypes validates input", {
    # Test with non-SNPData object
    # Verify error when input is not SNPData object
    expect_error(
        remove_na_clonotypes("not_snpdata"),
        "Input must be a SNPData object"
    )
})

test_that("all remove functions handle edge cases", {
    # Setup - Create data with all cells/SNPs meeting removal criteria

    # Test remove_doublets with all doublets
    all_doublets_info <- data.frame(
        cell_id = c("cell_1", "cell_2"),
        donor = c("doublet", "doublet"),
        clonotype = c("clonotype_1", "clonotype_2"),
        stringsAsFactors = FALSE
    )

    snp_data_all_doublets <- SNPData(
        alt_count = test_alt_count,
        ref_count = test_ref_count,
        snp_info = test_snp_info,
        sample_info = all_doublets_info
    )

    # Verify all doublets are removed, leaving empty object
    all_doublets_removed <- remove_doublets(snp_data_all_doublets)
    expect_equal(ncol(all_doublets_removed), 0)

    # Test remove_na_genes with all NA genes
    all_na_genes_info <- data.frame(
        snp_id = c("snp_1", "snp_2", "snp_3"),
        pos = c(100, 200, 300),
        gene_name = c(NA, NA, NA),
        stringsAsFactors = FALSE
    )

    snp_data_all_na_genes <- SNPData(
        alt_count = test_alt_count,
        ref_count = test_ref_count,
        snp_info = all_na_genes_info,
        sample_info = test_sample_info
    )

    # Verify all NA genes are removed, leaving empty object
    all_na_genes_removed <- remove_na_genes(snp_data_all_na_genes)
    expect_equal(nrow(all_na_genes_removed), 0)
})

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
test_barcode_info <- data.frame(
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
        barcode_info = test_barcode_info
    )
}

# ==============================================================================

test_that("remove_doublets returns SNPData object", {
    snp_data <- create_test_snpdata()

    filtered_data <- remove_doublets(snp_data)

    # Verify filtered data is still SNPData object
    expect_s4_class(filtered_data, "SNPData")
})

test_that("remove_doublets removes doublet cells", {
    snp_data <- create_test_snpdata()

    filtered_data <- remove_doublets(snp_data)

    # Check that doublet cell was removed (only 1 cell should remain)
    expect_equal(ncol(filtered_data), 1)
})

test_that("remove_doublets keeps only non-doublet cells", {
    snp_data <- create_test_snpdata()

    filtered_data <- remove_doublets(snp_data)

    # Verify remaining cell is the non-doublet donor
    remaining_barcode_info <- get_barcode_info(filtered_data)
    expect_equal(remaining_barcode_info$donor, "donor_1")
})

test_that("remove_doublets correctly subsets count matrices", {
    snp_data <- create_test_snpdata()

    filtered_data <- remove_doublets(snp_data)

    # Check that matrices were properly subsetted to match remaining cells
    expect_equal(dim(alt_count(filtered_data)), c(3, 1))
    # Verify ref_count matrix dimensions match alt_count
    expect_equal(dim(ref_count(filtered_data)), c(3, 1))
})

test_that("remove_doublets preserves SNP info", {
    snp_data <- create_test_snpdata()

    filtered_data <- remove_doublets(snp_data)

    # Verify SNP info unchanged (row subsetting not affected)
    expect_equal(nrow(get_snp_info(filtered_data)), 3)
})

test_that("remove_doublets handles missing donor column", {
    barcode_info_no_donor <- data.frame(
        cell_id = c("cell_1", "cell_2"),
        clonotype = c("clonotype_1", "clonotype_2"),
        stringsAsFactors = FALSE
    )

    snp_data_no_donor <- SNPData(
        alt_count = test_alt_count,
        ref_count = test_ref_count,
        snp_info = test_snp_info,
        barcode_info = barcode_info_no_donor
    )

    # Verify warning is issued when donor column is missing
    expect_warning(
        result <- remove_doublets(snp_data_no_donor),
        "No 'donor' column found in barcode_info, returning original object"
    )

    # Check that original number of columns is preserved
    expect_equal(ncol(result), ncol(snp_data_no_donor))
    # Verify original number of rows is preserved
    expect_equal(nrow(result), nrow(snp_data_no_donor))
})

test_that("remove_doublets removes NA donors when drop_na is TRUE", {
    barcode_info_with_na <- data.frame(
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
        barcode_info = barcode_info_with_na
    )

    filtered_drop_na <- remove_doublets(snp_data_with_na, drop_na = TRUE)

    # Verify NA donor cells are removed when drop_na = TRUE
    expect_equal(ncol(filtered_drop_na), 2)
})

test_that("remove_doublets keeps NA donors when drop_na is FALSE", {
    barcode_info_with_na <- data.frame(
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
        barcode_info = barcode_info_with_na
    )

    filtered_keep_na <- remove_doublets(snp_data_with_na, drop_na = FALSE)

    # Verify NA donor cells are kept when drop_na = FALSE
    expect_equal(ncol(filtered_keep_na), 3)
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

test_that("remove_na_genes returns SNPData object", {
    snp_data <- create_test_snpdata()

    filtered_data <- remove_na_genes(snp_data)

    # Verify filtered data is still SNPData object
    expect_s4_class(filtered_data, "SNPData")
})

test_that("remove_na_genes removes SNPs with NA gene names", {
    snp_data <- create_test_snpdata()

    filtered_data <- remove_na_genes(snp_data)

    # Check that SNP with NA gene_name was removed (2 SNPs should remain)
    expect_equal(nrow(filtered_data), 2)
})

test_that("remove_na_genes keeps only SNPs with valid gene names", {
    snp_data <- create_test_snpdata()

    filtered_data <- remove_na_genes(snp_data)

    remaining_snp_info <- get_snp_info(filtered_data)
    # Verify remaining SNPs don't have NA gene names
    expect_false(any(is.na(remaining_snp_info$gene_name)))
    # Check that correct genes are retained
    expect_equal(remaining_snp_info$gene_name, c("GENE1", "GENE3"))
})

test_that("remove_na_genes correctly subsets count matrices", {
    snp_data <- create_test_snpdata()

    filtered_data <- remove_na_genes(snp_data)

    # Check that matrices were properly subsetted to match remaining SNPs
    expect_equal(dim(alt_count(filtered_data)), c(2, 2))
    # Verify ref_count matrix dimensions match alt_count
    expect_equal(dim(ref_count(filtered_data)), c(2, 2))
})

test_that("remove_na_genes preserves barcode info", {
    snp_data <- create_test_snpdata()

    filtered_data <- remove_na_genes(snp_data)

    # Verify barcode info structure is unchanged (column subsetting not affected)
    expect_equal(ncol(get_barcode_info(filtered_data)), ncol(get_barcode_info(snp_data)))
})

test_that("remove_na_genes handles missing gene column", {
    snp_info_no_gene <- data.frame(
        snp_id = c("snp_1", "snp_2", "snp_3"),
        pos = c(100, 200, 300),
        stringsAsFactors = FALSE
    )

    snp_data_no_gene <- SNPData(
        alt_count = test_alt_count,
        ref_count = test_ref_count,
        snp_info = snp_info_no_gene,
        barcode_info = test_barcode_info
    )

    # Verify warning is issued when gene column is missing
    expect_warning(
        result <- remove_na_genes(snp_data_no_gene),
        "No 'gene_name' column found in snp_info, returning original object"
    )

    # Check that original number of rows is preserved
    expect_equal(nrow(result), nrow(snp_data_no_gene))
    # Verify original number of columns is preserved
    expect_equal(ncol(result), ncol(snp_data_no_gene))
})

test_that("remove_na_genes handles custom gene column", {
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
        barcode_info = test_barcode_info
    )

    filtered_data <- remove_na_genes(snp_data_custom, gene_col = "custom_gene")

    # Verify filtering worked with custom column
    expect_equal(nrow(filtered_data), 2)
    remaining_snp_info <- get_snp_info(filtered_data)
    # Check that NA genes are removed from custom column
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

test_that("remove_na_clonotypes returns SNPData object", {
    snp_data <- create_test_snpdata()

    filtered_data <- remove_na_clonotypes(snp_data)

    # Verify filtered data is still SNPData object
    expect_s4_class(filtered_data, "SNPData")
})

test_that("remove_na_clonotypes removes cells with NA clonotypes", {
    snp_data <- create_test_snpdata()

    filtered_data <- remove_na_clonotypes(snp_data)

    # Check that cell with NA clonotype was removed (1 cell should remain)
    expect_equal(ncol(filtered_data), 1)
})

test_that("remove_na_clonotypes keeps only cells with valid clonotypes", {
    snp_data <- create_test_snpdata()

    filtered_data <- remove_na_clonotypes(snp_data)

    remaining_barcode_info <- get_barcode_info(filtered_data)
    # Verify remaining cells don't have NA clonotypes
    expect_false(any(is.na(remaining_barcode_info$clonotype)))
    # Check that correct clonotype is retained
    expect_equal(remaining_barcode_info$clonotype, "clonotype_1")
})

test_that("remove_na_clonotypes correctly subsets count matrices", {
    snp_data <- create_test_snpdata()

    filtered_data <- remove_na_clonotypes(snp_data)

    # Check that matrices were properly subsetted to match remaining cells
    expect_equal(dim(alt_count(filtered_data)), c(3, 1))
    # Verify ref_count matrix dimensions match alt_count
    expect_equal(dim(ref_count(filtered_data)), c(3, 1))
})

test_that("remove_na_clonotypes preserves SNP info", {
    snp_data <- create_test_snpdata()

    filtered_data <- remove_na_clonotypes(snp_data)

    # Verify SNP info unchanged (row subsetting not affected)
    expect_equal(nrow(get_snp_info(filtered_data)), 3)
})

test_that("remove_na_clonotypes handles missing clonotype column", {
    barcode_info_no_clonotype <- data.frame(
        cell_id = c("cell_1", "cell_2"),
        donor = c("donor_1", "donor_2"),
        stringsAsFactors = FALSE
    )

    snp_data_no_clonotype <- SNPData(
        alt_count = test_alt_count,
        ref_count = test_ref_count,
        snp_info = test_snp_info,
        barcode_info = barcode_info_no_clonotype
    )

    # Verify warning is issued when clonotype column is missing
    expect_warning(
        result <- remove_na_clonotypes(snp_data_no_clonotype),
        "No 'clonotype' column found.*add_barcode_metadata"
    )

    # Check that original number of columns is preserved
    expect_equal(ncol(result), ncol(snp_data_no_clonotype))
    # Verify original number of rows is preserved
    expect_equal(nrow(result), nrow(snp_data_no_clonotype))
})

test_that("remove_na_clonotypes handles custom clonotype column", {
    barcode_info_custom <- data.frame(
        cell_id = c("cell_1", "cell_2"),
        donor = c("donor_1", "donor_2"),
        custom_clonotype = c("clonotype_1", NA),
        stringsAsFactors = FALSE
    )

    snp_data_custom <- SNPData(
        alt_count = test_alt_count,
        ref_count = test_ref_count,
        snp_info = test_snp_info,
        barcode_info = barcode_info_custom
    )

    filtered_data <- remove_na_clonotypes(snp_data_custom, clonotype_col = "custom_clonotype")

    # Verify filtering worked with custom column
    expect_equal(ncol(filtered_data), 1)
    remaining_barcode_info <- get_barcode_info(filtered_data)
    # Check that NA clonotypes are removed
    expect_false(any(is.na(remaining_barcode_info$custom_clonotype)))
})

test_that("remove_na_clonotypes validates input", {
    # Test with non-SNPData object
    # Verify error when input is not SNPData object
    expect_error(
        remove_na_clonotypes("not_snpdata"),
        "Input must be a SNPData object"
    )
})

test_that("remove_doublets handles all doublets edge case", {
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
        barcode_info = all_doublets_info
    )

    all_doublets_removed <- remove_doublets(snp_data_all_doublets)

    # Verify all doublets are removed, leaving empty object
    expect_equal(ncol(all_doublets_removed), 0)
})

test_that("remove_na_genes handles all NA genes edge case", {
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
        barcode_info = test_barcode_info
    )

    all_na_genes_removed <- remove_na_genes(snp_data_all_na_genes)

    # Verify all NA genes are removed, leaving empty object
    expect_equal(nrow(all_na_genes_removed), 0)
})

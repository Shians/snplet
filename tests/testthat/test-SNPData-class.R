# ==============================================================================
# Test Suite: SNPData S4 Class
# Description: Tests for SNPData class constructor, accessors, and core methods
# ==============================================================================

library(testthat)
library(Matrix)

# ------------------------------------------------------------------------------
# Test Data Setup
# ------------------------------------------------------------------------------

# Create test matrices
test_alt_count <- Matrix::Matrix(matrix(c(1, 2, 3, 4), nrow = 2, ncol = 2))
dimnames(test_alt_count) <- list(c("snp_1", "snp_2"), c("cell_1", "cell_2"))
test_ref_count <- Matrix::Matrix(matrix(c(5, 6, 7, 8), nrow = 2, ncol = 2))
dimnames(test_ref_count) <- list(c("snp_1", "snp_2"), c("cell_1", "cell_2"))

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
    stringsAsFactors = FALSE
)

# ==============================================================================

test_that("SNPData() creates valid S4 object with correct class", {
    snp_data <- SNPData(
        ref_count = test_ref_count,
        alt_count = test_alt_count,
        snp_info = test_snp_info,
        barcode_info = test_barcode_info
    )

    # Verify SNPData object is created successfully
    expect_s4_class(snp_data, "SNPData")
})

test_that("SNPData() sets correct dimensions", {
    snp_data <- SNPData(
        ref_count = test_ref_count,
        alt_count = test_alt_count,
        snp_info = test_snp_info,
        barcode_info = test_barcode_info
    )

    # Check that dim() method returns correct dimensions
    expect_equal(dim(snp_data), c(2, 2))
    # Verify nrow() method works correctly
    expect_equal(nrow(snp_data), 2)
    # Verify ncol() method works correctly
    expect_equal(ncol(snp_data), 2)
})

test_that("SNPData() sets row and column names from metadata", {
    snp_data <- SNPData(
        ref_count = test_ref_count,
        alt_count = test_alt_count,
        snp_info = test_snp_info,
        barcode_info = test_barcode_info
    )

    # Check that SNP IDs are used as row names
    expect_equal(rownames(snp_data), c("snp_1", "snp_2"))
    # Check that cell IDs are used as column names
    expect_equal(colnames(snp_data), c("cell_1", "cell_2"))
})

test_that("ref_count() returns reference count matrix with correct values and names", {
    snp_data <- SNPData(
        ref_count = test_ref_count,
        alt_count = test_alt_count,
        snp_info = test_snp_info,
        barcode_info = test_barcode_info
    )

    ref_count_matrix <- ref_count(snp_data)

    # Verify ref_count accessor returns correct matrix values
    expect_equal(as.matrix(ref_count_matrix), as.matrix(test_ref_count))
    # Check that ref_count matrix has correct row names
    expect_equal(rownames(ref_count_matrix), c("snp_1", "snp_2"))
    # Check that ref_count matrix has correct column names
    expect_equal(colnames(ref_count_matrix), c("cell_1", "cell_2"))
})

test_that("alt_count() returns alternate count matrix with correct values and names", {
    snp_data <- SNPData(
        ref_count = test_ref_count,
        alt_count = test_alt_count,
        snp_info = test_snp_info,
        barcode_info = test_barcode_info
    )

    alt_count_matrix <- alt_count(snp_data)

    # Verify alt_count accessor returns correct matrix values
    expect_equal(as.matrix(alt_count_matrix), as.matrix(test_alt_count))
    # Check that alt_count matrix has correct row names
    expect_equal(rownames(alt_count_matrix), c("snp_1", "snp_2"))
    # Check that alt_count matrix has correct column names
    expect_equal(colnames(alt_count_matrix), c("cell_1", "cell_2"))
})

test_that("get_snp_info() returns SNP metadata with computed coverage metrics", {
    snp_data <- SNPData(
        ref_count = test_ref_count,
        alt_count = test_alt_count,
        snp_info = test_snp_info,
        barcode_info = test_barcode_info
    )

    snp_info <- get_snp_info(snp_data)

    # Verify SNP IDs are preserved from input
    expect_equal(snp_info$snp_id, test_snp_info$snp_id)
    # Verify SNP positions are preserved from input
    expect_equal(snp_info$pos, test_snp_info$pos)
    # Check that coverage column is automatically added
    expect_true("coverage" %in% colnames(snp_info))
    # Check that non_zero_samples column is automatically added
    expect_true("non_zero_samples" %in% colnames(snp_info))
    # Verify coverage calculation: rowSums of alt_count + ref_count
    expect_equal(as.numeric(snp_info$coverage), c(16, 20))
    # Verify non_zero_samples count: all samples have counts
    expect_equal(as.numeric(snp_info$non_zero_samples), c(2, 2))
})

test_that("get_barcode_info() returns barcode metadata with computed library size metrics", {
    snp_data <- SNPData(
        ref_count = test_ref_count,
        alt_count = test_alt_count,
        snp_info = test_snp_info,
        barcode_info = test_barcode_info
    )

    barcode_info <- get_barcode_info(snp_data)

    # Verify cell IDs are preserved from input
    expect_equal(barcode_info$cell_id, test_barcode_info$cell_id)
    # Verify donor information is preserved from input
    expect_equal(barcode_info$donor, test_barcode_info$donor)
    # Verify clonotype information is preserved from input
    expect_equal(barcode_info$clonotype, test_barcode_info$clonotype)
    # Check that library_size column is automatically added
    expect_true("library_size" %in% colnames(barcode_info))
    # Check that non_zero_snps column is automatically added
    expect_true("non_zero_snps" %in% colnames(barcode_info))
    # Verify library_size calculation: colSums of alt_count + ref_count
    expect_equal(as.numeric(barcode_info$library_size), c(14, 22))
    # Verify non_zero_snps count: all SNPs have counts
    expect_equal(as.numeric(barcode_info$non_zero_snps), c(2, 2))
})

test_that("[() subsetting by numeric index returns SNPData object with correct dimensions", {
    snp_data <- SNPData(
        alt_count = test_alt_count,
        ref_count = test_ref_count,
        snp_info = test_snp_info,
        barcode_info = test_barcode_info
    )

    subset_data <- snp_data[1, 1]

    # Verify subsetting by index returns SNPData object
    expect_s4_class(subset_data, "SNPData")
    # Check that single element subset has correct dimensions
    expect_equal(dim(subset_data), c(1, 1))
})

test_that("[() subsetting by name returns SNPData object with correct dimensions", {
    snp_data <- SNPData(
        alt_count = test_alt_count,
        ref_count = test_ref_count,
        snp_info = test_snp_info,
        barcode_info = test_barcode_info
    )

    subset_data <- snp_data["snp_1", "cell_1"]

    # Verify subsetting by name returns SNPData object
    expect_s4_class(subset_data, "SNPData")
    # Check that named subset has correct dimensions
    expect_equal(dim(subset_data), c(1, 1))
})

test_that("[() ignores drop parameter and always returns SNPData object", {
    snp_data <- SNPData(
        alt_count = test_alt_count,
        ref_count = test_ref_count,
        snp_info = test_snp_info,
        barcode_info = test_barcode_info
    )

    subset_data <- snp_data[1, 1, drop = TRUE]

    # Verify drop=TRUE is ignored and still returns SNPData object
    expect_s4_class(subset_data, "SNPData")
    # Check that drop=TRUE doesn't affect dimensions
    expect_equal(dim(subset_data), c(1, 1))
})

test_that("[() subsetting to single row returns SNPData object", {
    snp_data <- SNPData(
        alt_count = test_alt_count,
        ref_count = test_ref_count,
        snp_info = test_snp_info,
        barcode_info = test_barcode_info
    )

    subset_data <- snp_data[1, ]

    # Verify single row subset returns SNPData object
    expect_s4_class(subset_data, "SNPData")
    # Check that single row subset has correct dimensions
    expect_equal(dim(subset_data), c(1, 2))
})

test_that("[() subsetting to single column returns SNPData object", {
    snp_data <- SNPData(
        alt_count = test_alt_count,
        ref_count = test_ref_count,
        snp_info = test_snp_info,
        barcode_info = test_barcode_info
    )

    subset_data <- snp_data[, 1]

    # Verify single column subset returns SNPData object
    expect_s4_class(subset_data, "SNPData")
    # Check that single column subset has correct dimensions
    expect_equal(dim(subset_data), c(2, 1))
})

test_that("coverage() returns sum of alt_count and ref_count matrices", {
    snp_data <- SNPData(
        alt_count = test_alt_count,
        ref_count = test_ref_count,
        snp_info = test_snp_info,
        barcode_info = test_barcode_info
    )
    expected_coverage <- Matrix::Matrix(
        matrix(c(6, 8, 10, 12), nrow = 2, ncol = 2),
        dimnames = list(c("snp_1", "snp_2"), c("cell_1", "cell_2"))
    )

    result <- coverage(snp_data)

    # Verify coverage calculation: alt_count + ref_count
    expect_equal(result, expected_coverage)
})
test_that("SNPData() auto-generates snp_id column when missing from snp_info", {
    test_snp_info_no_id <- data.frame(pos = c(100, 200))
    test_barcode_info_no_id <- data.frame(
        donor = c("donor_1", "donor_1"),
        clonotype = c("clonotype_1", "clonotype_2")
    )

    snp_data <- SNPData(
        alt_count = test_alt_count,
        ref_count = test_ref_count,
        snp_info = test_snp_info_no_id,
        barcode_info = test_barcode_info_no_id
    )

    # Verify snp_id column was auto-generated when missing
    expect_true("snp_id" %in% colnames(get_snp_info(snp_data)))
    # Check that auto-generated SNP IDs follow expected pattern
    expect_equal(get_snp_info(snp_data)$snp_id, c("snp_1", "snp_2"))
})

test_that("SNPData() auto-generates cell_id column when missing from barcode_info", {
    test_snp_info_no_id <- data.frame(pos = c(100, 200))
    test_barcode_info_no_id <- data.frame(
        donor = c("donor_1", "donor_1"),
        clonotype = c("clonotype_1", "clonotype_2")
    )

    snp_data <- SNPData(
        alt_count = test_alt_count,
        ref_count = test_ref_count,
        snp_info = test_snp_info_no_id,
        barcode_info = test_barcode_info_no_id
    )

    # Verify cell_id column was auto-generated when missing
    expect_true("cell_id" %in% colnames(get_barcode_info(snp_data)))
    # Check that auto-generated cell IDs follow expected pattern
    expect_equal(get_barcode_info(snp_data)$cell_id, c("cell_1", "cell_2"))
})

test_that("SNPData() throws error when alt_count and ref_count have mismatched row dimensions", {
    wrong_dim_alt_count <- Matrix::Matrix(matrix(c(1, 2, 3), nrow = 3, ncol = 1))

    # Verify error when alt_count and ref_count have different row counts
    expect_error(
        SNPData(
            ref_count = test_ref_count,
            alt_count = wrong_dim_alt_count,
            snp_info = test_snp_info,
            barcode_info = test_barcode_info
        ),
        "nrow\\(alt_count\\) == nrow\\(ref_count\\)"
    )
})

test_that("SNPData() throws error when snp_info rows don't match matrix rows", {
    wrong_dim_snp_info <- data.frame(snp_id = c("snp_1", "snp_2", "snp_3"))

    # Verify error when snp_info rows don't match matrix rows
    expect_error(
        SNPData(
            ref_count = test_ref_count,
            alt_count = test_alt_count,
            snp_info = wrong_dim_snp_info,
            barcode_info = test_barcode_info
        ),
        "nrow\\(ref_count\\) == nrow\\(snp_info\\)"
    )
})

test_that("SNPData() throws error when barcode_info rows don't match matrix columns", {
    wrong_dim_barcode_info <- data.frame(cell_id = c("cell_1", "cell_2", "cell_3"))

    # Verify error when barcode_info rows don't match matrix columns
    expect_error(
        SNPData(
            alt_count = test_alt_count,
            ref_count = test_ref_count,
            snp_info = test_snp_info,
            barcode_info = wrong_dim_barcode_info
        ),
        "ncol\\(alt_count\\) == nrow\\(barcode_info\\)"
    )
})

test_that("get_sample_info() returns same result as get_barcode_info()", {
    snp_data <- SNPData(
        ref_count = test_ref_count,
        alt_count = test_alt_count,
        snp_info = test_snp_info,
        barcode_info = test_barcode_info
    )

    # Verify get_sample_info is an alias for get_barcode_info
    expect_equal(get_sample_info(snp_data), get_barcode_info(snp_data))
})

test_that("barcode_info<- replaces barcode_info with valid data", {
    snp_data <- SNPData(
        ref_count = test_ref_count,
        alt_count = test_alt_count,
        snp_info = test_snp_info,
        barcode_info = test_barcode_info
    )

    new_barcode_info <- data.frame(
        cell_id = c("cell_1", "cell_2"),
        donor = c("donor_2", "donor_2"),
        clonotype = c("clonotype_3", "clonotype_4"),
        stringsAsFactors = FALSE
    )

    barcode_info(snp_data) <- new_barcode_info

    # Verify barcode_info was updated
    expect_equal(get_barcode_info(snp_data)$donor, c("donor_2", "donor_2"))
    # Verify clonotype was updated
    expect_equal(get_barcode_info(snp_data)$clonotype, c("clonotype_3", "clonotype_4"))
})

test_that("barcode_info<- throws error when dimensions don't match", {
    snp_data <- SNPData(
        ref_count = test_ref_count,
        alt_count = test_alt_count,
        snp_info = test_snp_info,
        barcode_info = test_barcode_info
    )

    wrong_dim_barcode_info <- data.frame(
        cell_id = c("cell_1", "cell_2", "cell_3"),
        donor = c("donor_1", "donor_1", "donor_1")
    )

    # Verify error when number of rows doesn't match
    expect_error(
        barcode_info(snp_data) <- wrong_dim_barcode_info,
        "Number of rows in barcode_info must match number of columns in count matrices"
    )
})

test_that("barcode_info<- throws error when cell_id column is missing", {
    snp_data <- SNPData(
        ref_count = test_ref_count,
        alt_count = test_alt_count,
        snp_info = test_snp_info,
        barcode_info = test_barcode_info
    )

    no_cell_id <- data.frame(
        donor = c("donor_1", "donor_1"),
        clonotype = c("clonotype_1", "clonotype_2")
    )

    # Verify error when cell_id column is missing
    expect_error(
        barcode_info(snp_data) <- no_cell_id,
        "barcode_info must contain a 'cell_id' column"
    )
})

test_that("barcode_info<- throws error when cell_id doesn't match matrix column names", {
    snp_data <- SNPData(
        ref_count = test_ref_count,
        alt_count = test_alt_count,
        snp_info = test_snp_info,
        barcode_info = test_barcode_info
    )

    mismatched_cell_id <- data.frame(
        cell_id = c("cell_3", "cell_4"),
        donor = c("donor_1", "donor_1")
    )

    # Verify error when cell_id doesn't match column names
    expect_error(
        barcode_info(snp_data) <- mismatched_cell_id,
        "barcode_info\\$cell_id must match column names of count matrices"
    )
})

test_that("snp_info<- replaces snp_info with valid data", {
    snp_data <- SNPData(
        ref_count = test_ref_count,
        alt_count = test_alt_count,
        snp_info = test_snp_info,
        barcode_info = test_barcode_info
    )

    new_snp_info <- data.frame(
        snp_id = c("snp_1", "snp_2"),
        pos = c(150, 250),
        chr = c("chr1", "chr2"),
        stringsAsFactors = FALSE
    )

    snp_info(snp_data) <- new_snp_info

    # Verify snp_info was updated
    expect_equal(get_snp_info(snp_data)$pos, c(150, 250))
    # Verify chr column was added
    expect_equal(get_snp_info(snp_data)$chr, c("chr1", "chr2"))
})

test_that("snp_info<- throws error when dimensions don't match", {
    snp_data <- SNPData(
        ref_count = test_ref_count,
        alt_count = test_alt_count,
        snp_info = test_snp_info,
        barcode_info = test_barcode_info
    )

    wrong_dim_snp_info <- data.frame(
        snp_id = c("snp_1", "snp_2", "snp_3"),
        pos = c(100, 200, 300)
    )

    # Verify error when number of rows doesn't match
    expect_error(
        snp_info(snp_data) <- wrong_dim_snp_info,
        "Number of rows in snp_info must match number of rows in count matrices"
    )
})

test_that("snp_info<- throws error when snp_id column is missing", {
    snp_data <- SNPData(
        ref_count = test_ref_count,
        alt_count = test_alt_count,
        snp_info = test_snp_info,
        barcode_info = test_barcode_info
    )

    no_snp_id <- data.frame(
        pos = c(100, 200),
        chr = c("chr1", "chr2")
    )

    # Verify error when snp_id column is missing
    expect_error(
        snp_info(snp_data) <- no_snp_id,
        "snp_info must contain a 'snp_id' column"
    )
})

test_that("snp_info<- throws error when snp_id doesn't match matrix row names", {
    snp_data <- SNPData(
        ref_count = test_ref_count,
        alt_count = test_alt_count,
        snp_info = test_snp_info,
        barcode_info = test_barcode_info
    )

    mismatched_snp_id <- data.frame(
        snp_id = c("snp_3", "snp_4"),
        pos = c(100, 200)
    )

    # Verify error when snp_id doesn't match row names
    expect_error(
        snp_info(snp_data) <- mismatched_snp_id,
        "snp_info\\$snp_id must match row names of count matrices"
    )
})

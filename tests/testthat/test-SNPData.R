library(testthat)
library(Matrix)

# Helper function to compare objects without names
expect_equal_unnamed <- function(object, expected, ...) {
    expect_equal(unname(object), unname(expected), ...)
}

# Create test data
test_alt_count <- Matrix::Matrix(matrix(c(1, 2, 3, 4), nrow = 2, ncol = 2))
test_ref_count <- Matrix::Matrix(matrix(c(5, 6, 7, 8), nrow = 2, ncol = 2))
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

test_that("SNPData constructor works correctly", {
    # Test basic construction
    snp_data <- SNPData(
        ref_count = test_ref_count,
        alt_count = test_alt_count,
        snp_info = test_snp_info,
        sample_info = test_sample_info
    )
    expect_s4_class(snp_data, "SNPData")

    # Test dimension compatibility
    expect_equal(dim(snp_data), c(2, 2))
    expect_equal(nrow(snp_data), 2)
    expect_equal(ncol(snp_data), 2)

    # Test rownames and colnames
    expect_equal(rownames(snp_data), c("snp_1", "snp_2"))
    expect_equal(colnames(snp_data), c("cell_1", "cell_2"))
})

test_that("SNPData accessors work correctly", {
    snp_data <- SNPData(
        ref_count = test_ref_count,
        alt_count = test_alt_count,
        snp_info = test_snp_info,
        sample_info = test_sample_info
    )

    # Test ref_count accessor
    ref_count_matrix <- ref_count(snp_data)
    expect_equal_unnamed(ref_count_matrix, test_ref_count)
    expect_equal(rownames(ref_count_matrix), c("snp_1", "snp_2"))
    expect_equal(colnames(ref_count_matrix), c("cell_1", "cell_2"))

    # Test alt_count accessor
    alt_count_matrix <- alt_count(snp_data)
    expect_equal_unnamed(alt_count_matrix, test_alt_count)
    expect_equal(rownames(alt_count_matrix), c("snp_1", "snp_2"))
    expect_equal(colnames(alt_count_matrix), c("cell_1", "cell_2"))

    # Test get_snp_info accessor
    snp_info <- get_snp_info(snp_data)
    expect_equal(snp_info$snp_id, test_snp_info$snp_id)
    expect_equal(snp_info$pos, test_snp_info$pos)
    expect_true("coverage" %in% colnames(snp_info))
    expect_true("non_zero_samples" %in% colnames(snp_info))
    expect_equal(snp_info$coverage, c(16, 20))  # rowSums of alt_count + ref_count
    expect_equal(snp_info$non_zero_samples, c(2, 2))  # all samples have counts

    # Test get_sample_info accessor
    sample_info <- get_sample_info(snp_data)
    expect_equal(sample_info$cell_id, test_sample_info$cell_id)
    expect_equal(sample_info$donor_id, test_sample_info$donor_id)
    expect_equal(sample_info$clonotype_id, test_sample_info$clonotype_id)
    expect_true("library_size" %in% colnames(sample_info))
    expect_true("non_zero_snps" %in% colnames(sample_info))
    expect_equal(sample_info$library_size, c(14, 22))  # colSums of alt_count + ref_count
    expect_equal(sample_info$non_zero_snps, c(2, 2))  # all SNPs have counts
})

test_that("SNPData subsetting works correctly", {
    snp_data <- SNPData(
        alt_count = test_alt_count,
        ref_count = test_ref_count,
        snp_info = test_snp_info,
        sample_info = test_sample_info
    )

    # Test subsetting by index
    subset_data <- snp_data[1, 1]
    expect_s4_class(subset_data, "SNPData")
    expect_equal(dim(subset_data), c(1, 1))

    # Test subsetting by name
    subset_data <- snp_data["snp_1", "cell_1"]
    expect_s4_class(subset_data, "SNPData")
    expect_equal(dim(subset_data), c(1, 1))

    # Test that drop parameter is ignored
    subset_data <- snp_data[1, 1, drop = TRUE]
    expect_s4_class(subset_data, "SNPData")
    expect_equal(dim(subset_data), c(1, 1))

    # Test subsetting to single row
    subset_data <- snp_data[1, ]
    expect_s4_class(subset_data, "SNPData")
    expect_equal(dim(subset_data), c(1, 2))

    # Test subsetting to single column
    subset_data <- snp_data[, 1]
    expect_s4_class(subset_data, "SNPData")
    expect_equal(dim(subset_data), c(2, 1))
})

test_that("SNPData coverage calculations work correctly", {
    snp_data <- SNPData(
        alt_count = test_alt_count,
        ref_count = test_ref_count,
        snp_info = test_snp_info,
        sample_info = test_sample_info
    )

    # Test coverage method
    expected_coverage <- Matrix::Matrix(matrix(c(6, 8, 10, 12), nrow = 2, ncol = 2))
    expect_equal_unnamed(coverage(snp_data), expected_coverage)

    # Test ref_fraction method
    expected_ref_fraction <- Matrix::Matrix(matrix(c(5/6, 6/8, 7/10, 8/12), nrow = 2, ncol = 2))
    expect_equal_unnamed(ref_fraction(snp_data), expected_ref_fraction)

    # Test major_allele_frac method
    expected_major_allele_frac <- Matrix::Matrix(matrix(c(0.5 + abs(5/6 - 0.5), 0.5 + abs(6/8 - 0.5),
                                                    0.5 + abs(7/10 - 0.5), 0.5 + abs(8/12 - 0.5)),
                                           nrow = 2, ncol = 2))
    expect_equal_unnamed(major_allele_frac(snp_data), expected_major_allele_frac)
})

test_that("SNPData count data frame methods work correctly", {
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
    expect_true(all(c("snp_id", "cell_id", "ref_count", "alt_count", "total_count", "ref_ratio", "maf") %in% colnames(cell_df)))

    # Test donor_count_df method
    donor_df <- donor_count_df(snp_data)
    expect_s3_class(donor_df, "tbl_df")
    expect_equal(nrow(donor_df), 2)  # 2 SNPs x 1 donor
    expect_true(all(c("snp_id", "donor", "ref_count", "alt_count", "total_count", "ref_ratio", "maf") %in% colnames(donor_df)))

    # Test clonotype_count_df method
    clonotype_df <- clonotype_count_df(snp_data)
    expect_s3_class(clonotype_df, "tbl_df")
    expect_equal(nrow(clonotype_df), 4)  # 2 SNPs x 2 clonotypes
    expect_true(all(c("snp_id", "clonotype", "ref_count", "alt_count", "total_count", "ref_ratio", "maf", "donor") %in% colnames(clonotype_df)))
})

test_that("SNPData handles missing IDs correctly", {
    # Create data without snp_id and cell_id
    test_snp_info_no_id <- data.frame(pos = c(100, 200))
    test_sample_info_no_id <- data.frame(
        donor_id = c("donor_1", "donor_1"),
        clonotype_id = c("clonotype_1", "clonotype_2")
    )

    snp_data <- SNPData(
        alt_count = test_alt_count,
        ref_count = test_ref_count,
        snp_info = test_snp_info_no_id,
        sample_info = test_sample_info_no_id
    )

    # Check that IDs were automatically assigned
    expect_true("snp_id" %in% colnames(get_snp_info(snp_data)))
    expect_true("cell_id" %in% colnames(get_sample_info(snp_data)))
    expect_equal(get_snp_info(snp_data)$snp_id, c("snp_1", "snp_2"))
    expect_equal(get_sample_info(snp_data)$cell_id, c("cell_1", "cell_2"))
})

test_that("SNPData validates input dimensions", {
    # Test mismatched dimensions
    wrong_dim_alt_count <- Matrix::Matrix(matrix(c(1, 2, 3), nrow = 3, ncol = 1))
    expect_error(SNPData(
        ref_count = test_ref_count,
        alt_count = wrong_dim_alt_count,
        snp_info = test_snp_info,
        sample_info = test_sample_info
    ),
    "nrow\\(alt_count\\) == nrow\\(ref_count\\)")

    wrong_dim_snp_info <- data.frame(snp_id = c("snp_1", "snp_2", "snp_3"))
    expect_error(SNPData(
        ref_count = test_ref_count,
        alt_count = test_alt_count,
        snp_info = wrong_dim_snp_info,
        sample_info = test_sample_info
    ),
    "nrow\\(ref_count\\) == nrow\\(snp_info\\)")

    wrong_dim_sample_info <- data.frame(cell_id = c("cell_1", "cell_2", "cell_3"))
    expect_error(SNPData(
        alt_count = test_alt_count,
        ref_count = test_ref_count,
        snp_info = test_snp_info,
        sample_info = wrong_dim_sample_info
    ),
    "ncol\\(alt_count\\) == nrow\\(sample_info\\)")
})

test_that("filter_snps and filter_barcodes check for missing columns", {
    snp_data <- SNPData(
        alt_count = test_alt_count,
        ref_count = test_ref_count,
        snp_info = test_snp_info,
        sample_info = test_sample_info
    )

    # filter_snps: should error if column does not exist
    expect_error(
        filter_snps(snp_data, not_a_column > 0),
        "The following columns are not present in snp_info or parent environment: not_a_column"
    )

    # filter_barcodes: should error if column does not exist
    expect_error(
        filter_barcodes(snp_data, not_a_column > 0),
        "The following columns are not present in sample_info or parent environment: not_a_column"
    )

    # filter_snps: should not error if column exists
    expect_s4_class(filter_snps(snp_data, pos > 100), "SNPData")

    # filter_barcodes: should not error if column exists
    expect_s4_class(filter_barcodes(snp_data, donor == 'donor_1'), "SNPData")
})

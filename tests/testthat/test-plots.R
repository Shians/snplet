# ==============================================================================
# Test Suite: Distribution Plots
# Description: Tests for library size, SNP coverage, and MAF distribution plots
# ==============================================================================

library(testthat)
library(Matrix)
library(ggplot2)

# ------------------------------------------------------------------------------
# Test Data Setup
# ------------------------------------------------------------------------------

create_test_snpdata_for_plots <- function() {
    alt_count <- Matrix::Matrix(matrix(c(1, 2, 3, 4), nrow = 2, ncol = 2))
    ref_count <- Matrix::Matrix(matrix(c(5, 6, 7, 8), nrow = 2, ncol = 2))
    snp_info <- data.frame(snp_id = c("snp_1", "snp_2"), pos = c(100, 200))
    barcode_info <- data.frame(cell_id = c("cell_1", "cell_2"))

    SNPData(
        ref_count = ref_count,
        alt_count = alt_count,
        snp_info = snp_info,
        barcode_info = barcode_info
    )
}

# ==============================================================================

test_that("plot_lib_size_distribution() returns a ggplot object", {
    snp_data <- create_test_snpdata_for_plots()

    p <- plot_lib_size_distribution(snp_data)

    # Verify a ggplot object is returned
    expect_s3_class(p, "ggplot")
    # Verify the plot builds without error
    expect_no_error(ggplot_build(p))
})

test_that("plot_lib_size_distribution() maps x to library_size", {
    snp_data <- create_test_snpdata_for_plots()

    p <- plot_lib_size_distribution(snp_data)

    # Verify the x aesthetic is mapped to library_size
    expect_equal(rlang::as_label(p$mapping$x), "library_size")
})

test_that("plot_snp_cov_distribution() returns a ggplot object", {
    snp_data <- create_test_snpdata_for_plots()

    p <- plot_snp_cov_distribution(snp_data)

    # Verify a ggplot object is returned
    expect_s3_class(p, "ggplot")
    # Verify the plot builds without error
    expect_no_error(ggplot_build(p))
})

test_that("plot_snp_cov_distribution() maps x to coverage", {
    snp_data <- create_test_snpdata_for_plots()

    p <- plot_snp_cov_distribution(snp_data)

    # Verify the x aesthetic is mapped to coverage
    expect_equal(rlang::as_label(p$mapping$x), "coverage")
})

test_that("plot_maf_distribution() returns a ggplot object", {
    df <- data.frame(maf = c(0.1, 0.2, 0.3, 0.4, 0.5))

    p <- plot_maf_distribution(df)

    # Verify a ggplot object is returned
    expect_s3_class(p, "ggplot")
    # Verify the plot builds without error
    expect_no_error(ggplot_build(p))
})

test_that("plot_maf_distribution() maps x to maf", {
    df <- data.frame(maf = c(0.1, 0.2, 0.3))

    p <- plot_maf_distribution(df)

    # Verify the x aesthetic is mapped to maf
    expect_equal(rlang::as_label(p$mapping$x), "maf")
})

# ==============================================================================
# Test Suite: Genomic Track Plots
# Description: Tests for gene annotation, MAF, and MAF p-value track plots
# ==============================================================================

library(testthat)
library(ggplot2)

# ==============================================================================

test_that("plot_gene_anno_track() returns a ggplot object", {
    gene_anno <- data.frame(
        pos = c(1000, 2000, 3000),
        gene_name = c("GENE1", "GENE2", "GENE3")
    )

    p <- plot_gene_anno_track(gene_anno, x_range = c(500, 3500))

    # Verify a ggplot object is returned
    expect_s3_class(p, "ggplot")
    # Verify the plot builds without error
    expect_no_error(ggplot_build(p))
})

test_that("plot_maf_track() returns a ggplot object", {
    allele_count_df <- data.frame(
        pos = c(1000, 2000, 3000),
        maf = c(0.05, 0.2, 0.4),
        donor_id = c("donor1", "donor1", "donor2")
    )

    p <- plot_maf_track(allele_count_df, facet = donor_id, x_range = c(500, 3500))

    # Verify a ggplot object is returned
    expect_s3_class(p, "ggplot")
    # Verify the plot builds without error
    expect_no_error(ggplot_build(p))
})

test_that("plot_maf_track() maps y to maf", {
    allele_count_df <- data.frame(
        pos = c(1000, 2000),
        maf = c(0.05, 0.2),
        donor_id = c("donor1", "donor2")
    )

    p <- plot_maf_track(allele_count_df, facet = donor_id, x_range = c(500, 2500))

    # Verify the y aesthetic is mapped to maf
    expect_equal(rlang::as_label(p$mapping$y), "maf")
})

test_that("plot_maf_pval_track() returns a ggplot object", {
    allele_counts_df <- data.frame(
        pos = c(1000, 2000, 3000),
        adj_p_val = c(0.01, 0.2, 0.5),
        donor_id = c("donor1", "donor1", "donor2")
    )

    p <- plot_maf_pval_track(allele_counts_df, facet = donor_id, x_range = c(500, 3500))

    # Verify a ggplot object is returned
    expect_s3_class(p, "ggplot")
    # Verify the plot builds without error
    expect_no_error(ggplot_build(p))
})

test_that("plot_maf_pval_track() maps y to -log10(adj_p_val)", {
    allele_counts_df <- data.frame(
        pos = c(1000, 2000),
        adj_p_val = c(0.01, 0.5),
        donor_id = c("donor1", "donor2")
    )

    p <- plot_maf_pval_track(allele_counts_df, facet = donor_id, x_range = c(500, 2500))

    # Verify the y aesthetic reflects the -log10(adj_p_val) transform
    expect_equal(rlang::as_label(p$mapping$y), "-log10(adj_p_val)")
})

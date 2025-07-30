# ==============================================================================
# Test Suite: Data Import Functions
# Description: Tests for cellSNP data import, export, and utility functions
# ==============================================================================

library(testthat)
library(Matrix)

# ==============================================================================

test_that("get_example_snpdata works correctly", {
    # Setup - Load example data
    snp_data <- get_example_snpdata()
    
    # Test basic object structure
    expect_s4_class(snp_data, "SNPData")
    expect_true(nrow(snp_data) > 0)
    expect_true(ncol(snp_data) > 0)
    
    # Test SNP info structure
    snp_info <- get_snp_info(snp_data)
    expected_snp_cols <- c("snp_id", "chrom", "pos")
    expect_true(all(expected_snp_cols %in% colnames(snp_info)))
    
    # Test sample info structure
    sample_info <- get_sample_info(snp_data)
    expect_true("cell_id" %in% colnames(sample_info))
})

test_that("read_vcf_base works correctly", {
    # Setup - Get example VCF file
    vcf_file <- system.file("extdata/example_snpdata/cellSNP.base.vcf.gz", package = "snplet")
    skip_if_not(file.exists(vcf_file), "Example VCF file not found")
    
    vcf_data <- read_vcf_base(vcf_file)
    
    # Test return type and structure
    expect_s3_class(vcf_data, "data.frame")
    
    # Test required columns
    expected_cols <- c("snp_id", "chrom", "pos", "id", "ref", "alt", "qual", "filter", "info")
    expect_true(all(expected_cols %in% colnames(vcf_data)))
    expect_equal(colnames(vcf_data)[1], "snp_id")
    
    # Test SNP ID generation
    expect_true(all(grepl("^snp_", vcf_data$snp_id)))
    expect_equal(vcf_data$snp_id[1], "snp_1")
    
    # Test data types
    expect_type(vcf_data$chrom, "character")
    expect_type(vcf_data$pos, "integer")
    expect_type(vcf_data$ref, "character")
    expect_type(vcf_data$alt, "character")
})

test_that("merge_cell_annotations works correctly", {
    # Setup - Create test data
    donor_info <- data.frame(
        cell = c("CELL1", "CELL2", "CELL3"),
        donor_id = c("donor1", "donor2", "donor1"),
        stringsAsFactors = FALSE
    )
    
    vdj_info <- data.frame(
        barcode = c("CELL1", "CELL2", "CELL4"),
        raw_clonotype_id = c("clonotype1", "clonotype2", "clonotype1"),
        other_col = c("A", "B", "C"),
        stringsAsFactors = FALSE
    )
    
    # Execute merge
    result <- merge_cell_annotations(
        donor_info = donor_info,
        vdj_info = vdj_info,
        barcode_column = "barcode",
        clonotype_column = "raw_clonotype_id"
    )
    
    # Test return structure
    expect_s3_class(result, "data.frame")
    required_cols <- c("cell_id", "donor", "clonotype")
    expect_true(all(required_cols %in% colnames(result)))
    
    # Test merge behavior (left join on donor_info)
    expect_equal(sort(result$barcode), c("CELL1", "CELL2", "CELL3"))
    expect_equal(nrow(result), 3)
    
    # Test specific merge results
    cell1_row <- result[result$barcode == "CELL1", ]
    expect_equal(cell1_row$donor, "donor1")
    expect_equal(cell1_row$clonotype, "clonotype1")
    
    cell3_row <- result[result$barcode == "CELL3", ]
    expect_equal(cell3_row$donor, "donor1")
    expect_true(is.na(cell3_row$clonotype))  # Not in VDJ data
})

test_that("import_cellsnp works with example data", {
    # Setup - Get example data file paths
    cellsnp_dir <- system.file("extdata/example_snpdata", package = "snplet")
    vdj_file <- system.file("extdata/example_snpdata/filtered_contig_annotations.csv", package = "snplet")
    gene_anno_file <- system.file("extdata/example_gene_anno.tsv", package = "snplet")
    vireo_file <- system.file("extdata/example_snpdata/donor_ids.tsv", package = "snplet")
    
    required_files <- c(cellsnp_dir, vdj_file, gene_anno_file, vireo_file)
    skip_if_not(all(file.exists(required_files)), "Example data files not found")
    
    # Load gene annotation
    gene_annotation <- readr::read_tsv(gene_anno_file, show_col_types = FALSE)
    
    # Execute import
    snp_data <- import_cellsnp(
        cellsnp_dir = cellsnp_dir,
        vdj_file = vdj_file,
        gene_annotation = gene_annotation,
        vireo_file = vireo_file
    )
    
    # Test return object
    expect_s4_class(snp_data, "SNPData")
    expect_true(nrow(snp_data) > 0)
    expect_true(ncol(snp_data) > 0)
    
    # Test matrix dimension consistency
    expect_equal(dim(ref_count(snp_data)), dim(alt_count(snp_data)))
    expect_equal(dim(ref_count(snp_data)), dim(oth_count(snp_data)))
    
    # Test metadata structure
    snp_info <- get_snp_info(snp_data)
    sample_info <- get_sample_info(snp_data)
    expect_equal(nrow(snp_info), nrow(snp_data))
    expect_equal(nrow(sample_info), ncol(snp_data))
    
    # Test required columns
    expected_snp_cols <- c("snp_id", "chrom", "pos", "ref", "alt")
    expected_sample_cols <- c("cell_id", "donor", "clonotype")
    expect_true(all(expected_snp_cols %in% colnames(snp_info)))
    expect_true(all(expected_sample_cols %in% colnames(sample_info)))
})

test_that("import_cellsnp validates gene_annotation input", {
    # Setup - Create invalid gene annotation (missing required columns)
    invalid_gene_anno <- data.frame(
        gene_id = c("gene1", "gene2"),
        other_col = c("A", "B")
    )
    
    # Test validation error
    expect_error(
        import_cellsnp(
            cellsnp_dir = "dummy",
            vdj_file = "dummy", 
            gene_annotation = invalid_gene_anno
        ),
        "gene_annotation is missing required columns"
    )
})

test_that("export_cellsnp creates output files", {
    # Setup - Get example data and create temp directory
    snp_data <- get_example_snpdata()
    temp_dir <- tempdir()
    out_dir <- file.path(temp_dir, "test_export")
    
    # Execute export
    export_cellsnp(snp_data, out_dir)
    
    # Test that all expected files were created
    expected_files <- c(
        "cellSNP.tag.AD.mtx",
        "cellSNP.tag.DP.mtx", 
        "cellSNP.tag.OTH.mtx",
        "cellSNP.base.vcf.gz",
        "donor_ids.tsv",
        "filtered_contig_annotations.csv"
    )
    
    for (file in expected_files) {
        expect_true(file.exists(file.path(out_dir, file)))
    }
    
    # Cleanup
    unlink(out_dir, recursive = TRUE)
})
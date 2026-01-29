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
    # Verify get_example_snpdata returns a valid SNPData object
    expect_s4_class(snp_data, "SNPData")
    # Verify example data has SNPs (rows > 0)
    expect_true(nrow(snp_data) > 0)
    # Verify example data has samples (columns > 0)
    expect_true(ncol(snp_data) > 0)

    # Test SNP info structure
    snp_info <- get_snp_info(snp_data)
    expected_snp_cols <- c("snp_id", "chrom", "pos")
    # Verify essential SNP info columns are present
    expect_true(all(expected_snp_cols %in% colnames(snp_info)))

    # Test sample info structure
    barcode_info <- get_barcode_info(snp_data)
    # Verify cell_id column is present in sample info
    expect_true("cell_id" %in% colnames(barcode_info))
})

test_that("read_vcf_base works correctly", {
    # Setup - Get example VCF file
    vcf_file <- system.file("extdata/example_snpdata/cellSNP.base.vcf.gz", package = "snplet")
    skip_if_not(file.exists(vcf_file), "Example VCF file not found")

    vcf_data <- read_vcf_base(vcf_file)

    # Test return type and structure
    # Verify read_vcf_base returns a data frame
    expect_s3_class(vcf_data, "data.frame")

    # Test required columns
    expected_cols <- c("snp_id", "chrom", "pos", "id", "ref", "alt", "qual", "filter", "info")
    # Verify all expected VCF columns are present
    expect_true(all(expected_cols %in% colnames(vcf_data)))
    # Verify snp_id is the first column
    expect_equal(colnames(vcf_data)[1], "snp_id")

    # Test SNP ID generation
    # Verify all SNP IDs follow standardized format "chr:pos:ref:alt"
    expect_true(all(grepl("^chr.+:[0-9]+:[ACGT]+:[ACGT]+$", vcf_data$snp_id)))
    # Verify first SNP ID matches expected format from VCF data
    expect_equal(vcf_data$snp_id[1], "chr1:1108258:A:C")

    # Test data types
    # Verify chromosome column is character type
    expect_type(vcf_data$chrom, "character")
    # Verify position column is integer type
    expect_type(vcf_data$pos, "integer")
    # Verify reference allele column is character type
    expect_type(vcf_data$ref, "character")
    # Verify alternate allele column is character type
    expect_type(vcf_data$alt, "character")
})

test_that("merge_cell_annotations works correctly with standard column names", {
    # Setup - Create test data with standard column names
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
    # Verify merge_cell_annotations returns a data frame
    expect_s3_class(result, "data.frame")
    required_cols <- c("cell_id", "donor", "clonotype")
    # Verify all required columns are present after merge
    expect_true(all(required_cols %in% colnames(result)))

    # Test merge behavior (left join on donor_info)
    # Verify only donor_info cells are included (left join behavior)
    expect_equal(sort(result$barcode), c("CELL1", "CELL2", "CELL3"))
    # Verify cell_id values are generated sequentially
    expect_equal(sort(result$cell_id), c("cell_1", "cell_2", "cell_3"))
    # Verify result has expected number of rows
    expect_equal(nrow(result), 3)

    # Test specific merge results
    cell1_row <- result[result$barcode == "CELL1", ]
    # Verify CELL1 has correct donor assignment
    expect_equal(cell1_row$donor, "donor1")
    # Verify CELL1 has correct clonotype assignment from VDJ merge
    expect_equal(cell1_row$clonotype, "clonotype1")

    cell3_row <- result[result$barcode == "CELL3", ]
    # Verify CELL3 has correct donor assignment
    expect_equal(cell3_row$donor, "donor1")
    # Verify CELL3 has NA clonotype (not present in VDJ data)
    expect_true(is.na(cell3_row$clonotype)) # Not in VDJ data
})

test_that("merge_cell_annotations handles cell_id column naming", {
    # Setup - Create test data with cell_id instead of cell
    donor_info <- data.frame(
        cell_id = c("CELL1", "CELL2", "CELL3"),
        donor_id = c("donor1", "donor2", "donor1"),
        stringsAsFactors = FALSE
    )

    vdj_info <- data.frame(
        barcode = c("CELL1", "CELL2", "CELL4"),
        raw_clonotype_id = c("clonotype1", "clonotype2", "clonotype1"),
        stringsAsFactors = FALSE
    )

    # Execute merge
    result <- merge_cell_annotations(
        donor_info = donor_info,
        vdj_info = vdj_info,
        barcode_column = "barcode",
        clonotype_column = "raw_clonotype_id"
    )

    # Test that cell_id column was renamed to barcode for internal processing
    # Verify merge_cell_annotations returns a data frame
    expect_s3_class(result, "data.frame")
    # Verify all required columns are present after merge
    expect_true(all(c("cell_id", "donor", "clonotype") %in% colnames(result)))
    # Verify correct number of rows after merge
    expect_equal(nrow(result), 3)

    # Test specific merge results
    cell1_row <- result[result$barcode == "CELL1", ]
    # Verify CELL1 has correct donor assignment with cell_id input
    expect_equal(cell1_row$donor, "donor1")
    # Verify CELL1 has correct clonotype assignment with cell_id input
    expect_equal(cell1_row$clonotype, "clonotype1")
})

test_that("merge_cell_annotations handles donor column without donor_id", {
    # Setup - Create test data with donor instead of donor_id
    donor_info <- data.frame(
        cell = c("CELL1", "CELL2", "CELL3"),
        donor = c("donor1", "donor2", "donor1"),
        stringsAsFactors = FALSE
    )

    vdj_info <- data.frame(
        barcode = c("CELL1", "CELL2", "CELL4"),
        raw_clonotype_id = c("clonotype1", "clonotype2", "clonotype1"),
        stringsAsFactors = FALSE
    )

    # Execute merge
    result <- merge_cell_annotations(
        donor_info = donor_info,
        vdj_info = vdj_info,
        barcode_column = "barcode",
        clonotype_column = "raw_clonotype_id"
    )

    # Test that donor column is preserved without renaming
    # Verify merge_cell_annotations returns a data frame
    expect_s3_class(result, "data.frame")
    # Verify all required columns are present after merge
    expect_true(all(c("cell_id", "donor", "clonotype") %in% colnames(result)))

    # Test specific merge results
    cell1_row <- result[result$barcode == "CELL1", ]
    # Verify CELL1 has correct donor assignment when donor column already exists
    expect_equal(cell1_row$donor, "donor1")
})

test_that("merge_cell_annotations handles mixed column naming scenarios", {
    # Setup - Create test data with cell_id and no donor_id (only donor)
    donor_info <- data.frame(
        cell_id = c("CELL1", "CELL2", "CELL3"),
        donor = c("donor1", "donor2", "donor1"),
        stringsAsFactors = FALSE
    )

    vdj_info <- data.frame(
        barcode = c("CELL1", "CELL2", "CELL4"),
        raw_clonotype_id = c("clonotype1", "clonotype2", "clonotype1"),
        stringsAsFactors = FALSE
    )

    # Execute merge
    result <- merge_cell_annotations(
        donor_info = donor_info,
        vdj_info = vdj_info,
        barcode_column = "barcode",
        clonotype_column = "raw_clonotype_id"
    )

    # Test mixed column naming scenario
    # Verify merge_cell_annotations returns a data frame
    expect_s3_class(result, "data.frame")
    # Verify all required columns are present after merge
    expect_true(all(c("cell_id", "donor", "clonotype") %in% colnames(result)))
    # Verify correct number of rows after merge
    expect_equal(nrow(result), 3)

    # Test specific merge results
    cell2_row <- result[result$barcode == "CELL2", ]
    # Verify CELL2 has correct donor assignment in mixed naming scenario
    expect_equal(cell2_row$donor, "donor2")
    # Verify CELL2 has correct clonotype assignment in mixed naming scenario
    expect_equal(cell2_row$clonotype, "clonotype2")
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
    # Verify import completes without error
    snp_data <- expect_no_error(import_cellsnp(
        cellsnp_dir = cellsnp_dir,
        gene_annotation = gene_annotation,
        vdj_file = vdj_file,
        vireo_file = vireo_file
    ))

    # Test return object
    # Verify import_cellsnp returns a valid SNPData object
    expect_s4_class(snp_data, "SNPData")
    # Verify imported data has SNPs (rows > 0)
    expect_true(nrow(snp_data) > 0)
    # Verify imported data has samples (columns > 0)
    expect_true(ncol(snp_data) > 0)

    # Test matrix dimension consistency
    # Verify ref_count and alt_count matrices have same dimensions
    expect_equal(dim(ref_count(snp_data)), dim(alt_count(snp_data)))
    # Verify ref_count and oth_count matrices have same dimensions
    expect_equal(dim(ref_count(snp_data)), dim(oth_count(snp_data)))

    # Test metadata structure
    snp_info <- get_snp_info(snp_data)
    barcode_info <- get_barcode_info(snp_data)
    # Verify SNP info rows match matrix rows
    expect_equal(nrow(snp_info), nrow(snp_data))
    # Verify sample info rows match matrix columns
    expect_equal(nrow(barcode_info), ncol(snp_data))

    # Test required columns
    expected_snp_cols <- c("snp_id", "chrom", "pos", "ref", "alt")
    expected_sample_cols <- c("cell_id", "donor", "clonotype")
    # Verify all expected SNP info columns are present
    expect_true(all(expected_snp_cols %in% colnames(snp_info)))
    # Verify all expected sample info columns are present
    expect_true(all(expected_sample_cols %in% colnames(barcode_info)))
})

test_that("import_cellsnp validates gene_annotation input", {
    # Setup - Create invalid gene annotation (missing required columns)
    invalid_gene_anno <- data.frame(
        gene_id = c("gene1", "gene2"),
        other_col = c("A", "B")
    )

    # Test validation error
    # Verify error when gene_annotation is missing required columns
    expect_error(
        import_cellsnp(
            cellsnp_dir = "dummy",
            gene_annotation = invalid_gene_anno,
            vdj_file = "dummy"
        ),
        "gene_annotation is missing required columns"
    )
})

test_that("merge_cell_annotations works without VDJ info", {
    # Setup - Create donor info only, no VDJ info
    donor_info <- data.frame(
        cell = c("CELL1", "CELL2", "CELL3"),
        donor_id = c("donor1", "donor2", "donor1"),
        stringsAsFactors = FALSE
    )

    # Execute merge without VDJ info
    result <- merge_cell_annotations(
        donor_info = donor_info,
        vdj_info = NULL
    )

    # Test return structure
    # Verify merge_cell_annotations returns a data frame when VDJ info is NULL
    expect_s3_class(result, "data.frame")
    required_cols <- c("cell_id", "donor", "clonotype")
    # Verify all required columns are present including clonotype (as NA)
    expect_true(all(required_cols %in% colnames(result)))

    # Test that clonotype is all NA
    # Verify all clonotype values are NA when no VDJ info provided
    expect_true(all(is.na(result$clonotype)))

    # Test merge behavior
    # Verify all donor cells are included
    expect_equal(sort(result$barcode), c("CELL1", "CELL2", "CELL3"))
    # Verify cell_id values are generated sequentially
    expect_equal(sort(result$cell_id), c("cell_1", "cell_2", "cell_3"))
    # Verify result has expected number of rows
    expect_equal(nrow(result), 3)

    # Test donor assignments
    cell1_row <- result[result$barcode == "CELL1", ]
    # Verify CELL1 has correct donor assignment without VDJ
    expect_equal(cell1_row$donor, "donor1")
})

test_that("import_cellsnp works without vdj_file", {
    # Setup - Get example data file paths without VDJ
    cellsnp_dir <- system.file("extdata/example_snpdata", package = "snplet")
    gene_anno_file <- system.file("extdata/example_gene_anno.tsv", package = "snplet")
    vireo_file <- system.file("extdata/example_snpdata/donor_ids.tsv", package = "snplet")

    required_files <- c(cellsnp_dir, gene_anno_file, vireo_file)
    skip_if_not(all(file.exists(required_files)), "Example data files not found")

    # Load gene annotation
    gene_annotation <- readr::read_tsv(gene_anno_file, show_col_types = FALSE)

    # Execute import without VDJ file
    # Verify import completes without error
    snp_data <- expect_no_error(import_cellsnp(
        cellsnp_dir = cellsnp_dir,
        gene_annotation = gene_annotation,
        vireo_file = vireo_file
    ))

    # Test return object
    # Verify import_cellsnp returns a valid SNPData object without VDJ
    expect_s4_class(snp_data, "SNPData")
    # Verify imported data has SNPs (rows > 0)
    expect_true(nrow(snp_data) > 0)
    # Verify imported data has samples (columns > 0)
    expect_true(ncol(snp_data) > 0)

    # Test metadata structure
    barcode_info <- get_barcode_info(snp_data)
    # Verify clonotype column exists
    expect_true("clonotype" %in% colnames(barcode_info))
    # Verify all clonotype values are NA when no VDJ file provided
    expect_true(all(is.na(barcode_info$clonotype)))

    # Test required columns
    expected_sample_cols <- c("cell_id", "donor")
    # Verify expected sample info columns are present
    expect_true(all(expected_sample_cols %in% colnames(barcode_info)))
})

test_that("import_cellsnp works without vdj_file or vireo_file", {
    # Setup - Get example data file paths with minimal inputs
    cellsnp_dir <- system.file("extdata/example_snpdata", package = "snplet")
    gene_anno_file <- system.file("extdata/example_gene_anno.tsv", package = "snplet")

    required_files <- c(cellsnp_dir, gene_anno_file)
    skip_if_not(all(file.exists(required_files)), "Example data files not found")

    # Load gene annotation
    gene_annotation <- readr::read_tsv(gene_anno_file, show_col_types = FALSE)

    # Execute import without VDJ or Vireo files
    # Verify import completes without error
    snp_data <- expect_no_error(import_cellsnp(
        cellsnp_dir = cellsnp_dir,
        gene_annotation = gene_annotation
    ))

    # Test return object
    # Verify import_cellsnp returns a valid SNPData object with minimal inputs
    expect_s4_class(snp_data, "SNPData")
    # Verify imported data has SNPs (rows > 0)
    expect_true(nrow(snp_data) > 0)
    # Verify imported data has samples (columns > 0)
    expect_true(ncol(snp_data) > 0)

    # Test metadata structure
    barcode_info <- get_barcode_info(snp_data)
    # Verify clonotype column exists even without VDJ
    expect_true("clonotype" %in% colnames(barcode_info))
    # Verify donor column exists even without Vireo
    expect_true("donor" %in% colnames(barcode_info))
    # Verify all clonotype values are NA when no VDJ file provided
    expect_true(all(is.na(barcode_info$clonotype)))
    # Verify all donor values are "donor0" when no Vireo file provided
    expect_true(all(barcode_info$donor == "donor0"))
})

test_that("export_cellsnp creates output files", {
    # Setup - Get example data and create temp directory
    snp_data <- get_example_snpdata()
    out_dir <- withr::local_tempdir()

    # Execute export
    # Verify export completes without error
    expect_no_error(export_cellsnp(snp_data, out_dir))

    # Test that all expected files were created
    expected_files <- c(
        "cellSNP.tag.AD.mtx",
        "cellSNP.tag.DP.mtx",
        "cellSNP.tag.OTH.mtx",
        "cellSNP.base.vcf.gz",
        "donor_ids.tsv",
        "filtered_contig_annotations.csv"
    )

    # Verify all expected output files were created
    expect_true(file.exists(file.path(out_dir, "cellSNP.tag.AD.mtx")))
    expect_true(file.exists(file.path(out_dir, "cellSNP.tag.DP.mtx")))
    expect_true(file.exists(file.path(out_dir, "cellSNP.tag.OTH.mtx")))
    expect_true(file.exists(file.path(out_dir, "cellSNP.base.vcf.gz")))
    expect_true(file.exists(file.path(out_dir, "donor_ids.tsv")))
    expect_true(file.exists(file.path(out_dir, "filtered_contig_annotations.csv")))
})

test_that("export_cellsnp skips VDJ export when clonotype missing", {
    # Setup - Get example data without VDJ and create temp directory
    cellsnp_dir <- system.file("extdata/example_snpdata", package = "snplet")
    gene_anno_file <- system.file("extdata/example_gene_anno.tsv", package = "snplet")

    required_files <- c(cellsnp_dir, gene_anno_file)
    skip_if_not(all(file.exists(required_files)), "Example data files not found")

    gene_annotation <- readr::read_tsv(gene_anno_file, show_col_types = FALSE)
    snp_data <- import_cellsnp(
        cellsnp_dir = cellsnp_dir,
        gene_annotation = gene_annotation
    )

    out_dir <- withr::local_tempdir()

    # Execute export
    # Verify export completes without error
    expect_no_error(export_cellsnp(snp_data, out_dir))

    # Test that VDJ file was NOT created
    vdj_file <- file.path(out_dir, "filtered_contig_annotations.csv")
    # Verify VDJ file is not created when clonotype info missing
    expect_false(file.exists(vdj_file))

    # Verify other expected output files were created
    expect_true(file.exists(file.path(out_dir, "cellSNP.tag.AD.mtx")))
    expect_true(file.exists(file.path(out_dir, "cellSNP.tag.DP.mtx")))
    expect_true(file.exists(file.path(out_dir, "cellSNP.tag.OTH.mtx")))
    expect_true(file.exists(file.path(out_dir, "cellSNP.base.vcf.gz")))
    expect_true(file.exists(file.path(out_dir, "donor_ids.tsv")))
})

# ==============================================================================
# Integration Tests for Optional VDJ Workflow
# ==============================================================================

test_that("complete workflow: import without VDJ then add clonotype data", {
    # Setup - Get example data file paths
    cellsnp_dir <- system.file("extdata/example_snpdata", package = "snplet")
    gene_anno_file <- system.file("extdata/example_gene_anno.tsv", package = "snplet")
    vireo_file <- system.file("extdata/example_snpdata/donor_ids.tsv", package = "snplet")

    required_files <- c(cellsnp_dir, gene_anno_file, vireo_file)
    skip_if_not(all(file.exists(required_files)), "Example data files not found")

    # Step 1: Import without VDJ
    gene_annotation <- readr::read_tsv(gene_anno_file, show_col_types = FALSE)
    snp_data <- import_cellsnp(
        cellsnp_dir = cellsnp_dir,
        gene_annotation = gene_annotation,
        vireo_file = vireo_file
    )

    # Verify initial state - clonotype exists but all NA
    barcode_info <- get_barcode_info(snp_data)
    # Verify clonotype column exists
    expect_true("clonotype" %in% colnames(barcode_info))
    # Verify all clonotype values are initially NA
    expect_true(all(is.na(barcode_info$clonotype)))

    # Verify clonotype functions error appropriately
    # Confirm error when trying to use clonotype functions without data
    expect_error(
        clonotype_count_df(snp_data),
        "All clonotype values are NA"
    )

    # Step 2: Add clonotype information manually
    # Create mock clonotype data for the first few cells
    clonotype_data <- data.frame(
        cell_id = barcode_info$cell_id[1:min(10, nrow(barcode_info))],
        clonotype = paste0("clonotype_", 1:min(10, nrow(barcode_info))),
        stringsAsFactors = FALSE
    )

    snp_data_with_clonotype <- add_barcode_metadata(
        snp_data,
        clonotype_data,
        join_by = "cell_id",
        overwrite = TRUE
    )

    # Verify clonotype data was added
    barcode_info_updated <- get_barcode_info(snp_data_with_clonotype)
    # Verify clonotype column still exists
    expect_true("clonotype" %in% colnames(barcode_info_updated))
    # Verify not all clonotypes are NA anymore
    expect_false(all(is.na(barcode_info_updated$clonotype)))

    # Verify first cells have clonotype data
    first_cells <- barcode_info_updated[1:min(10, nrow(barcode_info_updated)), ]
    # Verify first 10 cells now have clonotype assignments
    expect_true(all(!is.na(first_cells$clonotype)))

    # Step 3: Use clonotype functions on subset with clonotype data
    # Filter to cells with clonotype data
    snp_data_filtered <- filter_barcodes(
        snp_data_with_clonotype,
        !is.na(clonotype)
    )

    # Verify clonotype functions now work
    result <- clonotype_count_df(snp_data_filtered, test_maf = FALSE)
    # Verify clonotype_count_df returns valid data frame
    expect_s3_class(result, "data.frame")
    # Verify clonotype column is present in results
    expect_true("clonotype" %in% colnames(result))
    # Verify results are not empty
    expect_true(nrow(result) > 0)

    # Verify to_expr_matrix works
    expr_mat <- to_expr_matrix(snp_data_filtered, level = "clonotype")
    # Verify to_expr_matrix returns a matrix
    expect_true(is.matrix(expr_mat))
    # Verify matrix has columns (clonotypes)
    expect_true(ncol(expr_mat) > 0)
})

test_that("import with VDJ then export and re-import preserves clonotype", {
    # Setup - Import with VDJ
    cellsnp_dir <- system.file("extdata/example_snpdata", package = "snplet")
    vdj_file <- system.file("extdata/example_snpdata/filtered_contig_annotations.csv", package = "snplet")
    gene_anno_file <- system.file("extdata/example_gene_anno.tsv", package = "snplet")
    vireo_file <- system.file("extdata/example_snpdata/donor_ids.tsv", package = "snplet")

    required_files <- c(cellsnp_dir, vdj_file, gene_anno_file, vireo_file)
    skip_if_not(all(file.exists(required_files)), "Example data files not found")

    gene_annotation <- readr::read_tsv(gene_anno_file, show_col_types = FALSE)

    # Import with VDJ
    snp_data_original <- import_cellsnp(
        cellsnp_dir = cellsnp_dir,
        gene_annotation = gene_annotation,
        vdj_file = vdj_file,
        vireo_file = vireo_file
    )

    # Verify clonotype data present
    barcode_info_original <- get_barcode_info(snp_data_original)
    # Verify clonotype column exists in imported data
    expect_true("clonotype" %in% colnames(barcode_info_original))
    clonotypes_with_data <- sum(!is.na(barcode_info_original$clonotype))
    # Verify some cells have clonotype assignments
    expect_true(clonotypes_with_data > 0)

    # Export to temp directory
    out_dir <- withr::local_tempdir()
    # Verify export completes without error
    expect_no_error(export_cellsnp(snp_data_original, out_dir))

    # Verify VDJ file was created
    vdj_exported <- file.path(out_dir, "filtered_contig_annotations.csv")
    # Verify VDJ file was exported when clonotype data present
    expect_true(file.exists(vdj_exported))

    # Re-import
    snp_data_reimported <- import_cellsnp(
        cellsnp_dir = out_dir,
        gene_annotation = gene_annotation,
        vdj_file = vdj_exported,
        vireo_file = file.path(out_dir, "donor_ids.tsv")
    )

    # Verify clonotype data preserved
    barcode_info_reimported <- get_barcode_info(snp_data_reimported)
    # Verify clonotype column exists after re-import
    expect_true("clonotype" %in% colnames(barcode_info_reimported))

    # Compare number of cells with clonotype data
    clonotypes_reimported <- sum(!is.na(barcode_info_reimported$clonotype))
    # Verify same number of cells have clonotype data after roundtrip
    expect_equal(clonotypes_reimported, clonotypes_with_data)
})

test_that("import without VDJ, export, re-import maintains no clonotype state", {
    # Setup - Import without VDJ
    cellsnp_dir <- system.file("extdata/example_snpdata", package = "snplet")
    gene_anno_file <- system.file("extdata/example_gene_anno.tsv", package = "snplet")

    required_files <- c(cellsnp_dir, gene_anno_file)
    skip_if_not(all(file.exists(required_files)), "Example data files not found")

    gene_annotation <- readr::read_tsv(gene_anno_file, show_col_types = FALSE)

    # Import without VDJ
    snp_data_original <- import_cellsnp(
        cellsnp_dir = cellsnp_dir,
        gene_annotation = gene_annotation
    )

    # Verify all clonotypes are NA
    barcode_info_original <- get_barcode_info(snp_data_original)
    # Verify all clonotype values are NA when no VDJ imported
    expect_true(all(is.na(barcode_info_original$clonotype)))

    # Export to temp directory
    out_dir <- withr::local_tempdir()
    # Verify export completes without error
    expect_no_error(export_cellsnp(snp_data_original, out_dir))

    # Verify VDJ file was NOT created
    vdj_exported <- file.path(out_dir, "filtered_contig_annotations.csv")
    # Verify VDJ file is not exported when no clonotype data
    expect_false(file.exists(vdj_exported))

    # Re-import without VDJ file
    snp_data_reimported <- import_cellsnp(
        cellsnp_dir = out_dir,
        gene_annotation = gene_annotation,
        vireo_file = file.path(out_dir, "donor_ids.tsv")
    )

    # Verify all clonotypes still NA
    barcode_info_reimported <- get_barcode_info(snp_data_reimported)
    # Verify clonotype column exists after re-import
    expect_true("clonotype" %in% colnames(barcode_info_reimported))
    # Verify all clonotype values remain NA after roundtrip
    expect_true(all(is.na(barcode_info_reimported$clonotype)))
})

# ==============================================================================
# Error Handling Tests
# ==============================================================================

test_that("merge_cell_annotations errors when barcode_column not found in vdj_info", {
    # Setup - Create donor and VDJ info with mismatched column name
    donor_info <- data.frame(
        cell = c("CELL1", "CELL2", "CELL3"),
        donor_id = c("donor1", "donor2", "donor1"),
        stringsAsFactors = FALSE
    )

    vdj_info <- data.frame(
        wrong_barcode_name = c("CELL1", "CELL2", "CELL4"),
        raw_clonotype_id = c("clonotype1", "clonotype2", "clonotype1"),
        stringsAsFactors = FALSE
    )

    # Verify error when barcode_column is not found in vdj_info
    expect_error(
        merge_cell_annotations(
            donor_info = donor_info,
            vdj_info = vdj_info,
            barcode_column = "barcode",
            clonotype_column = "raw_clonotype_id"
        ),
        "Column barcode not found in VDJ annotation file"
    )
})

test_that("merge_cell_annotations errors when clonotype_column not found in vdj_info", {
    # Setup - Create donor and VDJ info with mismatched clonotype column name
    donor_info <- data.frame(
        cell = c("CELL1", "CELL2", "CELL3"),
        donor_id = c("donor1", "donor2", "donor1"),
        stringsAsFactors = FALSE
    )

    vdj_info <- data.frame(
        barcode = c("CELL1", "CELL2", "CELL4"),
        wrong_clonotype_name = c("clonotype1", "clonotype2", "clonotype1"),
        stringsAsFactors = FALSE
    )

    # Verify error when clonotype_column is not found in vdj_info
    expect_error(
        merge_cell_annotations(
            donor_info = donor_info,
            vdj_info = vdj_info,
            barcode_column = "barcode",
            clonotype_column = "raw_clonotype_id"
        ),
        "Column raw_clonotype_id not found in VDJ annotation file"
    )
})

test_that("merge_cell_annotations handles donor_info without donor column", {
    # Setup - Create donor info without donor or donor_id column
    donor_info <- data.frame(
        cell = c("CELL1", "CELL2", "CELL3"),
        other_column = c("A", "B", "C"),
        stringsAsFactors = FALSE
    )

    vdj_info <- data.frame(
        barcode = c("CELL1", "CELL2", "CELL4"),
        raw_clonotype_id = c("clonotype1", "clonotype2", "clonotype1"),
        stringsAsFactors = FALSE
    )

    # Execute merge
    result <- merge_cell_annotations(
        donor_info = donor_info,
        vdj_info = vdj_info,
        barcode_column = "barcode",
        clonotype_column = "raw_clonotype_id"
    )

    # Verify merge_cell_annotations returns a data frame
    expect_s3_class(result, "data.frame")
    # Verify donor column exists
    expect_true("donor" %in% colnames(result))
    # Verify all donor values are NA when donor column not in input
    expect_true(all(is.na(result$donor)))
    # Verify clonotype column exists and has values
    expect_true("clonotype" %in% colnames(result))
    # Verify clonotypes are assigned correctly
    expect_equal(result$clonotype[result$barcode == "CELL1"], "clonotype1")
})

test_that("merge_cell_annotations handles empty vdj_info data frame", {
    # Setup - Create donor info and empty VDJ info
    donor_info <- data.frame(
        cell = c("CELL1", "CELL2", "CELL3"),
        donor_id = c("donor1", "donor2", "donor1"),
        stringsAsFactors = FALSE
    )

    # Create empty vdj_info with correct columns but no rows
    vdj_info <- data.frame(
        barcode = character(0),
        raw_clonotype_id = character(0),
        stringsAsFactors = FALSE
    )

    # Execute merge
    result <- merge_cell_annotations(
        donor_info = donor_info,
        vdj_info = vdj_info,
        barcode_column = "barcode",
        clonotype_column = "raw_clonotype_id"
    )

    # Verify merge_cell_annotations returns a data frame
    expect_s3_class(result, "data.frame")
    # Verify all required columns exist
    expect_true(all(c("cell_id", "barcode", "donor", "clonotype") %in% colnames(result)))
    # Verify all clonotypes are NA when vdj_info is empty
    expect_true(all(is.na(result$clonotype)))
    # Verify donors are assigned correctly
    expect_equal(result$donor[result$barcode == "CELL1"], "donor1")
})

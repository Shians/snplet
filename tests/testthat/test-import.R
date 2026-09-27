# ==============================================================================
# Test Suite: Data Import Functions
# Description: Tests for cellSNP data import, export, and utility functions
# ==============================================================================

library(testthat)
library(Matrix)

# ==============================================================================

test_that("get_example_snpdata returns a valid SNPData object with data", {
    snp_data <- get_example_snpdata()

    # Verify get_example_snpdata returns a valid SNPData object
    expect_s4_class(snp_data, "SNPData")
    # Verify example data has SNPs (rows > 0)
    expect_true(nrow(snp_data) > 0)
    # Verify example data has samples (columns > 0)
    expect_true(ncol(snp_data) > 0)
})

test_that("get_example_snpdata includes essential SNP info columns", {
    snp_data <- get_example_snpdata()

    snp_info <- snp_info(snp_data)
    expected_snp_cols <- c("snp_id", "chrom", "pos")

    # Verify essential SNP info columns are present
    expect_true(all(expected_snp_cols %in% colnames(snp_info)))
})

test_that("get_example_snpdata includes cell_id in sample info", {
    snp_data <- get_example_snpdata()

    barcode_info <- barcode_info(snp_data)

    # Verify cell_id column is present in sample info
    expect_true("cell_id" %in% colnames(barcode_info))
})

test_that("import_cellsnp() stores the gene annotation it assigned gene names from", {
    snp_data <- get_example_snpdata()
    gene_anno <- readr::read_tsv(
        system.file("extdata/example_gene_anno.tsv", package = "snplet"),
        show_col_types = FALSE
    )
    stored <- gene_anno(snp_data)

    # Verify every gene of the import annotation is stored, strand included
    expect_setequal(stored$gene_name, gene_anno$gene_name)
    expect_true("strand" %in% colnames(stored))
    # Check the GFF attribute strings are not carried onto the object
    expect_false("attributes" %in% colnames(stored))
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

# Shared fixture paths and import call for the "works with example data" tests below
import_example_snpdata <- function(library_id = "lib1") {
    cellsnp_dir <- system.file("extdata/example_snpdata", package = "snplet")
    vdj_file <- system.file("extdata/example_snpdata/filtered_contig_annotations.csv", package = "snplet")
    gene_anno_file <- system.file("extdata/example_gene_anno.tsv", package = "snplet")
    vireo_folder <- system.file("extdata/example_snpdata", package = "snplet")

    required_files <- c(cellsnp_dir, vdj_file, gene_anno_file, file.path(vireo_folder, "donor_ids.tsv"))
    skip_if_not(all(file.exists(required_files)), "Example data files not found")

    gene_annotation <- readr::read_tsv(gene_anno_file, show_col_types = FALSE)

    import_cellsnp(
        cellsnp_dir = cellsnp_dir,
        gene_annotation = gene_annotation,
        vdj_file = vdj_file,
        vireo_folder = vireo_folder,
        library_id = library_id
    )
}

test_that("import_cellsnp() labels every cell with the supplied library_id", {
    snp_data <- import_example_snpdata(library_id = "lib_A")

    # Verify one cellSNP run yields one library label across all its cells
    expect_equal(unique(barcode_info(snp_data)$library_id), "lib_A")
})

test_that("import_cellsnp() leaves library_id as NA when not supplied", {
    cellsnp_dir <- system.file("extdata/example_snpdata", package = "snplet")
    gene_annotation <- data.frame(chrom = "chr1", start = 1, end = 1e9, gene_name = "dummy", strand = "+")

    # Verify a single-library workflow can omit library_id entirely
    snp_data <- expect_no_error(import_cellsnp(cellsnp_dir, gene_annotation))
    # Check that every cell is left with an NA library_id rather than a guess
    expect_true(all(is.na(barcode_info(snp_data)$library_id)))
})

test_that("import_cellsnp() rejects a library_id that is not a single string", {
    cellsnp_dir <- system.file("extdata/example_snpdata", package = "snplet")
    gene_annotation <- data.frame(chrom = "chr1", start = 1, end = 1e9, gene_name = "dummy", strand = "+")

    # Check that a vector of labels is refused, since one run is one library
    expect_error(
        import_cellsnp(cellsnp_dir, gene_annotation, library_id = c("lib_A", "lib_B")),
        "single string"
    )
})

test_that("import_cellsnp completes without error and returns a populated SNPData object", {
    # Verify import completes without error
    snp_data <- expect_no_error(import_example_snpdata())

    # Verify import_cellsnp returns a valid SNPData object
    expect_s4_class(snp_data, "SNPData")
    # Verify imported data has SNPs (rows > 0)
    expect_true(nrow(snp_data) > 0)
    # Verify imported data has samples (columns > 0)
    expect_true(ncol(snp_data) > 0)
})

test_that("import_cellsnp returns matrices with consistent dimensions", {
    snp_data <- import_example_snpdata()

    # Verify ref_count and alt_count matrices have same dimensions
    expect_equal(dim(ref_count(snp_data)), dim(alt_count(snp_data)))
    # Verify ref_count and oth_count matrices have same dimensions
    expect_equal(dim(ref_count(snp_data)), dim(oth_count(snp_data)))
})

test_that("import_cellsnp metadata row counts match matrix dimensions", {
    snp_data <- import_example_snpdata()

    snp_info <- snp_info(snp_data)
    barcode_info <- barcode_info(snp_data)

    # Verify SNP info rows match matrix rows
    expect_equal(nrow(snp_info), nrow(snp_data))
    # Verify sample info rows match matrix columns
    expect_equal(nrow(barcode_info), ncol(snp_data))
})

test_that("import_cellsnp metadata includes all expected columns", {
    snp_data <- import_example_snpdata()

    snp_info <- snp_info(snp_data)
    barcode_info <- barcode_info(snp_data)
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
            vdj_file = "dummy",
            library_id = "lib1"
        ),
        "gene_annotation is missing required columns"
    )
})

test_that("import_cellsnp validates gene_annotation strand values before reading any files", {
    star_strand <- data.frame(chrom = "chr1", start = 1, end = 1e9, gene_name = "dummy", strand = "*")

    # Verify an unstranded gene is refused up front, naming the argument and the value
    expect_error(
        import_cellsnp(cellsnp_dir = "dummy", gene_annotation = star_strand, library_id = "lib1"),
        "gene_annotation\\$strand must be .* found: \"\\*\""
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

# Shared fixture for the "import without vdj_file" tests below
import_snpdata_without_vdj <- function() {
    cellsnp_dir <- system.file("extdata/example_snpdata", package = "snplet")
    gene_anno_file <- system.file("extdata/example_gene_anno.tsv", package = "snplet")
    vireo_folder <- system.file("extdata/example_snpdata", package = "snplet")

    required_files <- c(cellsnp_dir, gene_anno_file, file.path(vireo_folder, "donor_ids.tsv"))
    skip_if_not(all(file.exists(required_files)), "Example data files not found")

    gene_annotation <- readr::read_tsv(gene_anno_file, show_col_types = FALSE)

    import_cellsnp(
        cellsnp_dir = cellsnp_dir,
        gene_annotation = gene_annotation,
        vireo_folder = vireo_folder,
        library_id = "lib1",
    )
}

test_that("import_cellsnp without vdj_file completes without error and returns a populated SNPData object", {
    # Verify import completes without error
    snp_data <- expect_no_error(import_snpdata_without_vdj())

    # Verify import_cellsnp returns a valid SNPData object without VDJ
    expect_s4_class(snp_data, "SNPData")
    # Verify imported data has SNPs (rows > 0)
    expect_true(nrow(snp_data) > 0)
    # Verify imported data has samples (columns > 0)
    expect_true(ncol(snp_data) > 0)
})

test_that("import_cellsnp without vdj_file leaves clonotype as NA", {
    snp_data <- import_snpdata_without_vdj()

    barcode_info <- barcode_info(snp_data)

    # Verify clonotype column exists
    expect_true("clonotype" %in% colnames(barcode_info))
    # Verify all clonotype values are NA when no VDJ file provided
    expect_true(all(is.na(barcode_info$clonotype)))
})

test_that("import_cellsnp without vdj_file includes expected sample columns", {
    snp_data <- import_snpdata_without_vdj()

    barcode_info <- barcode_info(snp_data)
    expected_sample_cols <- c("cell_id", "donor")

    # Verify expected sample info columns are present
    expect_true(all(expected_sample_cols %in% colnames(barcode_info)))
})

# Shared fixture for the "import without vdj_file or vireo_folder" tests below
import_snpdata_without_vdj_or_vireo <- function() {
    cellsnp_dir <- system.file("extdata/example_snpdata", package = "snplet")
    gene_anno_file <- system.file("extdata/example_gene_anno.tsv", package = "snplet")

    required_files <- c(cellsnp_dir, gene_anno_file)
    skip_if_not(all(file.exists(required_files)), "Example data files not found")

    gene_annotation <- readr::read_tsv(gene_anno_file, show_col_types = FALSE)

    import_cellsnp(
        cellsnp_dir = cellsnp_dir,
        gene_annotation = gene_annotation,
        library_id = "lib1",
    )
}

test_that("import_cellsnp without vdj_file or vireo_folder completes without error and returns a populated SNPData object", {
    # Verify import completes without error
    snp_data <- expect_no_error(import_snpdata_without_vdj_or_vireo())

    # Verify import_cellsnp returns a valid SNPData object with minimal inputs
    expect_s4_class(snp_data, "SNPData")
    # Verify imported data has SNPs (rows > 0)
    expect_true(nrow(snp_data) > 0)
    # Verify imported data has samples (columns > 0)
    expect_true(ncol(snp_data) > 0)
})

test_that("import_cellsnp without vdj_file or vireo_folder defaults clonotype and donor columns", {
    snp_data <- import_snpdata_without_vdj_or_vireo()

    barcode_info <- barcode_info(snp_data)

    # Verify clonotype column exists even without VDJ
    expect_true("clonotype" %in% colnames(barcode_info))
    # Verify donor column exists even without Vireo
    expect_true("donor" %in% colnames(barcode_info))
    # Verify all clonotype values are NA when no VDJ file provided
    expect_true(all(is.na(barcode_info$clonotype)))
    # Verify all donor values are "donor0" when no Vireo file provided
    expect_true(all(barcode_info$donor == "donor0"))
})

test_that("import_cellsnp(vireo_folder=) reads donor assignments when the folder has no genotype VCF", {
    cellsnp_dir <- system.file("extdata/example_snpdata", package = "snplet")
    gene_anno_file <- system.file("extdata/example_gene_anno.tsv", package = "snplet")
    vireo_folder <- system.file("extdata/example_snpdata", package = "snplet")

    required_files <- c(cellsnp_dir, gene_anno_file, file.path(vireo_folder, "donor_ids.tsv"))
    skip_if_not(all(file.exists(required_files)), "Example data files not found")
    # This folder is expected to have no GT_donors.vireo.vcf.gz -- that's the case under test
    skip_if(file.exists(file.path(vireo_folder, "GT_donors.vireo.vcf.gz")), "Unexpectedly found a genotype VCF")

    gene_annotation <- readr::read_tsv(gene_anno_file, show_col_types = FALSE)
    snp_data <- import_cellsnp(
        cellsnp_dir = cellsnp_dir,
        gene_annotation = gene_annotation,
        vireo_folder = vireo_folder,
        library_id = "lib1",
    )

    # Verify donor assignments were still read from donor_ids.tsv
    expect_true("donor" %in% colnames(barcode_info(snp_data)))
    # Verify donor_snp_info stays empty since there was no genotype VCF to parse
    expect_equal(nrow(donor_snp_info(snp_data)), 0)
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

    # Verify AD matrix file was created
    expect_true(file.exists(file.path(out_dir, "cellSNP.tag.AD.mtx")))
    # Verify DP matrix file was created
    expect_true(file.exists(file.path(out_dir, "cellSNP.tag.DP.mtx")))
    # Verify OTH matrix file was created
    expect_true(file.exists(file.path(out_dir, "cellSNP.tag.OTH.mtx")))
    # Verify base VCF file was created
    expect_true(file.exists(file.path(out_dir, "cellSNP.base.vcf.gz")))
    # Verify donor IDs file was created
    expect_true(file.exists(file.path(out_dir, "donor_ids.tsv")))
    # Verify VDJ contig annotations file was created
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
        gene_annotation = gene_annotation,
        library_id = "lib1",
    )

    out_dir <- withr::local_tempdir()

    # Execute export
    # Verify export completes without error
    expect_no_error(export_cellsnp(snp_data, out_dir))

    # Test that VDJ file was NOT created
    vdj_file <- file.path(out_dir, "filtered_contig_annotations.csv")
    # Verify VDJ file is not created when clonotype info missing
    expect_false(file.exists(vdj_file))

    # Verify AD matrix file was created
    expect_true(file.exists(file.path(out_dir, "cellSNP.tag.AD.mtx")))
    # Verify DP matrix file was created
    expect_true(file.exists(file.path(out_dir, "cellSNP.tag.DP.mtx")))
    # Verify OTH matrix file was created
    expect_true(file.exists(file.path(out_dir, "cellSNP.tag.OTH.mtx")))
    # Verify base VCF file was created
    expect_true(file.exists(file.path(out_dir, "cellSNP.base.vcf.gz")))
    # Verify donor IDs file was created
    expect_true(file.exists(file.path(out_dir, "donor_ids.tsv")))
})

# ==============================================================================
# Integration Tests for Optional VDJ Workflow
# ==============================================================================

test_that("complete workflow: import without VDJ then add clonotype data", {
    # Setup - Get example data file paths
    cellsnp_dir <- system.file("extdata/example_snpdata", package = "snplet")
    gene_anno_file <- system.file("extdata/example_gene_anno.tsv", package = "snplet")
    vireo_folder <- system.file("extdata/example_snpdata", package = "snplet")

    required_files <- c(cellsnp_dir, gene_anno_file, file.path(vireo_folder, "donor_ids.tsv"))
    skip_if_not(all(file.exists(required_files)), "Example data files not found")

    # Step 1: Import without VDJ
    gene_annotation <- readr::read_tsv(gene_anno_file, show_col_types = FALSE)
    snp_data <- import_cellsnp(
        cellsnp_dir = cellsnp_dir,
        gene_annotation = gene_annotation,
        vireo_folder = vireo_folder,
        library_id = "lib1",
    )

    # Verify initial state - clonotype exists but all NA
    barcode_info <- barcode_info(snp_data)
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
    barcode_info_updated <- barcode_info(snp_data_with_clonotype)
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
    vireo_folder <- system.file("extdata/example_snpdata", package = "snplet")

    required_files <- c(cellsnp_dir, vdj_file, gene_anno_file, file.path(vireo_folder, "donor_ids.tsv"))
    skip_if_not(all(file.exists(required_files)), "Example data files not found")

    gene_annotation <- readr::read_tsv(gene_anno_file, show_col_types = FALSE)

    # Import with VDJ
    snp_data_original <- import_cellsnp(
        cellsnp_dir = cellsnp_dir,
        gene_annotation = gene_annotation,
        vdj_file = vdj_file,
        vireo_folder = vireo_folder,
        library_id = "lib1",
    )

    # Verify clonotype data present
    barcode_info_original <- barcode_info(snp_data_original)
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
        vireo_folder = out_dir,
        library_id = "lib1",
    )

    # Verify clonotype data preserved
    barcode_info_reimported <- barcode_info(snp_data_reimported)
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
        gene_annotation = gene_annotation,
        library_id = "lib1",
    )

    # Verify all clonotypes are NA
    barcode_info_original <- barcode_info(snp_data_original)
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
        vireo_folder = out_dir,
        library_id = "lib1",
    )

    # Verify all clonotypes still NA
    barcode_info_reimported <- barcode_info(snp_data_reimported)
    # Verify clonotype column exists after re-import
    expect_true("clonotype" %in% colnames(barcode_info_reimported))
    # Verify all clonotype values remain NA after roundtrip
    expect_true(all(is.na(barcode_info_reimported$clonotype)))
})

test_that("export_cellsnp writes one annotation row per cell, in matrix column order", {
    # Setup - Example data has clonotypes, so both annotation files are written
    snp_data <- get_example_snpdata()
    out_dir <- withr::local_tempdir()

    export_cellsnp(snp_data, out_dir)
    cells <- barcode_info(snp_data)

    donor_df <- readr::read_tsv(file.path(out_dir, "donor_ids.tsv"), show_col_types = FALSE)
    # Verify donor_ids.tsv has exactly one row per cell rather than distinct pairs
    expect_equal(nrow(donor_df), ncol(snp_data))
    # Check barcodes are written in matrix column order
    expect_equal(donor_df$cell, cells$barcode)
    # Verify donor labels stay aligned with their cells
    expect_equal(donor_df$donor_id, cells$donor)
    # Ensure the exported directory records which library it came from
    expect_equal(donor_df$library_id, cells$library_id)
    # Ensure the object's own cell identifiers survive the export
    expect_equal(donor_df$cell_id, cells$cell_id)

    vdj_df <- readr::read_csv(file.path(out_dir, "filtered_contig_annotations.csv"), show_col_types = FALSE)
    # Verify VDJ annotations are also one row per cell in column order
    expect_equal(vdj_df$barcode, cells$barcode)
    # Confirm clonotypes stay aligned with their cells
    expect_equal(vdj_df$raw_clonotype_id, cells$clonotype)

    samples <- readr::read_tsv(
        file.path(out_dir, "cellSNP.samples.tsv"),
        col_names = "barcode",
        show_col_types = FALSE
    )
    # Ensure cellSNP.samples.tsv stays single-column for external readers
    expect_equal(ncol(samples), 1L)
    # Verify the barcode list matches the matrix columns
    expect_equal(samples$barcode, cells$barcode)
})

test_that("export_cellsnp rejects a multi-library object", {
    # Setup - Relabel half the cells as a second library
    snp_data <- get_example_snpdata()
    cells <- barcode_info(snp_data)
    cells$library_id <- rep(c("lib1", "lib2"), length.out = nrow(cells))
    barcode_info(snp_data) <- cells

    out_dir <- withr::local_tempdir()

    # Verify export refuses rather than fusing barcodes shared across libraries
    expect_error(export_cellsnp(snp_data, out_dir), "cannot write a multi-library object")
    # Ensure the error names the libraries so the user knows what to split
    expect_error(export_cellsnp(snp_data, out_dir), "lib1, lib2")
    # Confirm nothing was written before the check failed
    expect_false(file.exists(file.path(out_dir, "cellSNP.tag.AD.mtx")))
})

test_that("export_cellsnp rejects repeated barcodes", {
    # Setup - Duplicate one barcode within a single library
    snp_data <- get_example_snpdata()
    cells <- barcode_info(snp_data)
    skip_if(nrow(cells) < 2, "Example data has too few cells")
    cells$barcode[2] <- cells$barcode[1]
    barcode_info(snp_data) <- cells

    out_dir <- withr::local_tempdir()

    # Verify export refuses barcodes the cellSNP format cannot tell apart
    expect_error(export_cellsnp(snp_data, out_dir), "cannot write repeated barcodes")
    # Ensure the offending barcode is named in the message
    expect_error(export_cellsnp(snp_data, out_dir), cells$barcode[1], fixed = TRUE)
})

# ==============================================================================
# Test: as_singlecellexperiment()
# ==============================================================================

test_that("as_singlecellexperiment() wraps counts as assays with matching dimnames", {
    skip_if_not_installed("SingleCellExperiment")
    snp_data <- get_example_snpdata()

    sce <- as_singlecellexperiment(snp_data)

    # Verify the container class and dimensions match the SNPData object
    expect_s4_class(sce, "SingleCellExperiment")
    expect_equal(dim(sce), dim(snp_data))
    # Verify ref/alt assays are present and identical to the source matrices
    expect_equal(as.matrix(SummarizedExperiment::assay(sce, "ref")), as.matrix(ref_count(snp_data)))
    expect_equal(as.matrix(SummarizedExperiment::assay(sce, "alt")), as.matrix(alt_count(snp_data)))
    # Confirm rows/columns are dimnamed by snp_id/cell_id, matching snp_info/barcode_info
    expect_equal(rownames(sce), snp_info(snp_data)$snp_id)
    expect_equal(colnames(sce), barcode_info(snp_data)$cell_id)
})

test_that("as_singlecellexperiment() omits an all-zero oth assay", {
    skip_if_not_installed("SingleCellExperiment")
    ref <- Matrix::Matrix(matrix(c(5L, 3L, 2L, 8L), 2, 2), sparse = TRUE)
    alt <- Matrix::Matrix(matrix(c(1L, 2L, 4L, 1L), 2, 2), sparse = TRUE)
    oth <- Matrix::Matrix(matrix(0L, 2, 2), sparse = TRUE)
    snp_data <- SNPData(
        ref_count = ref,
        alt_count = alt,
        oth_count = oth,
        snp_info = data.frame(chrom = "chr1", pos = c(1L, 2L), ref = "A", alt = "G"),
        barcode_info = data.frame(barcode = c("c1", "c2"))
    )

    sce <- as_singlecellexperiment(snp_data)

    # Verify an all-zero oth_count is left out rather than carried as dead weight
    expect_false("oth" %in% names(SummarizedExperiment::assays(sce)))
})

test_that("as_singlecellexperiment() keeps a non-zero oth assay", {
    skip_if_not_installed("SingleCellExperiment")
    ref <- Matrix::Matrix(matrix(c(5L, 3L, 2L, 8L), 2, 2), sparse = TRUE)
    alt <- Matrix::Matrix(matrix(c(1L, 2L, 4L, 1L), 2, 2), sparse = TRUE)
    oth <- Matrix::Matrix(matrix(c(0L, 0L, 1L, 0L), 2, 2), sparse = TRUE)
    snp_data <- SNPData(
        ref_count = ref,
        alt_count = alt,
        oth_count = oth,
        snp_info = data.frame(chrom = "chr1", pos = c(1L, 2L), ref = "A", alt = "G"),
        barcode_info = data.frame(barcode = c("c1", "c2"))
    )

    sce <- as_singlecellexperiment(snp_data)

    # Verify a genuinely non-zero oth_count is carried through as its own assay
    expect_true("oth" %in% names(SummarizedExperiment::assays(sce)))
    expect_equal(as.matrix(SummarizedExperiment::assay(sce, "oth")), as.matrix(oth_count(snp_data)))
})

test_that("as_singlecellexperiment() carries barcode_info/snp_info columns into colData/rowData", {
    skip_if_not_installed("SingleCellExperiment")
    snp_data <- get_example_snpdata()

    sce <- as_singlecellexperiment(snp_data)
    col_data <- SummarizedExperiment::colData(sce)
    row_data <- SummarizedExperiment::rowData(sce)

    # Verify barcode_info columns (other than the cell_id used as dimnames) ride along
    expect_true(all(setdiff(colnames(barcode_info(snp_data)), "cell_id") %in% colnames(col_data)))
    # Verify snp_info columns (other than the snp_id used as dimnames) ride along
    expect_true(all(setdiff(colnames(snp_info(snp_data)), "snp_id") %in% colnames(row_data)))
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

# ==============================================================================
# Test Suite: Vireo genotype ingestion
# Description: .read_vireo_gt() parses a Vireo GT VCF into per-(SNP, donor)
#              zygosity calls, and import_cellsnp(vireo_folder=) wires it in.
# ==============================================================================

# A small multi-sample GT VCF matching Vireo's GT_donors.vireo.vcf.gz layout
# (FORMAT = GT:AD:DP:PL, one sample column per donor).
#   snp1: donor0 het (0/1), donor1 hom-ref (0/0)
#   snp2: donor0 hom-alt (1/1) with a confident PL, donor1 missing PL
#   snp3: not present in snp_info, so it must be filtered out
write_test_vireo_gt_vcf <- function(path) {
    header_lines <- c("##fileformat=VCFv4.2", "##source=snplet-test")
    col_line <- "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tdonor0\tdonor1"
    data_lines <- c(
        "chr1\t100\t.\tA\tG\t.\tPASS\t.\tGT:AD:DP:PL\t0/1:5:10:50,0,50\t0/0:0:10:0,30,300",
        "chr1\t200\t.\tC\tT\t.\tPASS\t.\tGT:AD:DP:PL\t1/1:10:10:300,30,0\t1/1:.:.:.",
        "chr1\t300\t.\tG\tA\t.\tPASS\t.\tGT:AD:DP:PL\t0/1:5:10:50,0,50\t0/1:5:10:50,0,50"
    )
    writeLines(c(header_lines, col_line, data_lines), path)
}

test_that(".read_vireo_gt() classifies zygosity from the GT field per donor", {
    gt_file <- withr::local_tempfile(fileext = ".vcf")
    write_test_vireo_gt_vcf(gt_file)
    snp_info <- data.frame(snp_id = c("chr1:100:A:G", "chr1:200:C:T"))

    calls <- snplet:::.read_vireo_gt(gt_file, snp_info)

    snp1_d0 <- dplyr::filter(calls, snp_id == "chr1:100:A:G", donor == "donor0")
    snp1_d1 <- dplyr::filter(calls, snp_id == "chr1:100:A:G", donor == "donor1")
    # Verify a 0/1 call is classified heterozygous
    expect_equal(snp1_d0$zygosity, "het")
    # Verify a 0/0 call is classified homozygous
    expect_equal(snp1_d1$zygosity, "hom")
    # Verify the source is stamped for every call
    expect_true(all(calls$zygosity_source == "vireo_gt"))
})

test_that(".read_vireo_gt() derives zygosity_gt_prob as the normalised posterior of the called genotype", {
    gt_file <- withr::local_tempfile(fileext = ".vcf")
    write_test_vireo_gt_vcf(gt_file)
    snp_info <- data.frame(snp_id = c("chr1:100:A:G", "chr1:200:C:T"))

    calls <- snplet:::.read_vireo_gt(gt_file, snp_info)
    snp2_d0 <- dplyr::filter(calls, snp_id == "chr1:200:C:T", donor == "donor0")

    # PL = 300,30,0 for genotypes 0/0, 0/1, 1/1; the called 1/1 genotype's
    # posterior is 10^0 / (10^-30 + 10^-3 + 10^0)
    expect_equal(snp2_d0$zygosity_gt_prob, 0.999000999, tolerance = 1e-6)
})

test_that(".read_vireo_gt() leaves zygosity_gt_prob NA when the PL field is missing", {
    gt_file <- withr::local_tempfile(fileext = ".vcf")
    write_test_vireo_gt_vcf(gt_file)
    snp_info <- data.frame(snp_id = c("chr1:100:A:G", "chr1:200:C:T"))

    calls <- snplet:::.read_vireo_gt(gt_file, snp_info)
    snp2_d1 <- dplyr::filter(calls, snp_id == "chr1:200:C:T", donor == "donor1")

    # donor1's PL field for snp2 is "." (missing)
    expect_true(is.na(snp2_d1$zygosity_gt_prob))
})

test_that(".read_vireo_gt() restricts calls to SNPs present in snp_info", {
    gt_file <- withr::local_tempfile(fileext = ".vcf")
    write_test_vireo_gt_vcf(gt_file)
    snp_info <- data.frame(snp_id = c("chr1:100:A:G", "chr1:200:C:T"))

    calls <- snplet:::.read_vireo_gt(gt_file, snp_info)

    # chr1:300:G:A is in the VCF but absent from snp_info, so it must be dropped
    expect_false("chr1:300:G:A" %in% calls$snp_id)
})

test_that("import_cellsnp(vireo_folder=) populates donor_snp_info from real Vireo genotypes", {
    cellsnp_dir <- testthat::test_path("..", "..", "example_data", "village", "cell_snp", "cellsnp_HealthyVillage")
    vireo_folder <- testthat::test_path("..", "..", "example_data", "village", "vireo", "vireo_HealthyVillage")
    skip_if_not(
        all(file.exists(
            cellsnp_dir,
            file.path(vireo_folder, "donor_ids.tsv"),
            file.path(vireo_folder, "GT_donors.vireo.vcf.gz")
        )),
        "Village example dataset not found (not part of the installed package)"
    )

    gene_annotation <- data.frame(chrom = "chr1", start = 1, end = 1e9, gene_name = "dummy", strand = "+")
    snp_data <- import_cellsnp(
        cellsnp_dir = cellsnp_dir,
        gene_annotation = gene_annotation,
        vireo_folder = vireo_folder,
        library_id = "lib1",
    )

    donor_snp_info <- donor_snp_info(snp_data)
    # Verify zygosity calls were populated from the real GT VCF
    expect_true(nrow(donor_snp_info) > 0)
    # Verify every call is sourced from Vireo genotypes
    expect_true(all(donor_snp_info$zygosity_source == "vireo_gt"))
    # Verify every call's snp_id matches a SNP actually present in the object
    expect_true(all(donor_snp_info$snp_id %in% snp_info(snp_data)$snp_id))
    # Verify the object's active zygosity source was set automatically from
    # the imported Vireo calls
    expect_equal(zygosity_source(snp_data), "vireo_gt")
})

test_that("import_cellsnp() leaves the molecules slot empty", {
    snp_data <- get_example_snpdata()

    # Verify import stores no molecule data: the SNP-to-gene map and BAM paths
    # now belong to phase_from_molecules(), which is given them directly
    expect_equal(nrow(molecule_calls(molecules(snp_data))), 0L)
    expect_equal(nrow(snp_gene_map(molecules(snp_data))), 0L)
})

# ------------------------------------------------------------------------------
# read_mtx(): MatrixMarket format variations
# ------------------------------------------------------------------------------

write_test_mtx <- function(lines, fileext = ".mtx") {
    mtx_file <- withr::local_tempfile(fileext = fileext, .local_envir = parent.frame())
    writeLines(lines, mtx_file)
    mtx_file
}

mtx_banner <- "%%MatrixMarket matrix coordinate integer general"

# The 3 x 4 matrix every well-formed general fixture below encodes
expected_mtx <- Matrix::sparseMatrix(i = c(1, 2, 3), j = c(1, 3, 4), x = c(5, 7, 2), dims = c(3, 4))

test_that("read_mtx() reads a space-separated coordinate file", {
    mtx_file <- write_test_mtx(c(mtx_banner, "3 4 3", "1 1 5", "2 3 7", "3 4 2"))

    result <- snplet:::read_mtx(mtx_file)

    # Verify the parsed matrix matches the encoded entries and dimensions
    expect_equal(result, expected_mtx)
})

test_that("read_mtx() reads a tab-separated coordinate file", {
    mtx_file <- write_test_mtx(c(mtx_banner, "3\t4\t3", "1\t1\t5", "2\t3\t7", "3\t4\t2"))

    result <- snplet:::read_mtx(mtx_file)

    # Verify tab separators parse identically to single spaces
    expect_equal(result, expected_mtx)
})

test_that("read_mtx() tolerates mixed, repeated, leading and trailing whitespace", {
    mtx_file <- write_test_mtx(c(mtx_banner, "  3  4\t3 ", "1 \t1   5", "\t2 3 7\t", "3 4 2  "))

    result <- snplet:::read_mtx(mtx_file)

    # Confirm irregular whitespace does not produce NA indices
    expect_equal(result, expected_mtx)
})

test_that("read_mtx() handles CRLF line endings, blank lines and many comment lines", {
    comments <- rep("% a comment line", 100)
    mtx_file <- write_test_mtx(c(
        paste0(mtx_banner, "\r"),
        comments,
        "",
        "3 4 3\r",
        "1 1 5\r",
        "",
        "2 3 7\r",
        "3 4 2\r"
    ))

    result <- snplet:::read_mtx(mtx_file)

    # Verify the size line is found beyond the initial header read and blank lines are skipped
    expect_equal(result, expected_mtx)
})

test_that("read_mtx() parses indices and values written in scientific notation", {
    mtx_file <- write_test_mtx(c(mtx_banner, "3e+00 4 3", "1 1 5", "2 3 7e+00", "3e0 4 2"))

    result <- snplet:::read_mtx(mtx_file)

    # Confirm scientific notation is read as whole numbers
    expect_equal(result, expected_mtx)
})

test_that("read_mtx() reads a gzip-compressed file", {
    mtx_file <- withr::local_tempfile(fileext = ".mtx.gz")
    gz_con <- gzfile(mtx_file, "w")
    writeLines(c(mtx_banner, "3 4 3", "1 1 5", "2 3 7", "3 4 2"), gz_con)
    close(gz_con)

    result <- snplet:::read_mtx(mtx_file)

    # Verify compressed input is decompressed transparently
    expect_equal(result, expected_mtx)
})

test_that("read_mtx() reads a file without a banner as coordinate real general", {
    mtx_file <- write_test_mtx(c("3 4 3", "1 1 5", "2 3 7", "3 4 2"))

    result <- snplet:::read_mtx(mtx_file)

    # Confirm the missing banner falls back to the general coordinate layout
    expect_equal(result, expected_mtx)
})

test_that("read_mtx() assigns a value of one to every entry of a pattern file", {
    mtx_file <- write_test_mtx(c(
        "%%MatrixMarket matrix coordinate pattern general",
        "3 4 3",
        "1 1",
        "2 3",
        "3 4"
    ))

    result <- snplet:::read_mtx(mtx_file)

    # Verify pattern entries become ones at the listed positions
    expect_equal(result, Matrix::sparseMatrix(i = c(1, 2, 3), j = c(1, 3, 4), x = 1, dims = c(3, 4)))
})

test_that("read_mtx() mirrors the stored triangle of a symmetric file", {
    mtx_file <- write_test_mtx(c(
        "%%MatrixMarket matrix coordinate real symmetric",
        "3 3 3",
        "1 1 4",
        "2 1 6",
        "3 2 8"
    ))

    result <- snplet:::read_mtx(mtx_file)

    # Verify off-diagonal entries appear in both triangles and the diagonal is not doubled
    expect_equal(as.matrix(result), matrix(c(4, 6, 0, 6, 0, 8, 0, 8, 0), nrow = 3))
})

test_that("read_mtx() negates the mirrored triangle of a skew-symmetric file", {
    mtx_file <- write_test_mtx(c(
        "%%MatrixMarket matrix coordinate real skew-symmetric",
        "2 2 1",
        "2 1 3"
    ))

    result <- snplet:::read_mtx(mtx_file)

    # Verify the upper triangle holds the negated lower-triangle entry
    expect_equal(as.matrix(result), matrix(c(0, 3, -3, 0), nrow = 2))
})

test_that("read_mtx() converts a dense array file to a sparse general matrix", {
    mtx_file <- write_test_mtx(c("%%MatrixMarket matrix array real general", "2 2", "1", "0", "0", "4"))

    result <- snplet:::read_mtx(mtx_file)

    # Verify array input is returned as a dgCMatrix
    expect_s4_class(result, "dgCMatrix")
    # Check the column-major values are placed correctly
    expect_equal(as.matrix(result), matrix(c(1, 0, 0, 4), nrow = 2))
})

test_that("read_mtx() errors when the file holds fewer entries than declared", {
    mtx_file <- write_test_mtx(c(mtx_banner, "3 4 3", "1 1 5", "2 3 7"))

    # Ensure a truncated file is reported rather than silently accepted
    expect_error(snplet:::read_mtx(mtx_file), "declares 3 entries but contains 2")
})

test_that("read_mtx() errors on a line with a missing field", {
    mtx_file <- write_test_mtx(c(mtx_banner, "3 4 3", "1 1 5", "2 3", "3 4 2"))

    # Ensure the malformed entry is located in the error message
    expect_error(snplet:::read_mtx(mtx_file), "malformed entry at entry 2")
})

test_that("read_mtx() errors on a line with surplus fields", {
    mtx_file <- write_test_mtx(c(mtx_banner, "3 4 3", "1 1 5", "2 3 7 9", "3 4 2"))

    # Ensure surplus fields are not silently dropped
    expect_error(snplet:::read_mtx(mtx_file), "malformed entry at entry 2")
})

test_that("read_mtx() errors on an index outside the declared dimensions", {
    mtx_file <- write_test_mtx(c(mtx_banner, "3 4 3", "1 1 5", "2 5 7", "3 4 2"))

    # Ensure out-of-range column indices are rejected
    expect_error(snplet:::read_mtx(mtx_file), "invalid index at entry 2")
})

test_that("read_mtx() errors on a complex matrix", {
    mtx_file <- write_test_mtx(c("%%MatrixMarket matrix coordinate complex general", "1 1 1", "1 1 2 3"))

    # Ensure unsupported fields are rejected with a clear message
    expect_error(snplet:::read_mtx(mtx_file), "unsupported field 'complex'")
})

test_that("read_mtx() errors on a malformed size line", {
    mtx_file <- write_test_mtx(c(mtx_banner, "3 4", "1 1 5"))

    # Ensure a coordinate size line missing nnz is reported
    expect_error(snplet:::read_mtx(mtx_file), "malformed size line")
})

# ==============================================================================
# import_cellsnp_libraries()
# ==============================================================================

example_cellsnp_dir <- function() {
    system.file("extdata/example_snpdata", package = "snplet")
}

example_gene_annotation <- function() {
    readr::read_tsv(system.file("extdata/example_gene_anno.tsv", package = "snplet"), show_col_types = FALSE)
}

# A second directory holding the same cellSNP-lite output, standing in for a
# repeat run of one library (a sheet may list each directory only once).
copy_example_cellsnp_dir <- function(env = parent.frame()) {
    copy_dir <- withr::local_tempdir(.local_envir = env)
    file.copy(list.files(example_cellsnp_dir(), full.names = TRUE), copy_dir)
    copy_dir
}

# Two cellSNP-lite directories splitting the example SNPs in half, standing in
# for runs over one BAM that pile up disjoint SNPs (e.g. split by chromosome).
export_example_snp_halves <- function(env = parent.frame()) {
    single <- import_cellsnp(example_cellsnp_dir(), example_gene_annotation(), library_id = "lib1")
    first_half <- seq_len(nrow(single)) <= nrow(single) / 2
    split_dirs <- c(withr::local_tempdir(.local_envir = env), withr::local_tempdir(.local_envir = env))
    export_cellsnp(single[first_half, ], split_dirs[1])
    export_cellsnp(single[!first_half, ], split_dirs[2])
    split_dirs
}

# Empty files standing in for BAMs: import only checks that the paths exist.
local_empty_bams <- function(file_names, env = parent.frame()) {
    bam_dir <- withr::local_tempdir(.local_envir = env)
    bam_paths <- file.path(bam_dir, file_names)
    file.create(bam_paths)
    bam_paths
}

test_that("import_cellsnp_libraries() with one row matches import_cellsnp()", {
    sheet <- tibble::tibble(cellsnp_dir = example_cellsnp_dir(), library_id = "lib1")

    from_sheet <- import_cellsnp_libraries(sheet, example_gene_annotation())
    direct <- import_cellsnp(example_cellsnp_dir(), example_gene_annotation(), library_id = "lib1")

    # Verify a one-run sheet yields the same counts as a direct import
    expect_equal(ref_count(from_sheet), ref_count(direct))
    # Verify the same cells, with the same metadata, are imported
    expect_equal(barcode_info(from_sheet), barcode_info(direct))
})

test_that("import_cellsnp_libraries() keeps cells from different libraries separate", {
    sheet <- tibble::tibble(
        cellsnp_dir = c(example_cellsnp_dir(), copy_example_cellsnp_dir()),
        library_id = c("lib1", "lib2"),
        donor_map = list(c(donor_lib1 = "donor0"), c(donor_lib2 = "donor0"))
    )
    single <- import_cellsnp(example_cellsnp_dir(), example_gene_annotation(), library_id = "lib1")

    combined <- import_cellsnp_libraries(sheet, example_gene_annotation())

    # Verify barcodes shared by chance across libraries stay as two cells each
    expect_equal(ncol(combined), 2 * ncol(single))
    # Confirm every cell carries the library label of its own row
    expect_equal(as.vector(table(barcode_info(combined)$library_id)), c(ncol(single), ncol(single)))
    # Verify each row's donor_map was applied before combining
    expect_setequal(donor_info(combined)$donor, c("donor_lib1", "donor_lib2"))
})

test_that("import_cellsnp_libraries() sums same-BAM runs that cover disjoint SNPs", {
    split_dirs <- export_example_snp_halves()
    bam_file <- local_empty_bams("lib1.bam")
    sheet <- tibble::tibble(cellsnp_dir = split_dirs, library_id = "lib1", bam_files = bam_file)
    single <- import_cellsnp(example_cellsnp_dir(), example_gene_annotation(), library_id = "lib1")

    combined <- expect_no_warning(import_cellsnp_libraries(sheet, example_gene_annotation()))

    # Verify the split runs rejoin into the original cells
    expect_equal(ncol(combined), ncol(single))
    # Confirm each read is counted once when the runs cover different SNPs
    expect_equal(sum(alt_count(combined)), sum(alt_count(single)))
})

test_that("import_cellsnp_libraries() errors when same-BAM runs share SNPs", {
    bam_file <- local_empty_bams("lib1.bam")
    sheet <- tibble::tibble(
        cellsnp_dir = c(example_cellsnp_dir(), copy_example_cellsnp_dir()),
        library_id = "lib1",
        bam_files = bam_file
    )

    # Ensure reads from one BAM cannot be counted twice at shared SNPs
    expect_error(
        import_cellsnp_libraries(sheet, example_gene_annotation()),
        "share a BAM file and SNPs"
    )
})

test_that("import_cellsnp_libraries() sums shared SNPs of runs over different BAMs", {
    bam_files <- local_empty_bams(c("lib1_seq1.bam", "lib1_seq2.bam"))
    sheet <- tibble::tibble(
        cellsnp_dir = c(example_cellsnp_dir(), copy_example_cellsnp_dir()),
        library_id = "lib1",
        bam_files = bam_files
    )
    single <- import_cellsnp(example_cellsnp_dir(), example_gene_annotation(), library_id = "lib1")

    combined <- expect_no_warning(import_cellsnp_libraries(sheet, example_gene_annotation()))

    # Verify a library sequenced twice contributes no new cells
    expect_equal(ncol(combined), ncol(single))
    # Confirm the second sequencing's reads are added to the shared cells
    expect_equal(sum(alt_count(combined)), 2 * sum(alt_count(single)))
})

test_that("import_cellsnp_libraries() warns when same-library runs share SNPs without BAM paths", {
    sheet <- tibble::tibble(
        cellsnp_dir = c(example_cellsnp_dir(), copy_example_cellsnp_dir()),
        library_id = "lib1"
    )

    # Ensure an unverifiable risk of double counting is reported
    expect_warning(
        import_cellsnp_libraries(sheet, example_gene_annotation()),
        "Runs of the same library share SNPs"
    )
})

test_that("import_cellsnp_libraries() errors when runs of a library disagree on a cell's donor", {
    bam_files <- local_empty_bams(c("lib1_seq1.bam", "lib1_seq2.bam"))
    sheet <- tibble::tibble(
        cellsnp_dir = c(example_cellsnp_dir(), copy_example_cellsnp_dir()),
        library_id = "lib1",
        bam_files = bam_files,
        donor_map = list(NULL, c(other_donor = "donor0"))
    )

    # Ensure a cell is not silently kept under the first run's donor label
    expect_error(
        import_cellsnp_libraries(sheet, example_gene_annotation()),
        "assigned different donors by different runs of the same library"
    )
})

test_that("import_cellsnp_libraries() does not store the sheet's BAM paths on the object", {
    bam_dir <- withr::local_tempdir()
    bam_paths <- file.path(bam_dir, c("lib1.bam", "lib2.bam"))
    file.create(bam_paths)
    sheet <- tibble::tibble(
        cellsnp_dir = c(example_cellsnp_dir(), copy_example_cellsnp_dir()),
        library_id = c("lib1", "lib2"),
        bam_files = bam_paths,
        donor_map = list(c(donor_lib1 = "donor0"), c(donor_lib2 = "donor0"))
    )

    combined <- import_cellsnp_libraries(sheet, example_gene_annotation())

    # Verify the paths serve only the repeat-run check: library_info records
    # libraries, not files
    expect_named(library_info(combined), c("library_id", "n_cells"))
})

test_that("import_cellsnp_libraries() errors on a BAM path that does not exist", {
    sheet <- tibble::tibble(
        cellsnp_dir = example_cellsnp_dir(),
        library_id = "lib1",
        bam_files = file.path(withr::local_tempdir(), "missing.bam")
    )

    # Ensure a mistyped path is caught before import rather than weakening the
    # repeat-run check
    expect_error(import_cellsnp_libraries(sheet, example_gene_annotation()), "Required file not found")
})

test_that("import_cellsnp_libraries() treats NA optional entries as absent", {
    sheet <- tibble::tibble(
        cellsnp_dir = example_cellsnp_dir(),
        library_id = "lib1",
        vdj_file = NA_character_,
        vireo_folder = NA_character_
    )

    combined <- import_cellsnp_libraries(sheet, example_gene_annotation())

    # Confirm an NA vdj_file imports no clonotypes rather than erroring
    expect_true(all(is.na(barcode_info(combined)$clonotype)))
})

test_that("import_cellsnp_libraries() errors when a donor label spans libraries", {
    sheet <- tibble::tibble(
        cellsnp_dir = c(example_cellsnp_dir(), copy_example_cellsnp_dir()),
        library_id = c("lib1", "lib2")
    )

    # Ensure the default donor0 label in two libraries is refused, not merged
    expect_error(import_cellsnp_libraries(sheet, example_gene_annotation()), "more than one library: donor0")
})

test_that("import_cellsnp_libraries() rejects a sheet without library_id", {
    sheet <- tibble::tibble(cellsnp_dir = example_cellsnp_dir())

    # Ensure the library label, which cell matching depends on, is required
    expect_error(import_cellsnp_libraries(sheet, example_gene_annotation()), "missing required column\\(s\\): library_id")
})

test_that("import_cellsnp_libraries() rejects an NA library_id", {
    sheet <- tibble::tibble(cellsnp_dir = example_cellsnp_dir(), library_id = NA_character_)

    # Ensure an unlabelled row cannot be combined
    expect_error(import_cellsnp_libraries(sheet, example_gene_annotation()), "Every row of sheet needs a library_id")
})

test_that("import_cellsnp_libraries() rejects an unrecognised column", {
    sheet <- tibble::tibble(cellsnp_dir = example_cellsnp_dir(), library_id = "lib1", vireo_dir = "vireo/")

    # Ensure a misspelt optional column is reported rather than silently ignored
    expect_error(import_cellsnp_libraries(sheet, example_gene_annotation()), "unrecognised column\\(s\\): vireo_dir")
})

test_that("import_cellsnp_libraries() rejects a directory listed twice", {
    sheet <- tibble::tibble(cellsnp_dir = rep(example_cellsnp_dir(), 2), library_id = "lib1")

    # Ensure one run cannot be imported twice and have its counts doubled
    expect_error(import_cellsnp_libraries(sheet, example_gene_annotation()), "same cellsnp_dir more than once")
})

test_that("import_cellsnp_libraries() rejects an empty sheet", {
    sheet <- tibble::tibble(cellsnp_dir = character(0), library_id = character(0))

    # Ensure a sheet with no runs is refused
    expect_error(import_cellsnp_libraries(sheet, example_gene_annotation()), "one row per cellSNP-lite run")
})

# ==============================================================================
# import_cellsnp() cell barcode alignment
# ==============================================================================

read_example_donor_ids <- function(cellsnp_dir) {
    readr::read_tsv(file.path(cellsnp_dir, "donor_ids.tsv"), col_types = readr::cols(.default = "c"))
}

write_example_donor_ids <- function(donor_ids, cellsnp_dir) {
    readr::write_tsv(donor_ids, file.path(cellsnp_dir, "donor_ids.tsv"))
}

test_that("import_cellsnp() matches Vireo donors to cells by barcode, not row order", {
    cellsnp_dir <- copy_example_cellsnp_dir()
    donor_ids <- read_example_donor_ids(cellsnp_dir)
    write_example_donor_ids(donor_ids[rev(seq_len(nrow(donor_ids))), ], cellsnp_dir)

    snp_data <- import_cellsnp(cellsnp_dir, example_gene_annotation(), vireo_folder = cellsnp_dir)

    imported <- barcode_info(snp_data)
    # Verify cells stay in cellSNP.samples.tsv order despite the reversed Vireo file
    expect_equal(imported$barcode, donor_ids$barcode)
    # Confirm each cell keeps the donor Vireo assigned to its own barcode
    expect_equal(imported$donor, donor_ids$donor_id)
})

test_that("import_cellsnp() errors when donor_ids.tsv lacks a cellSNP barcode", {
    cellsnp_dir <- copy_example_cellsnp_dir()
    donor_ids <- read_example_donor_ids(cellsnp_dir)
    write_example_donor_ids(donor_ids[-1, ], cellsnp_dir)

    # Ensure a Vireo file missing a cell is refused rather than shifted onto other cells
    expect_error(
        import_cellsnp(cellsnp_dir, example_gene_annotation(), vireo_folder = cellsnp_dir),
        "does not match the cellSNP-lite barcodes"
    )
})

test_that("import_cellsnp() errors when donor_ids.tsv names a barcode cellSNP did not count", {
    cellsnp_dir <- copy_example_cellsnp_dir()
    donor_ids <- read_example_donor_ids(cellsnp_dir)
    donor_ids$barcode[1] <- "NOTABARCODE"
    write_example_donor_ids(donor_ids, cellsnp_dir)

    # Ensure a Vireo file from a different cellSNP run is refused
    expect_error(
        import_cellsnp(cellsnp_dir, example_gene_annotation(), vireo_folder = cellsnp_dir),
        "NOTABARCODE"
    )
})

test_that("import_cellsnp() errors when donor_ids.tsv repeats a barcode", {
    cellsnp_dir <- copy_example_cellsnp_dir()
    donor_ids <- read_example_donor_ids(cellsnp_dir)
    write_example_donor_ids(rbind(donor_ids, donor_ids[1, ]), cellsnp_dir)

    # Ensure a cell with two Vireo rows is refused rather than duplicated
    expect_error(
        import_cellsnp(cellsnp_dir, example_gene_annotation(), vireo_folder = cellsnp_dir),
        "more than one row"
    )
})

test_that("import_cellsnp() errors when cellSNP.samples.tsv is missing", {
    cellsnp_dir <- copy_example_cellsnp_dir()
    file.remove(file.path(cellsnp_dir, "cellSNP.samples.tsv"))

    # Ensure the missing barcode file is reported by name
    expect_error(
        import_cellsnp(cellsnp_dir, example_gene_annotation()),
        "Required file not found: .*cellSNP.samples.tsv"
    )
})

test_that("import_cellsnp() errors when cellSNP.samples.tsv disagrees with the matrix column count", {
    cellsnp_dir <- copy_example_cellsnp_dir()
    samples_file <- file.path(cellsnp_dir, "cellSNP.samples.tsv")
    writeLines(readLines(samples_file)[-1], samples_file)

    # Ensure a barcode list that cannot label every matrix column is refused
    expect_error(
        import_cellsnp(cellsnp_dir, example_gene_annotation()),
        "55 barcodes but the count matrices have 56 columns"
    )
})

test_that("import_cellsnp() errors when cellSNP.samples.tsv repeats a barcode", {
    cellsnp_dir <- copy_example_cellsnp_dir()
    samples_file <- file.path(cellsnp_dir, "cellSNP.samples.tsv")
    barcodes <- readLines(samples_file)
    writeLines(c(barcodes[-1], barcodes[2]), samples_file)

    # Ensure repeated barcodes cannot silently give two columns the same cell
    expect_error(import_cellsnp(cellsnp_dir, example_gene_annotation()), "missing or repeated barcodes")
})

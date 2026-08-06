# ==============================================================================
# Test Suite: VCF/BAM External Tool Utilities
# Description: Tests for availability checks and system-tool wrappers used to
#              extract exon BEDs, filter VCFs, and compute BAM depth
# ==============================================================================

library(testthat)

# ==============================================================================
# Test: tool availability checks
# ==============================================================================

test_that("check_awk_available() errors when awk is not on PATH", {
    local_mocked_bindings(Sys.which = function(...) c(awk = ""), .package = "base")

    # Verify error is raised when awk cannot be found
    expect_error(check_awk_available(), "awk is not available")
})

test_that("check_awk_available() passes when awk is on PATH", {
    local_mocked_bindings(Sys.which = function(...) c(awk = "/usr/bin/awk"), .package = "base")

    # Verify no error is raised when awk is found
    expect_no_error(check_awk_available())
})

test_that("check_bcftools_available() errors when bcftools is not on PATH", {
    local_mocked_bindings(Sys.which = function(...) c(bcftools = ""), .package = "base")

    # Verify error is raised when bcftools cannot be found
    expect_error(check_bcftools_available(), "bcftools is not available")
})

test_that("check_bcftools_available() passes when bcftools is on PATH", {
    local_mocked_bindings(Sys.which = function(...) c(bcftools = "/usr/bin/bcftools"), .package = "base")

    # Verify no error is raised when bcftools is found
    expect_no_error(check_bcftools_available())
})

test_that("check_samtools_available() errors when samtools is not on PATH", {
    local_mocked_bindings(Sys.which = function(...) c(samtools = ""), .package = "base")

    # Verify error is raised when samtools cannot be found
    expect_error(check_samtools_available(), "samtools is not available")
})

test_that("check_samtools_available() passes when samtools is on PATH", {
    local_mocked_bindings(Sys.which = function(...) c(samtools = "/usr/bin/samtools"), .package = "base")

    # Verify no error is raised when samtools is found
    expect_no_error(check_samtools_available())
})

test_that("check_gunzip_available() errors when gunzip is not on PATH", {
    local_mocked_bindings(Sys.which = function(...) c(gunzip = ""), .package = "base")

    # Verify error is raised when gunzip cannot be found
    expect_error(check_gunzip_available(), "gunzip is not available")
})

test_that("check_gunzip_available() passes when gunzip is on PATH", {
    local_mocked_bindings(Sys.which = function(...) c(gunzip = "/usr/bin/gunzip"), .package = "base")

    # Verify no error is raised when gunzip is found
    expect_no_error(check_gunzip_available())
})

# ==============================================================================
# Test: gff_exons_to_bed
# ==============================================================================

test_that("gff_exons_to_bed() errors when awk is unavailable", {
    local_mocked_bindings(check_awk_available = function() stop("awk is not available."))

    # Verify error propagates from the awk availability check
    expect_error(gff_exons_to_bed("nonexistent.gff", withr::local_tempfile()), "awk is not available")
})

test_that("gff_exons_to_bed() errors when the GFF file does not exist", {
    skip_if_not(nzchar(Sys.which("awk")), "awk not available")

    # Verify error is raised for a missing input file (normalizePath's own
    # mustWork = TRUE check fires before the function's explicit file check)
    expect_error(
        gff_exons_to_bed("does_not_exist.gff", withr::local_tempfile()),
        "No such file or directory"
    )
})

test_that("gff_exons_to_bed() extracts exon coordinates to a BED file", {
    skip_if_not(nzchar(Sys.which("awk")), "awk not available")

    gff_file <- withr::local_tempfile(fileext = ".gff")
    writeLines(
        c(
            "chr1\tsource\tgene\t1\t1000\t.\t+\t.\tID=gene1",
            "chr1\tsource\texon\t100\t200\t.\t+\t.\tID=exon1",
            "chr1\tsource\texon\t300\t400\t.\t+\t.\tID=exon2"
        ),
        gff_file
    )
    output_bed <- withr::local_tempfile(fileext = ".bed")

    gff_exons_to_bed(gff_file, output_bed)

    # Verify the BED file was created
    expect_true(file.exists(output_bed))
    bed_lines <- readLines(output_bed)
    # Verify only exon rows were extracted, with 0-based start coordinates
    expect_equal(bed_lines, c("chr1\t99\t200", "chr1\t299\t400"))
})

# ==============================================================================
# Test: filter_vcf_by_bed
# ==============================================================================

test_that("filter_vcf_by_bed() errors when bcftools is unavailable", {
    local_mocked_bindings(check_bcftools_available = function() stop("bcftools is not available."))

    # Verify error propagates from the bcftools availability check
    expect_error(
        filter_vcf_by_bed("nonexistent.vcf", "nonexistent.bed", withr::local_tempfile()),
        "bcftools is not available"
    )
})

test_that("filter_vcf_by_bed() errors when the input VCF does not exist", {
    skip_if_not(nzchar(Sys.which("bcftools")), "bcftools not available")

    bed_file <- withr::local_tempfile(fileext = ".bed")
    writeLines("chr1\t0\t100", bed_file)

    # Verify error is raised for a missing input VCF (normalizePath's own
    # mustWork = TRUE check fires before the function's explicit file check)
    expect_error(
        filter_vcf_by_bed("does_not_exist.vcf", bed_file, withr::local_tempfile()),
        "No such file or directory"
    )
})

test_that("filter_vcf_by_bed() errors when the BED file does not exist", {
    skip_if_not(nzchar(Sys.which("bcftools")), "bcftools not available")

    vcf_file <- withr::local_tempfile(fileext = ".vcf")
    writeLines(c("##fileformat=VCFv4.2", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"), vcf_file)

    # Verify error is raised for a missing BED file (normalizePath's own
    # mustWork = TRUE check fires before the function's explicit file check)
    expect_error(
        filter_vcf_by_bed(vcf_file, "does_not_exist.bed", withr::local_tempfile()),
        "No such file or directory"
    )
})

test_that("filter_vcf_by_bed() filters a bgzipped, indexed VCF by BED regions", {
    tools_available <- nzchar(Sys.which("bcftools")) &&
        nzchar(Sys.which("bgzip")) &&
        nzchar(Sys.which("tabix"))
    skip_if_not(tools_available, "bcftools/bgzip/tabix not available")

    vcf_dir <- withr::local_tempdir()
    vcf_file <- file.path(vcf_dir, "test.vcf")
    writeLines(
        c(
            "##fileformat=VCFv4.2",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
            "chr1\t100\t.\tA\tT\t.\tPASS\t.",
            "chr1\t500\t.\tG\tC\t.\tPASS\t.",
            "chr2\t100\t.\tA\tT\t.\tPASS\t."
        ),
        vcf_file
    )
    system2(Sys.which("bgzip"), args = c("-f", vcf_file))
    vcf_gz_file <- paste0(vcf_file, ".gz")
    system2(Sys.which("tabix"), args = c("-p", "vcf", vcf_gz_file))

    bed_file <- file.path(vcf_dir, "regions.bed")
    writeLines("chr1\t0\t200", bed_file)
    output_vcf <- file.path(vcf_dir, "out.vcf")

    filter_vcf_by_bed(vcf_gz_file, bed_file, output_vcf)

    # Verify the filtered output file was created
    expect_true(file.exists(output_vcf))
    out_lines <- readLines(output_vcf)
    data_lines <- out_lines[!grepl("^#", out_lines)]
    # Verify only the variant within the BED region is retained
    expect_equal(data_lines, "chr1\t100\t.\tA\tT\t.\tPASS\t.")
})

# ==============================================================================
# Test: bam_depth_at_positions
# ==============================================================================

test_that("bam_depth_at_positions() errors when samtools is unavailable", {
    local_mocked_bindings(check_samtools_available = function() stop("samtools is not available."))

    # Verify error propagates from the samtools availability check
    expect_error(
        bam_depth_at_positions("nonexistent.bam", "nonexistent.vcf", withr::local_tempfile()),
        "samtools is not available"
    )
})

test_that("bam_depth_at_positions() errors when the BAM file does not exist", {
    skip_if_not(
        nzchar(Sys.which("samtools")) && nzchar(Sys.which("awk")) && nzchar(Sys.which("gunzip")),
        "samtools/awk/gunzip not available"
    )

    vcf_file <- withr::local_tempfile(fileext = ".vcf")
    writeLines(c("##fileformat=VCFv4.2", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"), vcf_file)

    # Verify error is raised for a missing BAM file (normalizePath's own
    # mustWork = TRUE check fires before the function's explicit file check)
    expect_error(
        bam_depth_at_positions("does_not_exist.bam", vcf_file, withr::local_tempfile()),
        "No such file or directory"
    )
})

test_that("bam_depth_at_positions() errors when the VCF file does not exist", {
    skip_if_not(
        nzchar(Sys.which("samtools")) && nzchar(Sys.which("awk")) && nzchar(Sys.which("gunzip")),
        "samtools/awk/gunzip not available"
    )

    bam_file <- withr::local_tempfile(fileext = ".bam")
    writeLines("not a real bam", bam_file)

    # Verify error is raised for a missing VCF file (normalizePath's own
    # mustWork = TRUE check fires before the function's explicit file check)
    expect_error(
        bam_depth_at_positions(bam_file, "does_not_exist.vcf", withr::local_tempfile()),
        "No such file or directory"
    )
})

test_that("bam_depth_at_positions() computes read depth at VCF positions from a BAM file", {
    tools_available <- nzchar(Sys.which("samtools")) &&
        nzchar(Sys.which("awk")) &&
        nzchar(Sys.which("gunzip"))
    skip_if_not(tools_available, "samtools/awk/gunzip not available")

    bam_dir <- withr::local_tempdir()
    sam_file <- file.path(bam_dir, "test.sam")
    writeLines(
        c(
            "@HD\tVN:1.6\tSO:coordinate",
            "@SQ\tSN:chr1\tLN:1000",
            "read1\t0\tchr1\t95\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII",
            "read2\t0\tchr1\t100\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII"
        ),
        sam_file
    )
    bam_file <- file.path(bam_dir, "test.bam")
    system2(Sys.which("samtools"), args = c("view", "-b", "-o", bam_file, sam_file))
    system2(Sys.which("samtools"), args = c("index", bam_file))

    vcf_file <- file.path(bam_dir, "positions.vcf")
    writeLines(
        c(
            "##fileformat=VCFv4.2",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
            "chr1\t100\t.\tA\tT\t.\tPASS\t."
        ),
        vcf_file
    )
    output_file <- file.path(bam_dir, "depth.txt")

    bam_depth_at_positions(bam_file, vcf_file, output_file)

    # Verify the depth output file was created
    expect_true(file.exists(output_file))
    depth_lines <- readLines(output_file)
    # Verify depth is reported at the VCF position, covered by both reads
    expect_equal(depth_lines, "chr1\t100\t2")
})

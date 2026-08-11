# ==============================================================================
# Test Suite: Molecule-level allele counting and read-backed phasing
# Description: Tests for extract_snp_calls(), molecule_snp_alleles(),
#              phase_snps(), and assign_snp_genes()
# ==============================================================================

library(testthat)
library(dplyr)

# ------------------------------------------------------------------------------
# Test Data Setup
# ------------------------------------------------------------------------------

# A small, real chrX BAM fixture (2603 alignments overlap chrX:30169685, all
# with CB/UB tags and mapq >= 20) used for the extract_snp_calls() integration
# tests below.
test_bam <- system.file("../../example_data/snp_counting/TASL.bam", package = "snplet")
if (!nzchar(test_bam) || !file.exists(test_bam)) {
    test_bam <- "../../example_data/snp_counting/TASL.bam"
}

# ==============================================================================
# Test: molecule_snp_alleles()
# ==============================================================================

test_that("molecule_snp_alleles() takes the majority call per (molecule, SNP)", {
    tallies <- tibble::tribble(
        ~barcode,
        ~umi,
        ~snp_id,
        ~allele,
        ~n_calls,
        "AAA",
        "UMI1",
        "snp1",
        "REF",
        3L,
        "AAA",
        "UMI1",
        "snp1",
        "ALT",
        1L,
        "AAA",
        "UMI1",
        "snp2",
        "OTH",
        2L
    )

    result <- molecule_snp_alleles(tallies)

    # Verify only the majority allele survives for the contested (molecule, SNP)
    expect_equal(result$allele[result$snp_id == "snp1"], "REF")
    # Check that a molecule/SNP call resolving to OTH is dropped entirely
    expect_false("snp2" %in% result$snp_id)
})

test_that("molecule_snp_alleles() drops OTH-only calls and keeps ties resolved deterministically", {
    tallies <- tibble::tribble(
        ~barcode,
        ~umi,
        ~snp_id,
        ~allele,
        ~n_calls,
        "BBB",
        "UMI2",
        "snp1",
        "REF",
        2L,
        "BBB",
        "UMI2",
        "snp1",
        "ALT",
        2L
    )

    result <- molecule_snp_alleles(tallies)

    # Verify a tie still resolves to exactly one row (with_ties = FALSE)
    expect_equal(nrow(result), 1)
    # Check the resolved allele is one of the tied candidates
    expect_true(result$allele %in% c("REF", "ALT"))
})

# ==============================================================================
# Test: molecule_read_strand()
# ==============================================================================

test_that("molecule_read_strand() takes the majority alignment strand per molecule", {
    reads <- tibble::tribble(
        ~barcode,
        ~umi,
        ~qname,
        ~strand,
        "AAA",
        "UMI1",
        "read1",
        "+",
        "AAA",
        "UMI1",
        "read2",
        "+",
        "AAA",
        "UMI1",
        "read3",
        "-"
    )

    result <- molecule_read_strand(reads)

    # Verify the majority strand ("+", 2 of 3 reads) wins for the molecule
    expect_equal(result$strand[result$barcode == "AAA" & result$umi == "UMI1"], "+")
    expect_equal(nrow(result), 1)
})

test_that("molecule_read_strand() resolves one row per distinct molecule", {
    reads <- tibble::tribble(
        ~barcode,
        ~umi,
        ~qname,
        ~strand,
        "AAA",
        "UMI1",
        "read1",
        "+",
        "BBB",
        "UMI2",
        "read2",
        "-"
    )

    result <- molecule_read_strand(reads)

    # Verify each distinct molecule gets its own row and strand
    expect_equal(nrow(result), 2)
    expect_equal(result$strand[result$barcode == "AAA"], "+")
    expect_equal(result$strand[result$barcode == "BBB"], "-")
})

# ==============================================================================
# Test: phase_snps()
# ==============================================================================

test_that("phase_snps() links two SNPs observed on the same molecules as 'same'", {
    per_snp <- tibble::tribble(
        ~barcode,
        ~umi,
        ~snp_id,
        ~allele,
        "c1",
        "u1",
        "snp_a",
        "REF",
        "c1",
        "u1",
        "snp_b",
        "REF",
        "c2",
        "u2",
        "snp_a",
        "REF",
        "c2",
        "u2",
        "snp_b",
        "REF",
        "c3",
        "u3",
        "snp_a",
        "ALT",
        "c3",
        "u3",
        "snp_b",
        "ALT",
        "c4",
        "u4",
        "snp_a",
        "REF",
        "c4",
        "u4",
        "snp_b",
        "REF",
        "c5",
        "u5",
        "snp_a",
        "ALT",
        "c5",
        "u5",
        "snp_b",
        "ALT"
    )

    result <- phase_snps(per_snp, min_molecules = 5, min_consistency = 0.9)

    # Verify both SNPs land in the same phase block
    expect_equal(dplyr::n_distinct(result$block), 1)
    # Check that co-occurring REF alleles are assigned the same H1 label
    h1 <- setNames(result$allele_on_h1, result$snp_id)
    expect_equal(h1[["snp_a"]], h1[["snp_b"]])
})

test_that("phase_snps() flips orientation for SNPs observed as 'opposite'", {
    per_snp <- tibble::tribble(
        ~barcode,
        ~umi,
        ~snp_id,
        ~allele,
        "c1",
        "u1",
        "snp_a",
        "REF",
        "c1",
        "u1",
        "snp_b",
        "ALT",
        "c2",
        "u2",
        "snp_a",
        "REF",
        "c2",
        "u2",
        "snp_b",
        "ALT",
        "c3",
        "u3",
        "snp_a",
        "ALT",
        "c3",
        "u3",
        "snp_b",
        "REF",
        "c4",
        "u4",
        "snp_a",
        "REF",
        "c4",
        "u4",
        "snp_b",
        "ALT",
        "c5",
        "u5",
        "snp_a",
        "ALT",
        "c5",
        "u5",
        "snp_b",
        "REF"
    )

    result <- phase_snps(per_snp, min_molecules = 5, min_consistency = 0.9)

    # Verify opposite-travelling alleles get different H1 labels
    h1 <- setNames(result$allele_on_h1, result$snp_id)
    expect_false(h1[["snp_a"]] == h1[["snp_b"]])
})

test_that("phase_snps() rejects an edge below min_molecules", {
    per_snp <- tibble::tribble(
        ~barcode,
        ~umi,
        ~snp_id,
        ~allele,
        "c1",
        "u1",
        "snp_a",
        "REF",
        "c1",
        "u1",
        "snp_b",
        "REF",
        "c2",
        "u2",
        "snp_a",
        "REF",
        "c2",
        "u2",
        "snp_b",
        "REF"
    )

    result <- phase_snps(per_snp, min_molecules = 5, min_consistency = 0.9)

    # Verify no edge is accepted when supporting molecules fall short of min_molecules
    expect_equal(nrow(result), 0)
})

test_that("phase_snps() rejects an edge below min_consistency", {
    per_snp <- tibble::tribble(
        ~barcode,
        ~umi,
        ~snp_id,
        ~allele,
        "c1",
        "u1",
        "snp_a",
        "REF",
        "c1",
        "u1",
        "snp_b",
        "REF",
        "c2",
        "u2",
        "snp_a",
        "REF",
        "c2",
        "u2",
        "snp_b",
        "REF",
        "c3",
        "u3",
        "snp_a",
        "REF",
        "c3",
        "u3",
        "snp_b",
        "ALT",
        "c4",
        "u4",
        "snp_a",
        "REF",
        "c4",
        "u4",
        "snp_b",
        "ALT",
        "c5",
        "u5",
        "snp_a",
        "REF",
        "c5",
        "u5",
        "snp_b",
        "REF"
    )
    # 3/5 "same", 2/5 "opposite" -- consistency 0.6, below the 0.9 threshold

    result <- phase_snps(per_snp, min_molecules = 5, min_consistency = 0.9)

    # Verify an inconsistent edge is not phased
    expect_equal(nrow(result), 0)
})

test_that("phase_snps() keeps disjoint SNP pairs in separate blocks", {
    per_snp <- dplyr::bind_rows(
        lapply(1:6, function(i) {
            tibble::tibble(
                barcode = paste0("c", i),
                umi = paste0("u", i),
                snp_id = c("snp_a", "snp_b"),
                allele = "REF"
            )
        }),
        lapply(1:6, function(i) {
            tibble::tibble(
                barcode = paste0("d", i),
                umi = paste0("v", i),
                snp_id = c("snp_x", "snp_y"),
                allele = "REF"
            )
        })
    )

    result <- phase_snps(per_snp, min_molecules = 5, min_consistency = 0.9)

    # Verify two unconnected SNP pairs form two separate phase blocks
    expect_equal(dplyr::n_distinct(result$block), 2)
    blocks <- split(result$snp_id, result$block)
    # Check that snp_a/snp_b never share a block with snp_x/snp_y
    expect_true(any(vapply(blocks, function(b) setequal(b, c("snp_a", "snp_b")), logical(1))))
    expect_true(any(vapply(blocks, function(b) setequal(b, c("snp_x", "snp_y")), logical(1))))
})

test_that("phase_snps() returns an empty tibble when no molecule spans multiple SNPs", {
    per_snp <- tibble::tribble(
        ~barcode,
        ~umi,
        ~snp_id,
        ~allele,
        "c1",
        "u1",
        "snp_a",
        "REF",
        "c2",
        "u2",
        "snp_b",
        "REF"
    )

    result <- phase_snps(per_snp)

    # Verify the result has the expected empty structure, not an error
    expect_equal(nrow(result), 0)
    expect_named(result, c("snp_id", "block", "allele_on_h1"))
})

# ==============================================================================
# Test: assign_snp_genes()
# ==============================================================================

test_that("assign_snp_genes() keeps a SNP overlapping exactly one gene", {
    snp_info <- tibble::tibble(snp_id = "snp1", chrom = "chr1", pos = 100L)
    gene_anno <- tibble::tibble(chrom = "chr1", start = 50L, end = 150L, gene_name = "GENE1", strand = "+")

    result <- assign_snp_genes(snp_info, gene_anno)

    # Verify the singly-overlapping SNP is assigned to its one gene
    expect_equal(result$gene_name[result$snp_id == "snp1"], "GENE1")
    # Verify a singly-overlapping SNP is not flagged as requiring strand resolution
    expect_false(result$ambiguous[result$snp_id == "snp1"])
})

test_that("assign_snp_genes() drops a SNP overlapping two genes on the same strand", {
    snp_info <- tibble::tibble(snp_id = "snp1", chrom = "chr1", pos = 100L)
    gene_anno <- tibble::tibble(
        chrom = c("chr1", "chr1"),
        start = c(50L, 90L),
        end = c(150L, 200L),
        gene_name = c("GENE1", "GENE2"),
        strand = c("+", "+")
    )

    result <- assign_snp_genes(snp_info, gene_anno)

    # Verify a SNP overlapping two same-strand gene bodies is excluded, not comma-joined
    expect_false("snp1" %in% result$snp_id)
})

test_that("assign_snp_genes() keeps both candidates of a SNP overlapping two opposite-strand genes", {
    snp_info <- tibble::tibble(snp_id = "snp1", chrom = "chr1", pos = 100L)
    gene_anno <- tibble::tibble(
        chrom = c("chr1", "chr1"),
        start = c(50L, 90L),
        end = c(150L, 200L),
        gene_name = c("GENE1", "GENE2"),
        strand = c("+", "-")
    )

    result <- assign_snp_genes(snp_info, gene_anno)

    # Verify both strand-resolvable candidates are retained, not dropped
    expect_setequal(result$gene_name, c("GENE1", "GENE2"))
    # Verify each candidate carries its own gene's strand
    expect_equal(result$gene_strand[result$gene_name == "GENE1"], "+")
    expect_equal(result$gene_strand[result$gene_name == "GENE2"], "-")
    # Verify both are flagged as requiring molecule-level strand resolution
    expect_true(all(result$ambiguous))
})

test_that("assign_snp_genes() drops candidates that still share a strand with another candidate", {
    # Three overlapping genes: GENE1/GENE2 share "+" (irresolvable even by strand),
    # GENE3 alone on "-" (resolvable).
    snp_info <- tibble::tibble(snp_id = "snp1", chrom = "chr1", pos = 100L)
    gene_anno <- tibble::tibble(
        chrom = rep("chr1", 3),
        start = c(50L, 60L, 90L),
        end = c(150L, 160L, 200L),
        gene_name = c("GENE1", "GENE2", "GENE3"),
        strand = c("+", "+", "-")
    )

    result <- assign_snp_genes(snp_info, gene_anno)

    # Verify only the strand-unique candidate survives
    expect_equal(result$gene_name, "GENE3")
    expect_true(result$ambiguous)
})

test_that("assign_snp_genes() drops a SNP overlapping no gene", {
    snp_info <- tibble::tibble(snp_id = "snp1", chrom = "chr1", pos = 1000L)
    gene_anno <- tibble::tibble(chrom = "chr1", start = 50L, end = 150L, gene_name = "GENE1", strand = "+")

    result <- assign_snp_genes(snp_info, gene_anno)

    # Verify a non-overlapping SNP does not appear in the result
    expect_false("snp1" %in% result$snp_id)
    expect_equal(nrow(result), 0)
})

test_that("assign_snp_genes() errors on missing required columns", {
    snp_info <- tibble::tibble(snp_id = "snp1", chrom = "chr1")
    gene_anno <- tibble::tibble(chrom = "chr1", start = 50L, end = 150L, gene_name = "GENE1", strand = "+")

    # Verify a clear error is raised when snp_info lacks 'pos'
    expect_error(assign_snp_genes(snp_info, gene_anno), "missing required column")
})

test_that("assign_snp_genes() errors when gene_anno lacks strand", {
    snp_info <- tibble::tibble(snp_id = "snp1", chrom = "chr1", pos = 100L)
    gene_anno <- tibble::tibble(chrom = "chr1", start = 50L, end = 150L, gene_name = "GENE1")

    # Verify a clear error is raised when gene_anno lacks 'strand'
    expect_error(assign_snp_genes(snp_info, gene_anno), "missing required column")
})

# ==============================================================================
# Test: extract_snp_calls() (integration, real BAM fixture)
# ==============================================================================

test_that("extract_snp_calls() reproduces Rsamtools::pileup()'s base count up to CIGAR no-calls", {
    skip_if_not(file.exists(test_bam), "TASL.bam fixture not found")
    skip_if_not(nzchar(Sys.which("samtools")) || TRUE) # exercised via Rsamtools either way

    snp_info <- tibble::tibble(chrom = "chrX", pos = 30169685L, ref = "A", alt = "C", snp_id = "test_snp")
    extracted <- extract_snp_calls(test_bam, snp_info, min_mapq = 20, min_baseq = 10)

    total_calls <- sum(extracted$tallies$n_calls)
    # Verify total resolved base calls sit at or below the number of overlapping
    # alignments (2603, independently confirmed via readGAlignments with the
    # same flag/mapq filters) -- some alignments overlap the window but do not
    # cover the exact base (intron/deletion), which is a genuine no-call
    expect_lte(total_calls, 2603)
    # Check this fixture position is well covered, not a fluke near-zero count
    expect_gt(total_calls, 2500)
    # Confirm no OTH calls at a position with no observed non-REF/ALT bases in
    # Rsamtools::pileup() at min_baseq = 10
    expect_false("OTH" %in% extracted$tallies$allele)
})

test_that("extract_snp_calls() restricts to the requested barcodes", {
    skip_if_not(file.exists(test_bam), "TASL.bam fixture not found")

    snp_info <- tibble::tibble(chrom = "chrX", pos = 30169685L, ref = "A", alt = "C", snp_id = "test_snp")
    full <- extract_snp_calls(test_bam, snp_info, min_mapq = 20, min_baseq = 10)
    one_barcode <- full$tallies$barcode[1]

    restricted <- extract_snp_calls(test_bam, snp_info, barcodes = one_barcode, min_mapq = 20, min_baseq = 10)

    # Verify only the requested barcode appears in the restricted tallies
    expect_true(all(restricted$tallies$barcode == one_barcode))
    # Check restricting barcodes strictly reduces (or matches) the call count
    expect_lte(sum(restricted$tallies$n_calls), sum(full$tallies$n_calls))
})

test_that("extract_snp_calls() errors on missing required snp_info columns", {
    snp_info <- tibble::tibble(chrom = "chrX", pos = 1L)

    # Verify a clear error names the missing columns rather than failing deep in Bioconductor code
    expect_error(extract_snp_calls("nonexistent.bam", snp_info), "missing required column")
})

test_that("extract_snp_calls() errors on empty snp_info", {
    snp_info <- tibble::tibble(
        snp_id = character(),
        chrom = character(),
        pos = integer(),
        ref = character(),
        alt = character()
    )

    # Verify an empty SNP set is rejected up front rather than producing a confusing downstream error
    expect_error(extract_snp_calls("nonexistent.bam", snp_info), "no rows")
})

# ==============================================================================
# Test: .orient_phase_blocks()
# ==============================================================================

test_that(".orient_phase_blocks() orients a block when all anchors agree H1 = X1", {
    phase <- tibble::tribble(
        ~snp_id,
        ~block,
        ~allele_on_h1,
        "anchor1",
        1L,
        "REF",
        "anchor2",
        1L,
        "REF",
        "bridge",
        1L,
        "ALT"
    )
    anchors <- tibble::tribble(
        ~snp_id,
        ~allele_on_x1_em,
        "anchor1",
        "REF",
        "anchor2",
        "REF"
    )

    result <- .orient_phase_blocks(phase, anchors, donor = "donor0")

    # Verify anchor SNPs keep their EM-agreeing molecule orientation
    expect_equal(result$allele_on_x1_molecule[result$snp_id == "anchor1"], "REF")
    # Check the non-anchor bridging SNP is still oriented via the block consensus
    expect_equal(result$allele_on_x1_molecule[result$snp_id == "bridge"], "ALT")
    # Confirm no SNP in a clean, all-agreeing block is flagged as conflicting
    expect_false(any(result$phase_conflict))
})

test_that(".orient_phase_blocks() flips the block when all anchors agree H1 = X2", {
    phase <- tibble::tribble(
        ~snp_id,
        ~block,
        ~allele_on_h1,
        "anchor1",
        1L,
        "ALT",
        "bridge",
        1L,
        "REF"
    )
    anchors <- tibble::tribble(
        ~snp_id,
        ~allele_on_x1_em,
        "anchor1",
        "REF"
    )
    # anchor1's H1 (ALT) disagrees with its own allele_on_x1_em (REF), so H1 = X2

    result <- .orient_phase_blocks(phase, anchors, donor = "donor0")

    # Verify the anchor's own molecule-derived value is the flip of its H1 label
    expect_equal(result$allele_on_x1_molecule[result$snp_id == "anchor1"], "REF")
    # Check the bridging SNP is flipped consistently with the same orientation
    expect_equal(result$allele_on_x1_molecule[result$snp_id == "bridge"], "ALT")
})

test_that(".orient_phase_blocks() leaves a zero-anchor block unoriented without flagging conflict", {
    phase <- tibble::tribble(
        ~snp_id,
        ~block,
        ~allele_on_h1,
        "snp1",
        1L,
        "REF",
        "snp2",
        1L,
        "ALT"
    )
    anchors <- tibble::tribble(~snp_id, ~allele_on_x1_em)[0, ]

    result <- .orient_phase_blocks(phase, anchors, donor = "donor0")

    # Verify no molecule-based orientation is assigned without a reachable anchor
    expect_true(all(is.na(result$allele_on_x1_molecule)))
    # Check the lack of an anchor is not itself treated as a conflict
    expect_false(any(result$phase_conflict))
})

test_that(".orient_phase_blocks() resolves via majority when one of >=3 anchors disagrees", {
    phase <- tibble::tribble(
        ~snp_id,
        ~block,
        ~allele_on_h1,
        "a1",
        1L,
        "REF",
        "a2",
        1L,
        "REF",
        "a3",
        1L,
        "ALT", # outlier: disagrees with a1/a2's implied orientation
        "bridge",
        1L,
        "REF"
    )
    anchors <- tibble::tribble(
        ~snp_id,
        ~allele_on_x1_em,
        "a1",
        "REF",
        "a2",
        "REF",
        "a3",
        "REF"
    )
    # a1: H1(REF)==EM(REF) -> H1 is X1. a2: same -> H1 is X1. a3: H1(ALT)!=EM(REF) -> H1 is X2 (outlier)

    result <- .orient_phase_blocks(phase, anchors, donor = "donor0")

    # Verify the majority orientation (H1 = X1) is used for the bridging SNP
    expect_equal(result$allele_on_x1_molecule[result$snp_id == "bridge"], "REF")
    # Check the majority-agreeing anchors are not flagged
    expect_false(result$phase_conflict[result$snp_id == "a1"])
    expect_false(result$phase_conflict[result$snp_id == "a2"])
    # Confirm the single dissenting anchor is flagged for inspection
    expect_true(result$phase_conflict[result$snp_id == "a3"])
    # Verify the outlier still receives the majority-consensus molecule orientation
    expect_equal(result$allele_on_x1_molecule[result$snp_id == "a3"], "ALT")
})

test_that(".orient_phase_blocks() refuses to resolve when exactly two anchors disagree", {
    phase <- tibble::tribble(
        ~snp_id,
        ~block,
        ~allele_on_h1,
        "a1",
        1L,
        "REF",
        "a2",
        1L,
        "ALT",
        "bridge",
        1L,
        "REF"
    )
    anchors <- tibble::tribble(
        ~snp_id,
        ~allele_on_x1_em,
        "a1",
        "REF",
        "a2",
        "REF"
    )
    # a1: H1 is X1. a2: H1(ALT) != EM(REF) -> H1 is X2. Two anchors, no majority.

    result <- .orient_phase_blocks(phase, anchors, donor = "donor0")

    # Verify no SNP in the unresolved component gets a molecule-based orientation
    expect_true(all(is.na(result$allele_on_x1_molecule)))
    # Check every SNP in the conflicting component is flagged, not just the two anchors
    expect_true(all(result$phase_conflict))
})

test_that(".orient_phase_blocks() refuses to resolve an even split among >=4 anchors", {
    phase <- tibble::tribble(
        ~snp_id,
        ~block,
        ~allele_on_h1,
        "a1",
        1L,
        "REF",
        "a2",
        1L,
        "REF",
        "a3",
        1L,
        "ALT",
        "a4",
        1L,
        "ALT"
    )
    anchors <- tibble::tribble(
        ~snp_id,
        ~allele_on_x1_em,
        "a1",
        "REF",
        "a2",
        "REF",
        "a3",
        "REF",
        "a4",
        "REF"
    )
    # a1,a2 imply H1=X1; a3,a4 imply H1=X2 -- an even 2-vs-2 split, not a single outlier

    result <- .orient_phase_blocks(phase, anchors, donor = "donor0")

    # Verify a 2-vs-2 split is treated as an unresolved conflict, not a forced majority
    expect_true(all(is.na(result$allele_on_x1_molecule)))
    expect_true(all(result$phase_conflict))
})

test_that(".orient_phase_blocks() keeps separate blocks independent", {
    phase <- tibble::tribble(
        ~snp_id,
        ~block,
        ~allele_on_h1,
        "a1",
        1L,
        "REF",
        "bridge1",
        1L,
        "ALT",
        "a2",
        2L,
        "ALT",
        "bridge2",
        2L,
        "ALT"
    )
    anchors <- tibble::tribble(
        ~snp_id,
        ~allele_on_x1_em,
        "a1",
        "REF",
        "a2",
        "REF"
    )
    # Block 1: a1 implies H1=X1. Block 2: a2's H1(ALT) != EM(REF), implies H1=X2.

    result <- .orient_phase_blocks(phase, anchors, donor = "donor0")

    # Verify block 1's bridging SNP is oriented per block 1's own anchor
    expect_equal(result$allele_on_x1_molecule[result$snp_id == "bridge1"], "ALT")
    # Check block 2's bridging SNP is oriented per block 2's own (independent) anchor
    expect_equal(result$allele_on_x1_molecule[result$snp_id == "bridge2"], "REF")
    # Confirm neither block's resolution was affected by the other
    expect_false(any(result$phase_conflict))
})

test_that(".orient_phase_blocks() gives an unlinked anchor its own singleton block", {
    # "lonely" never appears in phase at all -- phase_snps() never linked it
    # to any partner, e.g. a gene with only one heterozygous SNP.
    phase <- tibble::tribble(
        ~snp_id,
        ~block,
        ~allele_on_h1,
        "a1",
        1L,
        "REF",
        "a2",
        1L,
        "REF"
    )
    anchors <- tibble::tribble(
        ~snp_id,
        ~allele_on_x1_em,
        "a1",
        "REF",
        "a2",
        "REF",
        "lonely",
        "ALT"
    )

    result <- .orient_phase_blocks(phase, anchors, donor = "donor0")

    # Verify the unlinked anchor still appears in the result
    expect_true("lonely" %in% result$snp_id)
    # Check it carries its own EM-derived phase forward unchanged
    expect_equal(result$allele_on_x1_molecule[result$snp_id == "lonely"], "ALT")
    # Confirm it is not flagged as conflicting -- there is nothing to disagree with
    expect_false(result$phase_conflict[result$snp_id == "lonely"])
    # Verify its synthesized block id is distinct from every real phase_snps() block
    lonely_block <- result$phase_block[result$snp_id == "lonely"]
    expect_true(lonely_block < 0)
    expect_false(lonely_block %in% phase$block)
})

test_that(".orient_phase_blocks() gives every unlinked anchor a distinct singleton block", {
    phase <- tibble::tibble(snp_id = character(), block = integer(), allele_on_h1 = character())
    anchors <- tibble::tribble(
        ~snp_id,
        ~allele_on_x1_em,
        "solo1",
        "REF",
        "solo2",
        "ALT"
    )

    result <- .orient_phase_blocks(phase, anchors, donor = "donor0")

    # Verify both unlinked anchors are present with their own phase
    expect_equal(result$allele_on_x1_molecule[result$snp_id == "solo1"], "REF")
    expect_equal(result$allele_on_x1_molecule[result$snp_id == "solo2"], "ALT")
    # Check they never share a synthesized block id with each other
    expect_equal(dplyr::n_distinct(result$phase_block), 2)
})

# ==============================================================================
# Test: add_molecule_phase()
# ==============================================================================

# Builds a minimal, two-donor SNPData with XCI diagnostics injected directly
# (rather than fitting the EM), mirroring make_hap_fixture() in
# test-haplotype-expression.R. snp1 is EM-informative in donor0 (anchor);
# snp2 is a het SNP the EM never touched (allele_on_x1 NA), the kind of SNP
# add_molecule_phase() exists to rescue.
make_phase_fixture <- function() {
    ref <- rbind(
        c(0L, 0L, 10L, 10L),
        c(5L, 5L, 5L, 5L)
    )
    alt <- rbind(
        c(10L, 10L, 0L, 0L),
        c(5L, 5L, 5L, 5L)
    )
    snp_info <- data.frame(
        chrom = "chrX",
        pos = c(1000L, 2000L),
        ref = "A",
        alt = "G",
        stringsAsFactors = FALSE
    )
    barcode_info <- data.frame(
        barcode = paste0("cell", 1:4),
        donor = c("donor0", "donor0", "donor0", "donor0"),
        stringsAsFactors = FALSE
    )
    obj <- SNPData(
        ref_count = Matrix::Matrix(ref, sparse = TRUE),
        alt_count = Matrix::Matrix(alt, sparse = TRUE),
        snp_info = snp_info,
        barcode_info = barcode_info
    )
    snp_ids <- snp_info(obj)$snp_id

    donor_snp_diag <- data.frame(
        snp_id = snp_ids[1],
        donor = "donor0",
        xci_informative = TRUE,
        allele_on_x1 = "REF",
        zygosity = "het",
        zygosity_source = "vireo_gt",
        stringsAsFactors = FALSE
    )
    obj <- add_donor_snp_metadata(obj, donor_snp_diag, join_by = c("snp_id", "donor"), overwrite = TRUE)
    # snp2 is het too, but never reached by the EM -- no allele_on_x1 row for it.
    obj <- add_donor_snp_metadata(
        obj,
        data.frame(snp_id = snp_ids[2], donor = "donor0", zygosity = "het", zygosity_source = "vireo_gt"),
        join_by = c("snp_id", "donor"),
        overwrite = TRUE
    )
    obj <- add_barcode_metadata(
        obj,
        data.frame(cell_id = barcode_info(obj)$cell_id, active_x = c("X2", "X2", "X1", "X1")),
        join_by = "cell_id",
        overwrite = TRUE
    )
    list(obj = obj, snp_ids = snp_ids)
}

# add_molecule_phase() insists its BAM files exist and are indexed before it
# extracts anything, so tests that mock the extraction still need real paths.
# The files stay empty: nothing ever reads them.
local_fake_bam <- function(name = "fake.bam", env = parent.frame()) {
    dir <- withr::local_tempdir(.local_envir = env)
    bam <- file.path(dir, name)
    file.create(bam, paste0(bam, ".bai"))
    # add_molecule_phase() resolves the paths it is given, so return the
    # resolved form: on macOS the temp directory is reached through a symlink
    # and the two spellings would not compare equal.
    normalizePath(bam)
}

# Stubs out the BAM-reading step with an empty extraction, so the surrounding
# validation and control flow can be exercised without a real BAM. Records
# every file it was called for, so tests can assert which BAMs were opened.
local_mocked_extraction <- function(called_for, env = parent.frame()) {
    testthat::local_mocked_bindings(
        extract_snp_calls = function(bam_file, snp_info, barcodes = NULL, ...) {
            assign(
                "calls",
                c(get("calls", envir = called_for), bam_file),
                envir = called_for
            )
            list(
                tallies = tibble::tibble(
                    barcode = character(),
                    umi = character(),
                    snp_id = character(),
                    allele = character(),
                    n_calls = integer()
                ),
                reads = tibble::tibble(
                    barcode = character(),
                    umi = character(),
                    qname = character(),
                    strand = character()
                )
            )
        },
        .infer_bam_strand_orientation = function(bam_file, ...) {
            list(orientation = "sense", n_ts_reads = 500L, concordance = 1, n_scanned = 5000L)
        },
        .package = "snplet",
        .env = env
    )
}

# A place for local_mocked_extraction() to record the BAMs it was called for.
new_call_log <- function() {
    log_env <- new.env(parent = emptyenv())
    log_env$calls <- character()
    log_env
}

test_that(".pool_donor_calls() sums a split molecule's reads before the allele vote", {
    per_file <- list(
        list(
            tallies = tibble::tibble(
                barcode = "c1",
                umi = "u1",
                snp_id = "snp1",
                allele = "REF",
                n_calls = 3L
            ),
            molecule_strand = tibble::tibble(barcode = "c1", umi = "u1", transcript_strand = "+")
        ),
        list(
            tallies = tibble::tibble(
                barcode = c("c1", "c1"),
                umi = c("u1", "u1"),
                snp_id = c("snp1", "snp1"),
                allele = c("REF", "ALT"),
                n_calls = c(1L, 2L)
            ),
            molecule_strand = tibble::tibble(barcode = "c1", umi = "u1", transcript_strand = "+")
        )
    )

    pooled <- .pool_donor_calls(per_file)

    # Verify REF's reads are summed across the two files (3 + 1) rather than
    # taken as the per-file maximum, which would have lost the first file's 3
    expect_equal(pooled$tallies$n_calls[pooled$tallies$allele == "REF"], 4L)
    # Confirm the pooled counts elect REF, the winner over all 6 reads, where
    # the largest single-file tally would have elected ALT
    expect_equal(molecule_snp_alleles(pooled$tallies)$allele, "REF")
})

test_that(".pool_donor_calls() returns one strand row per molecule", {
    strand_table <- tibble::tibble(barcode = "c1", umi = "u1", transcript_strand = "+")
    per_file <- list(
        list(
            tallies = tibble::tibble(
                barcode = "c1",
                umi = "u1",
                snp_id = "snp1",
                allele = "REF",
                n_calls = 1L
            ),
            molecule_strand = strand_table
        ),
        list(
            tallies = tibble::tibble(
                barcode = "c1",
                umi = "u1",
                snp_id = "snp1",
                allele = "REF",
                n_calls = 1L
            ),
            molecule_strand = strand_table
        )
    )

    pooled <- .pool_donor_calls(per_file)

    # Ensure a molecule seen in both files yields a single strand row: the
    # left join onto the calls fans out if it does not
    expect_equal(nrow(pooled$molecule_strand), 1L)
    # Verify the agreed strand survives the majority vote
    expect_equal(pooled$molecule_strand$transcript_strand, "+")
})

test_that(".pool_donor_calls() resolves a cross-file strand disagreement to NA", {
    tallies <- tibble::tibble(barcode = "c1", umi = "u1", snp_id = "snp1", allele = "REF", n_calls = 1L)
    per_file <- list(
        list(
            tallies = tallies,
            molecule_strand = tibble::tibble(
                barcode = "c1",
                umi = "u1",
                transcript_strand = "+"
            )
        ),
        list(
            tallies = tallies,
            molecule_strand = tibble::tibble(
                barcode = "c1",
                umi = "u1",
                transcript_strand = "-"
            )
        )
    )

    pooled <- .pool_donor_calls(per_file)

    # Verify a tied disagreement becomes NA rather than an arbitrary pick,
    # leaving strand-ambiguous SNPs unresolved for that molecule
    expect_true(is.na(pooled$molecule_strand$transcript_strand))
})

test_that(".pool_donor_calls() lets a calibrated file outvote an uncalibrated one", {
    tallies <- tibble::tibble(barcode = "c1", umi = "u1", snp_id = "snp1", allele = "REF", n_calls = 1L)
    per_file <- list(
        list(
            tallies = tallies,
            molecule_strand = tibble::tibble(
                barcode = "c1",
                umi = "u1",
                transcript_strand = NA_character_
            )
        ),
        list(
            tallies = tallies,
            molecule_strand = tibble::tibble(
                barcode = "c1",
                umi = "u1",
                transcript_strand = "-"
            )
        )
    )

    pooled <- .pool_donor_calls(per_file)

    # Confirm an NA from a file whose orientation could not be inferred is an
    # absence of evidence, not a vote against the file that does know
    expect_equal(pooled$molecule_strand$transcript_strand, "-")
})

test_that("add_molecule_phase() errors when no XCI diagnostics are stored", {
    ref <- Matrix::Matrix(matrix(1L, 2, 2), sparse = TRUE)
    alt <- Matrix::Matrix(matrix(1L, 2, 2), sparse = TRUE)
    snp_info <- data.frame(chrom = "chrX", pos = c(1L, 2L), ref = "A", alt = "G")
    barcode_info <- data.frame(barcode = c("c1", "c2"), donor = c("donor0", "donor0"))
    obj <- SNPData(ref_count = ref, alt_count = alt, snp_info = snp_info, barcode_info = barcode_info)

    # Verify the function refuses to run before assign_xci() has stored a fit
    expect_error(add_molecule_phase(obj, bam_files = c(donor0 = "x.bam")), "Run assign_xci")
})

test_that("add_molecule_phase() falls back to the BAM paths recorded on the object", {
    fixture <- make_phase_fixture()
    obj <- add_barcode_metadata(
        fixture$obj,
        data.frame(cell_id = barcode_info(fixture$obj)$cell_id, library_id = "lib_A"),
        join_by = "cell_id",
        overwrite = TRUE
    )
    bam <- local_fake_bam()
    obj <- add_library_bams(obj, c(lib_A = bam))
    call_log <- new_call_log()
    local_mocked_extraction(call_log)

    add_molecule_phase(obj)

    # Verify a path recorded at import is used without being passed again,
    # which is the point of storing it on the object
    expect_identical(call_log$calls, bam)
})

test_that("add_molecule_phase() errors when no BAM paths are given or recorded", {
    fixture <- make_phase_fixture()

    # Verify the object says where the paths should come from rather than
    # failing on a missing argument
    expect_error(add_molecule_phase(fixture$obj), "none recorded on the object")
})

test_that("add_molecule_phase() errors on unrecognised library names in bam_files", {
    fixture <- make_phase_fixture()
    obj <- add_barcode_metadata(
        fixture$obj,
        data.frame(cell_id = barcode_info(fixture$obj)$cell_id, library_id = "lib_A"),
        join_by = "cell_id",
        overwrite = TRUE
    )

    # Verify a BAM keyed to a library absent from barcode_info is rejected up
    # front, rather than silently leaving that library's donors unphased
    expect_error(
        add_molecule_phase(obj, bam_files = c(lib_B = "x.bam")),
        "not found in barcode_info\\$library_id"
    )
})

test_that("add_molecule_phase() errors when one donor's cells span two libraries", {
    fixture <- make_phase_fixture()
    obj <- add_barcode_metadata(
        fixture$obj,
        data.frame(
            cell_id = barcode_info(fixture$obj)$cell_id,
            library_id = c("lib_A", "lib_A", "lib_B", "lib_B")
        ),
        join_by = "cell_id",
        overwrite = TRUE
    )

    # Verify a donor split across libraries is rejected: its BAM files are
    # looked up by library, so there is no single correct set to read
    expect_error(
        add_molecule_phase(obj, bam_files = c(lib_A = "x.bam", lib_B = "y.bam")),
        "more than one library"
    )
})

test_that("add_molecule_phase() errors when only some cells carry a library_id", {
    fixture <- make_phase_fixture()
    obj <- add_barcode_metadata(
        fixture$obj,
        data.frame(
            cell_id = barcode_info(fixture$obj)$cell_id,
            library_id = c("lib_A", "lib_A", NA, NA)
        ),
        join_by = "cell_id",
        overwrite = TRUE
    )

    # Verify a partly-labelled object is rejected rather than guessed at: the
    # unlabelled cells could belong to lib_A or to a library with no BAM at all
    expect_error(
        add_molecule_phase(obj, bam_files = c(lib_A = "x.bam")),
        "set for some cells but not others"
    )
})

test_that("add_molecule_phase() errors when an unlabelled object is given several libraries", {
    fixture <- make_phase_fixture()

    # Check that an object with no library labels cannot be handed more than
    # one library's BAMs, since nothing records which donor belongs to which
    expect_error(
        add_molecule_phase(fixture$obj, bam_files = c(lib_A = "x.bam", lib_B = "y.bam")),
        "must have exactly one entry"
    )
})

test_that("add_molecule_phase() errors when a library lists the same BAM twice", {
    fixture <- make_phase_fixture()
    bam <- local_fake_bam()

    # Verify a repeated file is rejected: it would be extracted twice and its
    # tallies summed, doubling every read behind that library's molecules
    expect_error(
        add_molecule_phase(fixture$obj, bam_files = list(lib_A = c(bam, bam))),
        "listed more than once"
    )
})

test_that("add_molecule_phase() errors when a BAM has no index", {
    fixture <- make_phase_fixture()
    dir <- withr::local_tempdir()
    unindexed <- file.path(dir, "unindexed.bam")
    file.create(unindexed)

    # Verify a missing index is a hard error: extraction seeks by index, and
    # without one the region-restricted scan degrades to reading the whole file
    expect_error(
        add_molecule_phase(fixture$obj, bam_files = c(lib_A = unindexed)),
        "no index"
    )
})

test_that("add_molecule_phase() accepts an object whose cells are all unlabelled", {
    fixture <- make_phase_fixture()
    bam <- local_fake_bam()
    call_log <- new_call_log()
    local_mocked_extraction(call_log)

    # Confirm an all-NA library_id is treated as one implicit library, so
    # objects imported without a library label still phase as before
    expect_no_error(add_molecule_phase(fixture$obj, bam_files = c(any_name = bam)))
    # Verify the single supplied BAM was the one opened for the only donor
    expect_identical(call_log$calls, bam)
})

test_that("add_molecule_phase() opens every BAM listed for a donor's library", {
    fixture <- make_phase_fixture()
    obj <- add_barcode_metadata(
        fixture$obj,
        data.frame(cell_id = barcode_info(fixture$obj)$cell_id, library_id = "lib_A"),
        join_by = "cell_id",
        overwrite = TRUE
    )
    dir <- withr::local_tempdir()
    bams <- file.path(dir, c("a.bam", "b.bam"))
    file.create(bams, paste0(bams, ".bai"))
    bams <- normalizePath(bams)
    call_log <- new_call_log()
    local_mocked_extraction(call_log)

    add_molecule_phase(obj, bam_files = list(lib_A = bams))

    # Verify both of the library's files were extracted for its donor, so a
    # molecule split across them can be pooled rather than half-counted
    expect_setequal(call_log$calls, bams)
})

test_that("add_molecule_phase() excludes the 'doublet' and 'unassigned' donor labels", {
    # A fixture where "doublet" is a real, valid donor label in barcode_info
    # (as Vireo emits it). Donors are now derived from the object rather than
    # named by the caller, so the exclusion has to be deliberate: nothing else
    # would stop these two being phased like any other donor.
    ref <- rbind(c(0L, 0L, 10L, 10L, 5L, 5L))
    alt <- rbind(c(10L, 10L, 0L, 0L, 5L, 5L))
    snp_info <- data.frame(chrom = "chrX", pos = 1000L, ref = "A", alt = "G", stringsAsFactors = FALSE)
    barcode_info <- data.frame(
        barcode = paste0("cell", 1:6),
        donor = c("donor0", "donor0", "donor0", "donor0", "doublet", "unassigned"),
        stringsAsFactors = FALSE
    )
    obj <- SNPData(
        ref_count = Matrix::Matrix(ref, sparse = TRUE),
        alt_count = Matrix::Matrix(alt, sparse = TRUE),
        snp_info = snp_info,
        barcode_info = barcode_info
    )
    snp_id <- snp_info(obj)$snp_id
    obj <- add_donor_snp_metadata(
        obj,
        data.frame(
            snp_id = snp_id,
            donor = "donor0",
            xci_informative = TRUE,
            allele_on_x1 = "REF",
            zygosity = "het",
            zygosity_source = "vireo_gt",
            stringsAsFactors = FALSE
        ),
        join_by = c("snp_id", "donor"),
        overwrite = TRUE
    )
    obj <- add_barcode_metadata(
        obj,
        data.frame(cell_id = barcode_info(obj)$cell_id, active_x = c("X2", "X2", "X1", "X1", NA, NA)),
        join_by = "cell_id",
        overwrite = TRUE
    )

    fake_tallies <- tibble::tibble(
        barcode = character(),
        umi = character(),
        snp_id = character(),
        allele = character(),
        n_calls = integer()
    )
    fake_reads <- tibble::tibble(
        barcode = character(),
        umi = character(),
        qname = character(),
        strand = character()
    )
    called_for <- character()
    testthat::local_mocked_bindings(
        extract_snp_calls = function(bam_file, snp_info, barcodes = NULL, ...) {
            called_for <<- c(called_for, bam_file)
            list(tallies = fake_tallies, reads = fake_reads)
        },
        .infer_bam_strand_orientation = function(bam_file, ...) {
            list(orientation = "sense", n_ts_reads = 500L, concordance = 1, n_scanned = 5000L)
        },
        .package = "snplet"
    )
    bam <- local_fake_bam()

    # The package test suite raises the global log threshold to FATAL (see
    # tests/testthat/setup.R), so capturing the WARN-level exclusion message
    # requires lowering it locally, per the pattern in test-SNPData-methods.R.
    log_file <- withr::local_tempfile(pattern = "add_molecule_phase_warn_", fileext = ".log")
    original_appender_name <- as.character(logger::log_appender())
    original_appender <- get(original_appender_name, asNamespace("logger"))
    original_threshold <- logger::log_threshold()
    logger::log_appender(logger::appender_file(log_file))
    logger::log_threshold(logger::WARN)
    withr::defer(logger::log_appender(original_appender))
    withr::defer(logger::log_threshold(original_threshold))

    add_molecule_phase(obj, bam_files = c(lib_A = bam))

    # Verify the log names the excluded non-donor labels
    log_lines <- readLines(log_file, warn = FALSE)
    expect_true(any(grepl("doublet, unassigned", log_lines)))
    # Check the BAM extraction step itself was only ever invoked once -- for
    # the real donor, never for "doublet" or "unassigned"
    expect_length(called_for, 1)
})

test_that("add_molecule_phase() returns x unchanged when every donor is doublet/unassigned", {
    ref <- rbind(c(5L, 5L))
    alt <- rbind(c(5L, 5L))
    snp_info <- data.frame(chrom = "chrX", pos = 1000L, ref = "A", alt = "G", stringsAsFactors = FALSE)
    barcode_info <- data.frame(
        barcode = c("cell1", "cell2"),
        donor = c("doublet", "unassigned"),
        stringsAsFactors = FALSE
    )
    obj <- SNPData(
        ref_count = Matrix::Matrix(ref, sparse = TRUE),
        alt_count = Matrix::Matrix(alt, sparse = TRUE),
        snp_info = snp_info,
        barcode_info = barcode_info
    )
    snp_id <- snp_info(obj)$snp_id
    # xci_informative needs at least one TRUE for .has_xci_diagnostics() to pass,
    # even though it is stored against a donor that will never be processed here.
    obj <- add_donor_snp_metadata(
        obj,
        data.frame(
            snp_id = snp_id,
            donor = "doublet",
            xci_informative = TRUE,
            allele_on_x1 = "REF",
            stringsAsFactors = FALSE
        ),
        join_by = c("snp_id", "donor"),
        overwrite = TRUE
    )

    # Verify no error, and the object is returned unchanged when nothing real is left to process
    result <- add_molecule_phase(obj, bam_files = c(lib_A = "fake.bam"))
    expect_identical(donor_snp_info(result), donor_snp_info(obj))
})

test_that("add_molecule_phase() phases a donor's split BAMs as if they were one file", {
    fixture <- make_phase_fixture()
    snp_ids <- fixture$snp_ids

    # Five molecules, each spanning both SNPs on the same haplotype. The same
    # reads are presented two ways: whole in one file, or split 2/3 across two
    # files of the same library, which is the arrangement pooling exists for.
    calls <- tibble::tibble(
        barcode = rep(paste0("c", 1:5), each = 2),
        umi = rep(paste0("u", 1:5), each = 2),
        snp_id = rep(snp_ids, times = 5),
        allele = rep(c("REF", "ALT"), times = 5)
    )
    reads <- tibble::tibble(
        barcode = paste0("c", 1:5),
        umi = paste0("u", 1:5),
        qname = paste0("read", 1:5),
        strand = "+"
    )
    whole <- dplyr::mutate(calls, n_calls = 5L)
    first_part <- dplyr::mutate(calls, n_calls = 2L)
    second_part <- dplyr::mutate(calls, n_calls = 3L)

    dir <- withr::local_tempdir()
    bams <- file.path(dir, c("whole.bam", "part_a.bam", "part_b.bam"))
    file.create(bams, paste0(bams, ".bai"))
    bams <- normalizePath(bams)

    testthat::local_mocked_bindings(
        extract_snp_calls = function(bam_file, ...) {
            tallies <- switch(
                basename(bam_file),
                part_a.bam = first_part,
                part_b.bam = second_part,
                whole
            )
            list(tallies = tallies, reads = reads)
        },
        .infer_bam_strand_orientation = function(bam_file, ...) {
            list(orientation = "sense", n_ts_reads = 500L, concordance = 1, n_scanned = 5000L)
        },
        .package = "snplet"
    )

    single <- add_molecule_phase(fixture$obj, bam_files = c(lib_A = bams[1]))
    split <- add_molecule_phase(fixture$obj, bam_files = list(lib_A = bams[2:3]))

    # Verify splitting a donor's reads across two of its library's BAM files
    # changes nothing about the phase it ends up with
    expect_identical(donor_snp_info(split), donor_snp_info(single))
    # Confirm the pooled molecule calls match too, so each split molecule was
    # reassembled with all of its reads rather than counted once per file
    expect_equal(
        dplyr::arrange(attr(split, "molecule_calls"), barcode, umi, snp_id),
        dplyr::arrange(attr(single, "molecule_calls"), barcode, umi, snp_id)
    )
})

test_that("add_molecule_phase() records the strand calibration of every BAM it scanned", {
    fixture <- make_phase_fixture()
    snp_ids <- fixture$snp_ids
    tallies <- tibble::tibble(
        barcode = rep(paste0("c", 1:5), each = 2),
        umi = rep(paste0("u", 1:5), each = 2),
        snp_id = rep(snp_ids, times = 5),
        allele = rep(c("REF", "ALT"), times = 5),
        n_calls = 5L
    )
    reads <- tibble::tibble(
        barcode = paste0("c", 1:5),
        umi = paste0("u", 1:5),
        qname = paste0("read", 1:5),
        strand = "+"
    )
    bam <- local_fake_bam()
    testthat::local_mocked_bindings(
        extract_snp_calls = function(...) list(tallies = tallies, reads = reads),
        .infer_bam_strand_orientation = function(bam_file, ...) {
            list(orientation = "antisense", n_ts_reads = 500L, concordance = 0.99, n_scanned = 5000L)
        },
        .package = "snplet"
    )

    calibration <- attr(add_molecule_phase(fixture$obj, bam_files = c(lib_A = bam)), "bam_calibration")

    # Verify the orientation applied to this file's molecules is recorded on
    # the object, so a misread strand can be diagnosed after the fact
    expect_equal(calibration$bam_file, bam)
    expect_equal(calibration$orientation, "antisense")
    expect_equal(calibration$concordance, 0.99)
})

test_that("add_molecule_phase() never overwrites an existing EM-derived allele_on_x1", {
    fixture <- make_phase_fixture()
    snp_ids <- fixture$snp_ids

    # Mock the BAM-extraction step so this test needs no real BAM file: molecules
    # phase snp1 (the EM anchor) and snp2 (the EM-uninformative SNP) as the same
    # haplotype, oriented so H1 = X1 -- allowing snp2 to be rescued.
    fake_tallies <- tibble::tribble(
        ~barcode,
        ~umi,
        ~snp_id,
        ~allele,
        ~n_calls,
        "c1",
        "u1",
        snp_ids[1],
        "REF",
        5L,
        "c1",
        "u1",
        snp_ids[2],
        "ALT",
        5L,
        "c2",
        "u2",
        snp_ids[1],
        "REF",
        5L,
        "c2",
        "u2",
        snp_ids[2],
        "ALT",
        5L,
        "c3",
        "u3",
        snp_ids[1],
        "REF",
        5L,
        "c3",
        "u3",
        snp_ids[2],
        "ALT",
        5L,
        "c4",
        "u4",
        snp_ids[1],
        "REF",
        5L,
        "c4",
        "u4",
        snp_ids[2],
        "ALT",
        5L,
        "c5",
        "u5",
        snp_ids[1],
        "REF",
        5L,
        "c5",
        "u5",
        snp_ids[2],
        "ALT",
        5L
    )
    fake_reads <- tibble::tibble(
        barcode = paste0("c", 1:5),
        umi = paste0("u", 1:5),
        qname = paste0("read", 1:5),
        strand = "+"
    )
    testthat::local_mocked_bindings(
        extract_snp_calls = function(...) list(tallies = fake_tallies, reads = fake_reads),
        .infer_bam_strand_orientation = function(bam_file, ...) {
            list(orientation = "sense", n_ts_reads = 500L, concordance = 1, n_scanned = 5000L)
        },
        .package = "snplet"
    )

    result <- add_molecule_phase(fixture$obj, bam_files = c(lib_A = local_fake_bam()))
    donor_snp_info <- donor_snp_info(result)

    # Verify the pre-existing EM-derived phase for the anchor SNP is unchanged
    expect_equal(
        donor_snp_info$allele_on_x1[donor_snp_info$snp_id == snp_ids[1] & donor_snp_info$donor == "donor0"],
        "REF"
    )
    # Check the EM-uninformative SNP is rescued with a molecule-derived phase
    expect_equal(
        donor_snp_info$allele_on_x1[donor_snp_info$snp_id == snp_ids[2] & donor_snp_info$donor == "donor0"],
        "ALT"
    )
    # Confirm the new provenance columns are present and consistent
    expect_false(donor_snp_info$phase_conflict[donor_snp_info$snp_id == snp_ids[1] & donor_snp_info$donor == "donor0"])
    expect_equal(
        donor_snp_info$phase_source[donor_snp_info$snp_id == snp_ids[2] & donor_snp_info$donor == "donor0"],
        "read_backed"
    )
})

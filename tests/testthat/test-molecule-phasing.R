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
    gene_anno <- tibble::tibble(chrom = "chr1", start = 50L, end = 150L, gene_name = "GENE1")

    result <- assign_snp_genes(snp_info, gene_anno)

    # Verify the singly-overlapping SNP is assigned to its one gene
    expect_equal(result$gene_name[result$snp_id == "snp1"], "GENE1")
})

test_that("assign_snp_genes() drops a SNP overlapping two genes", {
    snp_info <- tibble::tibble(snp_id = "snp1", chrom = "chr1", pos = 100L)
    gene_anno <- tibble::tibble(
        chrom = c("chr1", "chr1"),
        start = c(50L, 90L),
        end = c(150L, 200L),
        gene_name = c("GENE1", "GENE2")
    )

    result <- assign_snp_genes(snp_info, gene_anno)

    # Verify a SNP overlapping two gene bodies is excluded, not comma-joined
    expect_false("snp1" %in% result$snp_id)
})

test_that("assign_snp_genes() drops a SNP overlapping no gene", {
    snp_info <- tibble::tibble(snp_id = "snp1", chrom = "chr1", pos = 1000L)
    gene_anno <- tibble::tibble(chrom = "chr1", start = 50L, end = 150L, gene_name = "GENE1")

    result <- assign_snp_genes(snp_info, gene_anno)

    # Verify a non-overlapping SNP does not appear in the result
    expect_false("snp1" %in% result$snp_id)
    expect_equal(nrow(result), 0)
})

test_that("assign_snp_genes() errors on missing required columns", {
    snp_info <- tibble::tibble(snp_id = "snp1", chrom = "chr1")
    gene_anno <- tibble::tibble(chrom = "chr1", start = 50L, end = 150L, gene_name = "GENE1")

    # Verify a clear error is raised when snp_info lacks 'pos'
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
    snp_ids <- get_snp_info(obj)$snp_id

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
        data.frame(cell_id = get_barcode_info(obj)$cell_id, active_x = c("X2", "X2", "X1", "X1")),
        join_by = "cell_id",
        overwrite = TRUE
    )
    list(obj = obj, snp_ids = snp_ids)
}

test_that("add_molecule_phase() errors when no XCI diagnostics are stored", {
    ref <- Matrix::Matrix(matrix(1L, 2, 2), sparse = TRUE)
    alt <- Matrix::Matrix(matrix(1L, 2, 2), sparse = TRUE)
    snp_info <- data.frame(chrom = "chrX", pos = c(1L, 2L), ref = "A", alt = "G")
    barcode_info <- data.frame(barcode = c("c1", "c2"), donor = c("donor0", "donor0"))
    obj <- SNPData(ref_count = ref, alt_count = alt, snp_info = snp_info, barcode_info = barcode_info)

    # Verify the function refuses to run before assign_xci() has stored a fit
    expect_error(add_molecule_phase(obj, bam_files = c(donor0 = "x.bam")), "Run assign_xci")
})

test_that("add_molecule_phase() errors on unrecognised donor names in bam_files", {
    fixture <- make_phase_fixture()

    # Verify a BAM keyed to a donor absent from barcode_info is rejected up front
    expect_error(
        add_molecule_phase(fixture$obj, bam_files = c(not_a_donor = "x.bam")),
        "not found in barcode_info\\$donor"
    )
})

test_that("add_molecule_phase() ignores 'doublet' and 'unassigned' entries in bam_files", {
    # A fixture where "doublet" is a real, valid donor label in barcode_info
    # (as Vireo emits it) -- so the unknown_donors check alone wouldn't catch
    # it, and the exclusion has to be deliberate, not incidental.
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
    snp_id <- get_snp_info(obj)$snp_id
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
        data.frame(cell_id = get_barcode_info(obj)$cell_id, active_x = c("X2", "X2", "X1", "X1", NA, NA)),
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
    called_for <- character()
    testthat::local_mocked_bindings(
        extract_snp_calls = function(bam_file, snp_info, barcodes = NULL, ...) {
            called_for <<- c(called_for, bam_file)
            list(tallies = fake_tallies, reads = tibble::tibble())
        },
        .package = "snplet"
    )

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

    add_molecule_phase(obj, bam_files = c(donor0 = "fake.bam", doublet = "fake.bam", unassigned = "fake.bam"))

    # Verify the log names the excluded non-donor labels
    log_lines <- readLines(log_file, warn = FALSE)
    expect_true(any(grepl("doublet, unassigned", log_lines)))
    # Check the BAM extraction step itself was only ever invoked once -- for
    # the real donor, never for "doublet" or "unassigned"
    expect_length(called_for, 1)
})

test_that("add_molecule_phase() returns x unchanged when only doublet/unassigned are given", {
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
    snp_id <- get_snp_info(obj)$snp_id
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
    result <- add_molecule_phase(obj, bam_files = c(doublet = "fake.bam", unassigned = "fake.bam"))
    expect_identical(get_donor_snp_info(result), get_donor_snp_info(obj))
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
    testthat::local_mocked_bindings(
        extract_snp_calls = function(...) list(tallies = fake_tallies, reads = tibble::tibble()),
        .package = "snplet"
    )

    result <- add_molecule_phase(fixture$obj, bam_files = c(donor0 = "fake.bam"))
    donor_snp_info <- get_donor_snp_info(result)

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

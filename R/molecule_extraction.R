# ==============================================================================
# Molecule-level allele counting from long-read BAMs
#
# cellSNP-lite counts *reads* at *SNPs*, which over-counts a gene's expression
# whenever a single molecule spans several of its heterozygous SNPs: each SNP
# then casts an independent allele vote for what is really one transcript.
# The functions here count molecules instead. BAM reads already carry a
# corrected cell barcode (CB) and UMI (UB) tag, so a molecule's identity needs
# no inference: (CB, UB) is the counting unit, and every SNP it covers is an
# independent observation of the same haplotype. Voting across those SNPs
# both removes the double-counting and yields read-backed phase for free.
# ==============================================================================

#' Extract the aligned base of every read at every target SNP
#'
#' Reads the region(s) spanned by `snp` from an indexed BAM and, for every
#' alignment overlapping a target position, records the base and quality
#' called there together with the read's identity, the information needed
#' to later group calls into molecules and phase blocks. Equivalent in
#' principle to `GenomicAlignments::pileLettersAt()`, but retains the read of
#' origin, which that function discards.
#'
#' A SNP falling in a deletion or an intron of a given read simply yields no
#' row for that read: a genuine no-call rather than a coerced base. Only
#' `M`/`=`/`X` CIGAR operations are treated as aligned to the reference.
#'
#' @param bam_file Character scalar, required. Path to an indexed BAM with
#'   `CB` and `UB` tags.
#' @param snp_info A data.frame/tibble, required, with one row per target SNP
#'   and columns `snp_id`, `chrom`, `pos`, `ref` (reference allele), and `alt`
#'   (alternate allele). Restrict to biallelic SNVs before calling: there is
#'   no unambiguous REF/ALT base to compare an aligned base against for a
#'   multi-allelic site or an indel.
#' @param barcodes Character vector, optional (default `NULL`, keeping every
#'   barcode carrying a `CB` and `UB` tag). Cell barcodes to retain; typically
#'   the barcodes of one donor, since a pooled sample dilutes the allele
#'   fraction of any site heterozygous in only some donors.
#' @param min_mapq Integer (default 20). Minimum mapping quality.
#' @param min_baseq Integer (default 10). Minimum base quality at the SNP position.
#' @param chunk_size Integer (default 20000). Alignments processed per chunk,
#'   to bound memory when expression is uneven across the region.
#' @param merge_gap Integer (default 10000). Alignments are fetched over SNP
#'   positions merged into windows within this many bases of each other,
#'   rather than over each SNP individually, since a `which` with many ranges
#'   would otherwise return a read once per overlapping range, duplicating
#'   any read spanning several nearby SNPs.
#' @param threads Integer (default 4). Threads for the `samtools` pre-filter
#'   (see `prefilter_bam` below); gains plateau at 4.
#'
#' @return A list of two tibbles:
#'   \describe{
#'     \item{tallies}{One row per (`barcode`, `umi`, `snp_id`, `allele`)
#'       combination, with `n_calls` the number of reads of that molecule
#'       agreeing on that allele ("REF", "ALT", or "OTH") at that SNP.}
#'     \item{reads}{One row per (`barcode`, `umi`, `qname`), the distinct
#'       reads behind each molecule, with each read's alignment `strand`
#'       (`"+"`/`"-"`), for `molecule_read_strand()` to resolve into a single
#'       strand per molecule.}
#'   }
#'
#' @family molecule-level allele counting functions
#' @export
extract_snp_calls <- function(
    bam_file,
    snp_info,
    barcodes = NULL,
    min_mapq = 20L,
    min_baseq = 10L,
    chunk_size = 20000L,
    merge_gap = 10000L,
    threads = 4L
) {
    required_cols <- c("snp_id", "chrom", "pos", "ref", "alt")
    missing_cols <- setdiff(required_cols, colnames(snp_info))
    if (length(missing_cols) > 0) {
        stop("snp_info is missing required column(s): ", paste(missing_cols, collapse = ", "))
    }
    if (nrow(snp_info) == 0) {
        stop("snp_info has no rows")
    }

    snp_gr <- plyranges::as_granges(snp_info, seqnames = chrom, start = pos, end = pos)

    # GenomicRanges::reduce must be qualified: purrr (pulled in by furrr) also
    # exports reduce(), and under future the environment of an exported function
    # global is rebound, which changes how an unqualified call resolves.
    windows <- GenomicRanges::reduce(snp_gr, min.gapwidth = merge_gap)
    logger::log_info("fetching reads over {length(windows)} merged window(s)")

    what <- c("qname", "seq", "qual", "mapq")
    # Secondary and supplementary alignments would multiply-count a single
    # molecule, which is exactly what this function exists to prevent.
    bam_flag <- Rsamtools::scanBamFlag(
        isSecondaryAlignment = FALSE,
        isSupplementaryAlignment = FALSE,
        isUnmappedQuery = FALSE,
        isDuplicate = FALSE
    )

    prefiltered <- .prefilter_bam(bam_file, windows, barcodes, min_mapq, threads)
    if (is.null(prefiltered)) {
        param <- Rsamtools::ScanBamParam(which = windows, what = what, tag = c("CB", "UB"), flag = bam_flag)
        galn <- GenomicAlignments::readGAlignments(bam_file, param = param, use.names = FALSE)
    } else {
        # The reads are already restricted to the windows, so no `which` is
        # needed and the intermediate BAM needs no index.
        on.exit(unlink(prefiltered), add = TRUE)
        param <- Rsamtools::ScanBamParam(what = what, tag = c("CB", "UB"), flag = bam_flag)
        galn <- GenomicAlignments::readGAlignments(prefiltered, param = param, use.names = FALSE)
    }

    # A read overlapping two windows is returned twice; keep one copy.
    galn <- galn[!duplicated(S4Vectors::mcols(galn)$qname)]

    keep <- !is.na(S4Vectors::mcols(galn)$CB) &
        !is.na(S4Vectors::mcols(galn)$UB) &
        S4Vectors::mcols(galn)$mapq >= min_mapq
    if (!is.null(barcodes)) {
        keep <- keep & S4Vectors::mcols(galn)$CB %in% barcodes
    }
    galn <- galn[keep]
    if (length(galn) == 0) {
        stop("No alignments passed filtering")
    }

    # findOverlaps() requires the position seqlevels to match the BAM's exactly.
    Seqinfo::seqlevels(snp_gr) <- Seqinfo::seqlevels(galn)

    chunks <- split(seq_along(galn), ceiling(seq_along(galn) / chunk_size))
    logger::log_info("processing {length(galn)} alignments in {length(chunks)} chunk(s)")

    per_chunk <- purrr::map(chunks, function(chunk_idx) {
        chunk_galn <- galn[chunk_idx]
        calls <- .map_bases_to_reads(chunk_galn, snp_gr, snp_info$snp_id) %>%
            dplyr::filter(base_quality >= min_baseq) %>%
            dplyr::mutate(
                barcode = S4Vectors::mcols(chunk_galn)$CB[read_idx],
                umi = S4Vectors::mcols(chunk_galn)$UB[read_idx],
                qname = S4Vectors::mcols(chunk_galn)$qname[read_idx],
                strand = as.character(BiocGenerics::strand(chunk_galn))[read_idx],
                allele = dplyr::case_when(
                    base == snp_info$ref[snp_idx] ~ "REF",
                    base == snp_info$alt[snp_idx] ~ "ALT",
                    TRUE ~ "OTH"
                )
            )

        # Collapsed before returning, so the (read x SNP) product for the whole
        # region never exists at once.
        list(
            tallies = dplyr::count(calls, barcode, umi, snp_id, allele, name = "n_calls"),
            reads = dplyr::distinct(calls, barcode, umi, qname, strand)
        )
    })

    list(
        tallies = dplyr::bind_rows(purrr::map(per_chunk, "tallies")) %>%
            dplyr::summarise(n_calls = sum(n_calls), .by = c(barcode, umi, snp_id, allele)),
        reads = dplyr::bind_rows(purrr::map(per_chunk, "reads")) %>%
            dplyr::distinct()
    )
}

#' Map each SNP position into query coordinates for every overlapping read
#'
#' Inlines the logic of `GenomicAlignments:::.pileLettersOnSingleRefAt` so the
#' read of origin survives, precisely what is needed to group calls into
#' molecules.
#'
#' @param galn A GAlignments object.
#' @param snp_gr A GRanges of target SNPs, width 1, same seqlevels as `galn`.
#' @param snp_ids Character vector of SNP identifiers, same length and order
#'   as `snp_gr`.
#'
#' @return A tibble with columns `snp_idx` (position in `snp_gr`/`snp_ids`),
#'   `read_idx` (position in `galn`), `snp_id`, `base` (the aligned base), and
#'   `base_quality` (Phred-scaled integer).
#'
#' @keywords internal
.map_bases_to_reads <- function(galn, snp_gr, snp_ids) {
    # cigarRangesAlong*Space() are formally deprecated in favour of the new
    # cigarillo package (GenomicAlignments >= 1.45.5) but still fully
    # functional; suppressed rather than migrated so behaviour doesn't shift
    # underneath a released Bioconductor API before cigarillo has matured.
    ref_ranges <- withCallingHandlers(
        GenomicAlignments::cigarRangesAlongReferenceSpace(
            GenomicAlignments::cigar(galn),
            pos = BiocGenerics::start(galn),
            ops = c("M", "=", "X")
        ),
        deprecatedWarning = function(w) invokeRestart("muffleWarning")
    )
    query_ranges <- withCallingHandlers(
        GenomicAlignments::cigarRangesAlongQuerySpace(
            GenomicAlignments::cigar(galn),
            ops = c("M", "=", "X")
        ),
        deprecatedWarning = function(w) invokeRestart("muffleWarning")
    )

    read_of_range <- IRanges::togroup(IRanges::PartitioningByWidth(ref_ranges))
    unlisted_ref <- unlist(ref_ranges, use.names = FALSE)
    unlisted_query <- unlist(query_ranges, use.names = FALSE)
    query_shift <- BiocGenerics::start(unlisted_ref) - BiocGenerics::start(unlisted_query)

    # Overlaps are computed per chromosome to keep the shift vectors aligned.
    galn_chrom <- as.character(Seqinfo::seqnames(galn))
    snp_chrom <- as.character(Seqinfo::seqnames(snp_gr))
    snp_pos <- BiocGenerics::start(snp_gr)

    per_chrom_calls <- purrr::map(unique(snp_chrom), function(chrom) {
        snp_idx_on_chrom <- which(snp_chrom == chrom)
        is_range_on_chrom <- galn_chrom[read_of_range] == chrom

        # One hit per (SNP, aligned range) overlap. Both index vectors are
        # local to the subsets passed in, so each is mapped back to its
        # position in the full snp_gr/unlisted_ref before being used.
        hits <- IRanges::findOverlaps(
            IRanges::IRanges(snp_pos[snp_idx_on_chrom], width = 1L),
            unlisted_ref[is_range_on_chrom]
        )
        range_idx <- which(is_range_on_chrom)[S4Vectors::subjectHits(hits)]
        snp_idx <- snp_idx_on_chrom[S4Vectors::queryHits(hits)]
        read_idx <- read_of_range[range_idx]
        # The SNP's reference position, shifted into the read's own coordinates
        # by the insertion/deletion offset accumulated up to its aligned block.
        query_pos <- snp_pos[snp_idx] - query_shift[range_idx]

        base_quality_chars <- as.character(
            Biostrings::subseq(S4Vectors::mcols(galn)$qual[read_idx], start = query_pos, width = 1L)
        )

        tibble::tibble(
            snp_idx = snp_idx,
            read_idx = read_idx,
            snp_id = snp_ids[snp_idx],
            base = as.character(
                Biostrings::subseq(S4Vectors::mcols(galn)$seq[read_idx], start = query_pos, width = 1L)
            ),
            # charToRaw() is not vectorised, so the single-character qualities are
            # concatenated and decoded in one pass.
            base_quality = as.integer(charToRaw(paste0(base_quality_chars, collapse = ""))) - 33L
        )
    })

    dplyr::bind_rows(per_chrom_calls)
}

#' Cut a BAM down to the reads that matter, using threaded samtools
#'
#' Rsamtools decompresses every read in a window before R can discard it;
#' `samtools` can apply the same barcode/UMI/quality criteria in threaded C
#' first, which is where nearly all of the speed comes from. This is a pure
#' accelerator: the filters applied afterwards in R remain the definition of
#' what is kept, so the two paths cannot diverge.
#'
#' @param bam_file Path to an indexed BAM.
#' @param windows A GRanges of merged fetch windows.
#' @param barcodes Character vector of cell barcodes to retain, or `NULL`.
#' @param min_mapq Minimum mapping quality.
#' @param threads Threads passed to `samtools view -@`.
#'
#' @return Path of an uncompressed temporary BAM, or `NULL` if `samtools` is
#'   unavailable, in which case the calling function falls back to reading
#'   the original BAM directly.
#'
#' @keywords internal
.prefilter_bam <- function(bam_file, windows, barcodes, min_mapq, threads) {
    if (!nzchar(Sys.which("samtools"))) {
        return(NULL)
    }

    bed_file <- tempfile(fileext = ".bed")
    utils::write.table(
        data.frame(
            as.character(Seqinfo::seqnames(windows)),
            BiocGenerics::start(windows) - 1L,
            BiocGenerics::end(windows)
        ),
        bed_file,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE,
        col.names = FALSE
    )
    out_bam <- tempfile(fileext = ".bam")
    on.exit(unlink(bed_file), add = TRUE)

    # 0xD04 = unmapped + duplicate + secondary + supplementary.
    view_args <- c("view", "-@", threads, "-M", "-L", bed_file, "-q", min_mapq, "-F", "0xD04")
    if (is.null(barcodes)) {
        view_args <- c(view_args, "-e", shQuote("[CB] && [UB]"))
    } else {
        barcode_file <- tempfile(fileext = ".txt")
        writeLines(barcodes, barcode_file)
        on.exit(unlink(barcode_file), add = TRUE)
        # samtools permits only one -d/-D tag filter, so a second tag requirement
        # has to go through the expression filter -e instead.
        view_args <- c(view_args, "-D", paste0("CB:", barcode_file), "-e", shQuote("[UB]"))
    }

    exit_status <- system2("samtools", c(view_args, "-u", "-o", out_bam, bam_file))
    if (exit_status != 0) {
        unlink(out_bam)
        warning("samtools pre-filter failed; falling back to reading the full BAM")
        return(NULL)
    }
    out_bam
}


#' Resolve the allele call for each (molecule, SNP)
#'
#' Duplicate reads of one molecule vote on the allele at each SNP. Only REF
#' and ALT are retained; OTH is a sequencing error at a known biallelic site
#' and carries no haplotype information, so it is excluded before the vote
#' rather than after -- a molecule read as 3 OTH and 1 REF still yields its
#' REF call instead of being discarded for having no majority allele.
#'
#' A molecule whose reads tie between REF and ALT has no majority to read off
#' and is dropped, matching how \code{\link{haplotype_expression_by_molecule}}
#' and \code{\link{molecule_haplotype_counts}} treat a tied haplotype vote.
#' Note that a BAM whose reads are already UMI-collapsed gives one read per
#' molecule, so neither case arises: the vote that matters there is across the
#' several SNPs a molecule spans, which those two functions take.
#'
#' @param tallies A tibble, required, as returned by
#'   `extract_snp_calls()$tallies`, with columns `barcode`, `umi`, `snp_id`,
#'   `allele`, and `n_calls`.
#'
#' @return A tibble with one row per (`barcode`, `umi`, `snp_id`) and columns
#'   `barcode`, `umi`, `snp_id`, `allele` (the majority call, "REF" or "ALT"),
#'   and `n_calls` (reads backing that call). A (molecule, SNP) with no REF or
#'   ALT read, or whose reads tie between the two, is absent from the result.
#'
#' @family molecule-level allele counting functions
#' @export
molecule_snp_alleles <- function(tallies) {
    # Excluded before the vote: OTH is never a candidate to win, so filtering
    # afterwards would let it take the argmax and discard a real REF/ALT call.
    informative <- dplyr::filter(tallies, allele %in% c("REF", "ALT"))
    # max() of a zero-row group warns before returning -Inf; short-circuit instead.
    if (nrow(informative) == 0) {
        return(dplyr::select(informative, barcode, umi, snp_id, allele, n_calls))
    }

    informative %>%
        # n_top must be computed before n_calls is redefined, or it would see
        # the scalar maximum instead of the original per-row counts.
        dplyr::summarise(
            n_top = sum(n_calls == max(n_calls)),
            allele = allele[which.max(n_calls)],
            n_calls = max(n_calls),
            .by = c(barcode, umi, snp_id)
        ) %>%
        # A tie has no majority to read off, so it is dropped, not resolved by
        # row order -- see molecule_read_strand()'s own tied vote for the same rule.
        dplyr::filter(n_top == 1) %>%
        dplyr::select(barcode, umi, snp_id, allele, n_calls)
}

#' Resolve the alignment strand of each molecule
#'
#' A molecule is one transcript, so every read behind it should agree on
#' alignment strand; a majority vote absorbs the rare mismapped or chimeric
#' read rather than letting one read decide. This is alignment strand only;
#' converting it to the strand of the original transcript (needed to
#' disambiguate a SNP overlapping genes on opposite strands) additionally
#' requires the BAM's sense/antisense orientation, since some demultiplexing
#' pipelines flip reads relative to the transcript (see
#' `.infer_bam_strand_orientation()`).
#'
#' A molecule whose reads tie between the two strands has no majority to read
#' off and reports `NA`, matching how \code{\link{molecule_snp_alleles}}
#' treats a tied allele vote. A tie means the evidence does not say which
#' strand the transcript came from, and a strand picked between two equal
#' options would be used downstream exactly as a resolved one:
#' `haplotype_expression_by_molecule()` attributes an ambiguous SNP's molecule
#' to whichever overlapping gene shares its strand, so a guess there sends the
#' molecule's counts to a gene it may not have come from, while `NA` correctly
#' withholds it. Note that a UMI-collapsed BAM gives one read per molecule and
#' so cannot tie.
#'
#' @param reads A tibble, required, as returned by `extract_snp_calls()$reads`,
#'   with columns `barcode`, `umi`, `qname`, `strand`.
#'
#' @return A tibble with one row per (`barcode`, `umi`) and columns
#'   `barcode`, `umi`, `strand` (the majority call: `"+"`, `"-"`, or `NA`
#'   where the molecule's reads tie between the two).
#'
#' @family molecule-level allele counting functions
#' @export
molecule_read_strand <- function(reads) {
    reads %>%
        dplyr::count(barcode, umi, strand, name = "n_reads") %>%
        # with_ties = TRUE keeps a tie as two rows so it can be detected below,
        # rather than resolved arbitrarily by row order.
        dplyr::slice_max(n_reads, n = 1, by = c(barcode, umi), with_ties = TRUE) %>%
        # A tie has no majority to read off, so it reports NA -- the same rule
        # molecule_snp_alleles() applies to its own tied vote.
        dplyr::summarise(
            strand = dplyr::if_else(dplyr::n() == 1L, strand[1], NA_character_),
            .by = c(barcode, umi)
        )
}

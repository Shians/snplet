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
    # Excluded before the vote, not after: OTH is a sequencing error at a known
    # biallelic site, so it is not a candidate to win. Filtering afterwards
    # would let it take the argmax and discard the molecule's REF/ALT reads
    # along with it, losing a call that was there to be made.
    informative <- dplyr::filter(tallies, allele %in% c("REF", "ALT"))
    # summarise() still evaluates its expressions against a zero-row input, and
    # max() of nothing warns before returning -Inf. The result is empty either
    # way, so the empty case returns early rather than emitting a warning that
    # says nothing about the data.
    if (nrow(informative) == 0) {
        return(dplyr::select(informative, barcode, umi, snp_id, allele, n_calls))
    }

    informative %>%
        # n_top is computed before n_calls is redefined: summarise() evaluates
        # its arguments in order and each one masks the column it names, so a
        # later reference to n_calls would see the scalar maximum and count
        # every observation as untied.
        dplyr::summarise(
            n_top = sum(n_calls == max(n_calls)),
            allele = allele[which.max(n_calls)],
            n_calls = max(n_calls),
            .by = c(barcode, umi, snp_id)
        ) %>%
        # A tied molecule is residual noise with no majority to read off, so it
        # is dropped rather than resolved by row order -- the rule
        # haplotype_expression_by_molecule() and molecule_haplotype_counts()
        # already apply to their own votes.
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
        # with_ties = TRUE so a tie survives as two rows to be recognised
        # below; dropping one arbitrarily here would resolve the molecule by
        # row order and there would be nothing left to detect.
        dplyr::slice_max(n_reads, n = 1, by = c(barcode, umi), with_ties = TRUE) %>%
        # A molecule whose reads split evenly between strands has no majority
        # to read off, so it reports NA rather than whichever strand sorted
        # first -- the rule molecule_snp_alleles() and .pool_donor_calls()
        # already apply to their own tied votes. NA is what the consumers
        # expect for "unresolvable": .molecule_gene_phase_calls() requires a
        # non-NA strand before attributing an ambiguous SNP's molecule to a
        # gene, so a tie excludes the molecule there instead of sending its
        # counts to an arbitrary one of two overlapping genes.
        dplyr::summarise(
            strand = dplyr::if_else(dplyr::n() == 1L, strand[1], NA_character_),
            .by = c(barcode, umi)
        )
}

#' Phase heterozygous SNPs directly from the molecules that span them
#'
#' Two het SNPs observed on the same molecule sit on the same physical
#' haplotype, so their alleles co-occur non-randomly: mostly REF/REF and
#' ALT/ALT if the two REF alleles share a haplotype ("same"), or REF/ALT and
#' ALT/REF if they do not ("opposite"). This is read-backed phasing, and long
#' reads supply it directly, without needing a statistical fit.
#'
#' Orientation propagates across accepted edges, flipping on "opposite"
#' links; each connected component becomes a phase block. Blocks cannot be
#' joined beyond the reach of a single molecule, so a fragmented result is
#' expected rather than a failure. Phasing runs over every SNP regardless of
#' gene assignment: whether two SNPs share a haplotype is a property of the
#' molecule, not of gene annotation, and an ambiguous (multi-gene) SNP can
#' still serve as a valid bridge between two unambiguous ones.
#'
#' @section Accepting an edge:
#' Each molecule spanning two het SNPs is one vote on whether their REF
#' alleles share a haplotype. Writing \code{e} for `error_rate`, the two
#' hypotheses predict opposite things of a molecule: under "same" it agrees
#' with probability \code{1 - e}, under "opposite" with probability \code{e}.
#' With \code{k} of \code{n} molecules agreeing, the log-likelihood ratio
#' between them reduces to
#'
#' \deqn{LLR = |2k - n| \cdot \log\left(\frac{1-e}{e}\right)}
#'
#' the binomial coefficient cancelling because it depends on the data alone.
#' The two factors read separately: \code{|2k - n|} is the \emph{margin},
#' agreements minus disagreements, and \code{log((1-e)/e)} is what one
#' molecule is worth in log-odds, fixed by \code{e}. Evidence is therefore
#' additive per molecule -- each agreeing molecule adds that weight, each
#' disagreeing one subtracts it -- and an edge is accepted once the total
#' reaches `min_llr`.
#'
#' This replaces the fraction cutoff used previously, which conflated how
#' clean an edge is with how much of it there is. A fraction is also coarse
#' at the small \code{n} most edges have: at `min_molecules = 5` a 0.9 cutoff
#' admits only a unanimous 5 of 5, rejecting 4 of 5 despite it carrying
#' \code{LLR = 8.8} at \code{e = 0.05}. The practical effect of the change is
#' sensitivity rather than accuracy: far more true edges are accepted, which
#' leaves blocks less fragmented, at a small and now-detectable rise in
#' wrongly oriented ones (see \sQuote{Internally inconsistent blocks}).
#'
#' The LLR treats every molecule as an independent observation, which holds
#' across cells but not within one. A cell contributes many transcripts of the
#' same two chromosome copies, and on the X the inactive copy is largely
#' silent, so a gene's molecules there are near-monoallelic: the SNPs agree
#' because one haplotype was sampled repeatedly, not because both were seen to
#' co-occur. Ambient RNA and an undetected doublet are likewise properties of
#' a barcode, flipping that cell's molecules together rather than one at a
#' time. Molecules from a single cell therefore inflate \code{n} without
#' adding proportionate evidence, and `min_cells` requires an edge to be
#' corroborated across barcodes before it is accepted. It is a floor on
#' independence, not a second likelihood: the LLR itself still counts
#' molecules.
#'
#' The count is taken over the cells backing the relation the edge is accepted
#' on, not over every cell spanning the pair. The two differ precisely when a
#' cell dissents, and it is that case the distinction matters for: four
#' molecules from one cell reading "same" alongside one from another reading
#' "opposite" spans two barcodes but rests, as far as "same" is concerned, on
#' a single cell. Counting the dissenter towards the relation it contradicts
#' would let it vouch for exactly the single-cell edge `min_cells` exists to
#' reject.
#'
#' Note that \code{e} is treated as independent across molecules, which
#' sequencing error is and contamination is not: ambient RNA is correlated
#' within a cell and a doublet flips many molecules together. Where
#' contamination rather than base-calling error dominates, the LLR is
#' optimistic, so set \code{e} from observed discordance rather than from a
#' sequencing-error prior, and treat `min_llr` as a threshold to calibrate
#' rather than a p-value.
#'
#' @section Internally inconsistent blocks:
#' Orientations are fixed by a spanning tree of each component, so any edge
#' that closes a cycle is not needed to phase its endpoints -- but it is an
#' independent prediction of the relation between them. Where such an edge
#' contradicts the orientations already assigned, no assignment of alleles to
#' two haplotypes can satisfy every accepted edge at once. That is impossible
#' for a diploid genome, so at least one edge in the block is wrong: most
#' often a spurious link from ambient RNA, a mismapped paralogue, or an
#' undetected doublet whose two genotypes are read as one.
#'
#' Such blocks are still returned, oriented by the spanning tree as before --
#' the contradiction says one edge is wrong, not which one, and dropping the
#' block would discard its majority of sound edges along with the bad one.
#' They are instead flagged \code{block_conflict = TRUE} with
#' \code{n_block_conflicts} giving the number of contradicting edges, so a
#' block built on an unreliable link can be inspected or excluded rather than
#' trusted silently. A block whose SNPs are all genuinely linked has no
#' conflicts at all, so any non-zero count is worth investigating.
#'
#' @param per_snp A tibble, required, as returned by `molecule_snp_alleles()`,
#'   with columns `barcode`, `umi`, `snp_id`, `allele`.
#' @param min_molecules Integer (default 5). Molecules required to accept a
#'   SNP pair as an edge, applied before `min_llr` so a large ratio resting on
#'   one or two molecules cannot qualify.
#' @param min_cells Integer, `>= 1` (default 2). Distinct cells whose molecules
#'   must back the relation an edge is accepted on, applied alongside
#'   `min_molecules`. Counted over the cells voting for that relation, not over
#'   every cell spanning the pair, so a cell whose molecules argue for the
#'   losing relation does not help its rival meet this floor. Molecules from
#'   one cell are not independent evidence; see \sQuote{Accepting an edge}. Set
#'   to 1 to count molecules only, restoring the previous behaviour.
#' @param error_rate Numeric, in `(0, 0.5)` (default 0.05). Probability that a
#'   single molecule reports the wrong relation, absorbing sequencing error,
#'   mismapping, ambient RNA and undetected doublets. Sets how much one
#'   molecule is worth as evidence; see \sQuote{Accepting an edge}.
#' @param min_llr Numeric, `>= 0` (default 3). Log-likelihood ratio an edge
#'   must reach to be accepted. Roughly, 3 is ~20:1 odds and 4.6 is ~100:1.
#' @param min_consistency Deprecated and ignored (default `NULL`). Edges are
#'   accepted on a likelihood ratio rather than a fraction; passing it warns
#'   and has no effect. Use `error_rate` and `min_llr`.
#'
#' @return A tibble with columns `snp_id`, `block` (integer phase-block id,
#'   unique within this call), `allele_on_h1` ("REF" or "ALT", the
#'   allele carried by haplotype 1 at that SNP; H1 is an arbitrary label local
#'   to each block, not oriented to X1/X2), `n_block_conflicts` (integer;
#'   edges in this SNP's block that contradict the orientation assigned to
#'   them, 0 for a consistent block), and `block_conflict` (logical;
#'   `n_block_conflicts > 0`). See \sQuote{Internally inconsistent blocks}.
#'   SNPs that could not be linked to any other SNP by an accepted edge are
#'   absent from the result.
#'
#' @family molecule-level allele counting functions
#' @export
phase_snps <- function(
    per_snp,
    min_molecules = 5L,
    min_cells = 2L,
    error_rate = 0.05,
    min_llr = 3,
    min_consistency = NULL
) {
    # Retained only so an existing call still runs; it no longer has any effect.
    if (!is.null(min_consistency)) {
        logger::log_warn(
            "phase_snps(min_consistency = ) is deprecated and ignored; edges are now accepted on a ",
            "likelihood ratio instead of a fraction. Use error_rate and min_llr (see ?phase_snps)."
        )
    }
    # error_rate is a per-molecule probability of a discordant observation, so
    # it must leave room for one: at 0 no disagreement is ever explicable and
    # the weight is infinite, at 0.5 the two hypotheses predict identical data
    # and the weight is 0, and beyond 0.5 the test reads backwards.
    if (!is.numeric(error_rate) || length(error_rate) != 1 || is.na(error_rate)) {
        stop("error_rate must be a single non-missing number.")
    }
    if (error_rate <= 0 || error_rate >= 0.5) {
        stop("error_rate must be in (0, 0.5); got ", error_rate, ".")
    }
    if (!is.numeric(min_llr) || length(min_llr) != 1 || is.na(min_llr) || min_llr < 0) {
        stop("min_llr must be a single non-missing number >= 0.")
    }
    if (!is.numeric(min_cells) || length(min_cells) != 1 || is.na(min_cells) || min_cells < 1) {
        stop("min_cells must be a single non-missing number >= 1.")
    }
    molecule_snps <- per_snp %>%
        dplyr::arrange(barcode, umi, snp_id) %>%
        dplyr::summarise(snps = list(snp_id), alleles = list(allele), .by = c(barcode, umi)) %>%
        dplyr::filter(lengths(snps) >= 2)

    if (nrow(molecule_snps) == 0) {
        return(tibble::tibble(
            snp_id = character(),
            block = integer(),
            allele_on_h1 = character(),
            n_block_conflicts = integer(),
            block_conflict = logical()
        ))
    }

    # Every molecule contributes one vote per pair of SNPs it spans: the two
    # alleles it read either agree ("same" -- both REF or both ALT, so the two
    # REF alleles share a haplotype) or they do not. `snp_pairs` has the two
    # rows of combn() naming the first and second member of each pair, so
    # first_of_pair/second_of_pair index the molecule's own SNPs and alleles in
    # lockstep.
    pair_votes <- purrr::map(seq_len(nrow(molecule_snps)), function(molecule_idx) {
        molecule_snp_ids <- molecule_snps$snps[[molecule_idx]]
        molecule_alleles <- molecule_snps$alleles[[molecule_idx]]
        snp_pairs <- utils::combn(length(molecule_snp_ids), 2)
        first_of_pair <- snp_pairs[1, ]
        second_of_pair <- snp_pairs[2, ]
        tibble::tibble(
            snp_a = molecule_snp_ids[first_of_pair],
            snp_b = molecule_snp_ids[second_of_pair],
            same = molecule_alleles[first_of_pair] == molecule_alleles[second_of_pair],
            # Carried through so an edge's support can be counted in cells as
            # well as molecules; see `min_cells`.
            barcode = molecule_snps$barcode[molecule_idx]
        )
    })

    # Evidence for an edge is the *margin* between agreeing and disagreeing
    # molecules, not the fraction agreeing: log(L_same / L_opposite) reduces to
    # (2 * n_same - n) * log((1 - error_rate) / error_rate), a net vote count
    # times a fixed per-molecule weight. A fraction conflates how clean an edge
    # is with how much of it there is, and at small n it is the coarser of the
    # two -- 4 of 5 agreeing is 0.8, below any useful fraction cutoff, while
    # carrying LLR 8.8 at error_rate = 0.05, which is decisive.
    weight_per_molecule <- log((1 - error_rate) / error_rate)
    edges <- dplyr::bind_rows(pair_votes) %>%
        dplyr::summarise(
            n = dplyr::n(),
            n_same = sum(same),
            # Counted per relation rather than over the pair as a whole: a cell
            # is only independent evidence *for* the relation its molecules
            # actually voted for. Counting every cell touching the pair would
            # let a dissenting cell satisfy `min_cells` on behalf of the
            # relation it argues against -- 4 molecules from one cell saying
            # "same" plus 1 from another saying "opposite" would read as
            # two-cell support for "same", which is exactly the single-cell
            # edge the threshold exists to reject.
            n_cells_same = dplyr::n_distinct(barcode[same]),
            n_cells_opposite = dplyr::n_distinct(barcode[!same]),
            .by = c(snp_a, snp_b)
        ) %>%
        dplyr::mutate(
            relation = dplyr::if_else(n_same >= n - n_same, "same", "opposite"),
            consistency = pmax(n_same, n - n_same) / n,
            # Cells backing the relation that won, so `min_cells` is a floor on
            # the independence of the evidence actually being accepted.
            n_cells = dplyr::if_else(relation == "same", n_cells_same, n_cells_opposite),
            # Sign records which relation is favoured, which `relation` already
            # holds; only the strength of the evidence is thresholded.
            llr = abs(2 * n_same - n) * weight_per_molecule
        ) %>%
        dplyr::filter(n >= min_molecules, n_cells >= min_cells, llr >= min_llr)

    if (nrow(edges) == 0) {
        return(tibble::tibble(
            snp_id = character(),
            block = integer(),
            allele_on_h1 = character(),
            n_block_conflicts = integer(),
            block_conflict = logical()
        ))
    }

    snp_ids <- sort(unique(c(edges$snp_a, edges$snp_b)))
    # orientation[snp] is 0 when that SNP's REF allele sits on H1 and 1 when it
    # sits on H2; it doubles as the "already visited" marker, NA meaning the
    # traversal has not reached this SNP yet.
    orientation <- stats::setNames(rep(NA_integer_, length(snp_ids)), snp_ids)
    block <- stats::setNames(rep(NA_integer_, length(snp_ids)), snp_ids)

    # Undirected adjacency, keyed by SNP: every edge appears twice, once walkable
    # from each end, so `neighbours_of_snp[[snp]]` lists everything reachable in
    # one step whichever end the traversal arrives from. `edge_idx` carries each
    # half-edge back to its row in `edges`, so a contradiction found while
    # walking can be attributed to the edge that caused it.
    neighbours_of_snp <- split(
        rbind(
            data.frame(to = edges$snp_b, flip = edges$relation == "opposite", edge_idx = seq_len(nrow(edges))),
            data.frame(to = edges$snp_a, flip = edges$relation == "opposite", edge_idx = seq_len(nrow(edges)))
        ),
        c(edges$snp_a, edges$snp_b)
    )

    # An edge reaching an already-oriented SNP is a cycle closing. The spanning
    # tree has already fixed both endpoints' orientations, so this edge is not
    # needed to phase anything -- but it is an independent prediction of the
    # relation between them, and it either agrees with the tree or it does not.
    # A disagreement means no assignment of alleles to two haplotypes can
    # satisfy every accepted edge at once, which is physically impossible for a
    # diploid genome and so evidence that one of the edges is wrong (a spurious
    # link from ambient RNA, a mismapped paralogue, or a doublet's two
    # genotypes read as one). Silently keeping the tree's answer would discard
    # exactly the signal that says the block is untrustworthy, so the conflicts
    # are counted and reported per block instead.
    is_edge_inconsistent <- rep(FALSE, nrow(edges))
    conflicts_per_block <- integer(0)

    # Breadth-first traversal of the edge graph. Each unvisited SNP seeds a new
    # connected component -- a phase block -- and is arbitrarily declared to
    # carry its REF allele on H1 (orientation 0); every SNP reachable from it
    # then inherits an orientation forced by the relations along the way. The
    # visited SNPs of one component are exactly one block, so the traversal
    # assigns block membership and relative phase in a single pass.
    block_id <- 0L
    for (seed_snp in snp_ids) {
        if (!is.na(orientation[seed_snp])) {
            next
        }
        block_id <- block_id + 1L
        orientation[seed_snp] <- 0L
        block[seed_snp] <- block_id
        snp_queue <- seed_snp
        while (length(snp_queue) > 0) {
            current_snp <- snp_queue[1]
            snp_queue <- snp_queue[-1]
            neighbours <- neighbours_of_snp[[current_snp]]
            if (is.null(neighbours)) {
                next
            }
            for (neighbour_idx in seq_len(nrow(neighbours))) {
                neighbour_snp <- neighbours$to[neighbour_idx]
                # An "opposite" edge flips the orientation across it, a "same"
                # edge carries it through unchanged -- so the orientation this
                # edge implies for the neighbour is the current SNP's, XORed
                # with the edge's flip.
                implied_orientation <- as.integer(
                    xor(orientation[[current_snp]] == 1L, neighbours$flip[neighbour_idx])
                )
                if (is.na(orientation[neighbour_snp])) {
                    orientation[neighbour_snp] <- implied_orientation
                    block[neighbour_snp] <- block_id
                    snp_queue <- c(snp_queue, neighbour_snp)
                } else if (orientation[[neighbour_snp]] != implied_orientation) {
                    # Each undirected edge is walked from both ends, so the same
                    # conflict is seen twice; `edge_idx` identifies the original
                    # row so it is only ever recorded once.
                    is_edge_inconsistent[neighbours$edge_idx[neighbour_idx]] <- TRUE
                }
            }
        }
    }

    conflicting_edges <- edges[is_edge_inconsistent, , drop = FALSE]
    if (nrow(conflicting_edges) > 0) {
        # Both endpoints of a conflicting edge are in one block by construction,
        # so either endpoint identifies the block the conflict belongs to.
        conflicts_per_block <- table(block[conflicting_edges$snp_a])
        logger::log_warn(
            "{nrow(conflicting_edges)} edge(s) contradict the phase they were assigned, in ",
            "{length(conflicts_per_block)} block(s); flagging those SNPs block_conflict = TRUE"
        )
    }

    conflicted_block_ids <- as.integer(names(conflicts_per_block))

    # Resolved before the tibble() call rather than inside it: tibble() builds
    # columns sequentially in its own scope, so a `block = ` column defined
    # there masks this `block` lookup vector for every argument after it.
    block_of_snp <- unname(block[snp_ids])
    n_conflicts_of_snp <- dplyr::coalesce(as.integer(conflicts_per_block[as.character(block_of_snp)]), 0L)

    tibble::tibble(
        snp_id = snp_ids,
        block = block_of_snp,
        allele_on_h1 = dplyr::if_else(unname(orientation[snp_ids]) == 0L, "REF", "ALT"),
        n_block_conflicts = n_conflicts_of_snp,
        block_conflict = block_of_snp %in% conflicted_block_ids
    )
}

#' Assign each SNP to the gene(s) whose molecules it can be attributed to
#'
#' Deliberately distinct from `add_snp_gene_names()`, which comma-joins every
#' overlapping gene into a single display label. A SNP overlapping two gene
#' bodies on the same strand cannot have its molecules attributed to either
#' one without guessing, and is dropped entirely. A SNP overlapping genes on
#' different strands is not actually ambiguous at the molecule level: a
#' read's own alignment strand (once corrected for the BAM's sense/antisense
#' orientation, see `.infer_bam_strand_orientation()`) picks out which gene's
#' transcript it came from. Such a SNP is kept with one row per
#' strand-resolvable candidate gene, flagged `ambiguous = TRUE`, for
#' `haplotype_expression_by_molecule()` to resolve per molecule using
#' `molecule_read_strand()`. A candidate gene sharing its strand with another
#' candidate at the same SNP is not resolvable even by strand and is dropped
#' from that SNP's candidates. Pre-evaluation on real data found this
#' recovers real SNPs (24 of 81 multi-gene overlaps in one donor's chrX<20Mb
#' het set were opposite-strand); exon/splice-based recovery of same-strand
#' overlaps was evaluated and rejected as not worth the added complexity for
#' 6x less yield.
#'
#' Any SNP can still take part in `phase_snps()` regardless of its gene
#' assignment here, since phasing does not depend on gene assignment.
#'
#' @param snp_info A data.frame/tibble, required, with columns `snp_id`,
#'   `chrom`, `pos`.
#' @param gene_anno A data.frame/tibble, required, with columns `chrom`,
#'   `start`, `end`, `gene_name`, `strand` (`"+"`/`"-"`), one row per gene body.
#'
#' @return A tibble with columns `snp_id`, `gene_name`, `gene_strand`, and
#'   `ambiguous` (`TRUE` if the SNP overlaps more than one gene overall, so
#'   this candidate must still be matched against a molecule's own strand
#'   before use; `FALSE` if `gene_name` is this SNP's only candidate and
#'   applies regardless of strand).
#'
#' @family molecule-level allele counting functions
#' @export
assign_snp_genes <- function(snp_info, gene_anno) {
    required_snp_cols <- c("snp_id", "chrom", "pos")
    missing_snp_cols <- setdiff(required_snp_cols, colnames(snp_info))
    if (length(missing_snp_cols) > 0) {
        stop("snp_info is missing required column(s): ", paste(missing_snp_cols, collapse = ", "))
    }
    required_gene_cols <- c("chrom", "start", "end", "gene_name", "strand")
    missing_gene_cols <- setdiff(required_gene_cols, colnames(gene_anno))
    if (length(missing_gene_cols) > 0) {
        stop("gene_anno is missing required column(s): ", paste(missing_gene_cols, collapse = ", "))
    }

    snp_gr <- plyranges::as_granges(snp_info, seqnames = chrom, start = pos, end = pos)
    gene_gr <- plyranges::as_granges(gene_anno, seqnames = chrom, start = start, end = end)

    hits <- IRanges::findOverlaps(snp_gr, gene_gr)
    snp_gene_hits <- tibble::tibble(
        snp_id = snp_info$snp_id[S4Vectors::queryHits(hits)],
        gene_name = gene_anno$gene_name[S4Vectors::subjectHits(hits)],
        gene_strand = gene_anno$strand[S4Vectors::subjectHits(hits)]
    ) %>%
        dplyr::distinct()

    genes_per_snp <- snp_gene_hits %>% dplyr::summarise(n_genes = dplyr::n_distinct(gene_name), .by = snp_id)

    # Two candidate genes sharing a strand at the same SNP are still
    # indistinguishable even once a molecule's strand is known, so they are
    # dropped from that strand's candidates; a strand contributing exactly
    # one gene survives as a molecule-resolvable candidate.
    snp_gene_hits %>%
        dplyr::mutate(n_genes_same_strand = dplyr::n(), .by = c(snp_id, gene_strand)) %>%
        dplyr::filter(n_genes_same_strand == 1) %>%
        dplyr::inner_join(genes_per_snp, by = "snp_id") %>%
        dplyr::transmute(snp_id, gene_name, gene_strand, ambiguous = n_genes > 1)
}


#' Add read-backed molecule phase to a SNPData object's donor SNP metadata
#'
#' Extracts molecule-level allele calls from each library's BAM files, phases them
#' with `phase_snps()`, orients the resulting blocks to X1/X2 against
#' `assign_xci()`'s already-stored EM phase (see `.orient_phase_blocks()`),
#' and writes the result into `donor_snp_info`. Additive only:
#' `assign_xci()`'s own phase calls are never replaced, and
#' `haplotype_expression()` is unaffected unless it is told to read the new
#' columns.
#'
#' @inheritSection assign_xci Phase is inferred from expression, not genotyped
#'
#' @param x A SNPData object, required, that has already been fit by
#'   `assign_xci()` or `assign_xci_by_clonotype()` (see
#'   `.has_xci_diagnostics()`), with a `donor` column in `barcode_info`. Each
#'   donor's cells must all carry the same `library_id`, since a donor's BAM
#'   files are looked up by its library; a donor spanning two libraries is an
#'   error. `library_id` may instead be unset for every cell, in which case the
#'   object is treated as a single library, but a mix of labelled and
#'   unlabelled cells is an error.
#' @param bam_files A named character vector or list, optional (default
#'   `NULL`, taking the paths recorded in `library_info(x)$bam_files` by
#'   `import_cellsnp()` or `add_library_bams()`, and erroring if none were),
#'   `library_id = path(s)`, giving the indexed BAM file or files holding that
#'   library's reads. Names must match `barcode_info(x)$library_id`, or, for an
#'   object with no library labels at all, be a single entry of any name.
#'   Several files under one library are pooled per molecule, so a molecule
#'   split across them votes once with all of its reads; listing the same file
#'   twice is an error, as it would double every count it contributes.
#'   Libraries whose donors are all `"doublet"` or `"unassigned"` contribute
#'   nothing, since neither is a real donor with its own genotype to phase
#'   against.
#' @param target_chrom Character scalar (default `"chrX"`, matching
#'   `assign_xci()`'s own restriction). Canonical chromosome to restrict
#'   het-SNP selection to.
#' @param min_mapq,min_baseq,threads Integer (defaults 20, 10, 4). Passed to
#'   `extract_snp_calls()`.
#' @param min_molecules,min_cells,error_rate,min_llr Integer/numeric (defaults
#'   5, 2, 0.05, 3). Passed to `phase_snps()`, which accepts a SNP pair as an
#'   edge on the likelihood ratio between the two haplotype hypotheses, backed
#'   by molecules from at least `min_cells` distinct cells; see its
#'   \sQuote{Accepting an edge} section.
#'
#' @return A SNPData object with `donor_snp_info` gaining the columns
#'   `allele_on_x1_em` (the EM's own phase, unchanged, given its own name for
#'   symmetry with `allele_on_x1_molecule`), `allele_on_x1_molecule`,
#'   `phase_block`, `phase_source` (`"read_backed_propagated"`,
#'   `"read_backed_anchor"` or `"em"`, recording how much molecule evidence
#'   stands behind each SNP's phase; only the first is derived from molecules
#'   spanning the SNP itself, and `startsWith(phase_source, "read_backed")`
#'   selects the two block-based cases together), `phase_conflict`, and a
#'   re-derived
#'   `allele_on_x1`: the EM's value where present (so nothing downstream that
#'   already reads `allele_on_x1` changes behaviour), else the molecule value,
#'   else `NA` where `phase_conflict` is `TRUE`. Also carries a
#'   `"molecule_calls"` attribute (a tibble with columns `donor`, `barcode`,
#'   `umi`, `snp_id`, `allele`, `transcript_strand` (`"+"`/`"-"`/`NA`, the
#'   molecule's inferred transcript strand, see
#'   `.infer_bam_strand_orientation()` and `molecule_read_strand()`, used by
#'   `haplotype_expression_by_molecule()` to resolve SNPs `assign_snp_genes()`
#'   flagged `ambiguous`), the per-donor `molecule_snp_alleles()` output
#'   already computed here), which `haplotype_expression_by_molecule()` reads
#'   straight off the object it is given: nothing needs passing by hand or
#'   re-extracting from the BAM. Being an attribute, it does not survive
#'   operations that rebuild the object, so subset before this call rather than
#'   after. A second attribute, `"bam_calibration"`,
#'   records one row per BAM file scanned with columns `bam_file`,
#'   `orientation` (`"sense"`/`"antisense"`/`NA` where it could not be
#'   inferred), `n_ts_reads`, `concordance`, and `n_scanned`, so the strand
#'   call applied to each file's molecules can be inspected after the fact.
#'
#' @family molecule-level allele counting functions
#' @family X-chromosome inactivation functions
#' @export
add_molecule_phase <- function(
    x,
    bam_files = NULL,
    target_chrom = "chrX",
    min_mapq = 20L,
    min_baseq = 10L,
    min_molecules = 5L,
    min_cells = 2L,
    error_rate = 0.05,
    min_llr = 3,
    threads = 4L
) {
    if (!.has_xci_diagnostics(x)) {
        stop("No stored XCI diagnostics found. Run assign_xci(x) first.")
    }
    if (is.null(bam_files)) {
        bam_files <- .stored_bam_files(x)
    }
    bam_files <- .as_library_bam_list(bam_files)

    snp_info <- snp_info(x)
    if (!"chrom_canonical" %in% colnames(snp_info)) {
        stop("No canonical chromosome names available; SNPData must be built with a 'chrom' column.")
    }
    x_chrom <- filter_snps(x, chrom_canonical == target_chrom)
    donor_snp_info <- donor_snp_info(x)

    # Deferred and cached rather than computed up front: donor_het_status_df()
    # requires a zygosity source, and a caller whose only donors are
    # "doublet"/"unassigned" (filtered out inside the helper before any donor
    # is ever queried) should see that short-circuit rather than an unrelated
    # zygosity error.
    het_status <- NULL
    snp_ids_fn <- function(donor_id) {
        if (is.null(het_status)) {
            het_status <<- donor_het_status_df(x_chrom) %>% dplyr::filter(zygosity == "het")
        }
        unique(het_status$snp_id[het_status$donor == donor_id])
    }
    extracted <- .extract_and_phase_donor_molecules(
        x_chrom,
        bam_files = bam_files,
        snp_ids_fn = snp_ids_fn,
        min_mapq = min_mapq,
        min_baseq = min_baseq,
        min_molecules = min_molecules,
        min_cells = min_cells,
        error_rate = error_rate,
        min_llr = min_llr,
        threads = threads
    )
    per_donor <- extracted$per_donor
    calibration <- extracted$calibration
    if (length(per_donor) == 0) {
        logger::log_warn("No phase blocks formed for any donor; SNPData unchanged.")
        return(x)
    }

    per_donor_phase <- purrr::imap(per_donor, function(donor_result, donor_id) {
        anchors <- donor_snp_info %>%
            dplyr::filter(donor == donor_id, !is.na(allele_on_x1)) %>%
            dplyr::transmute(snp_id, allele_on_x1_em = allele_on_x1)

        oriented <- .orient_phase_blocks(donor_result$phase, anchors, donor = donor_id)
        oriented$donor <- donor_id
        list(oriented = oriented, per_snp = donor_result$per_snp)
    })

    new_cols <- dplyr::bind_rows(purrr::map(per_donor_phase, "oriented"))
    molecule_calls <- dplyr::bind_rows(purrr::map(per_donor_phase, "per_snp"))
    if (nrow(new_cols) == 0) {
        logger::log_warn("No phase blocks could be oriented for any donor; SNPData unchanged.")
        return(x)
    }

    new_cols$zygosity_source <- zygosity_source(x)

    # Step 1 adds only brand-new columns, so unmatched existing rows (this
    # donor_snp_info row was never touched by phasing) pass through untouched
    # -- critically, this cannot clobber any pre-existing allele_on_x1.
    x <- add_donor_snp_metadata(
        x,
        new_cols[, c(
            "snp_id",
            "donor",
            "zygosity_source",
            "allele_on_x1_molecule",
            "phase_block",
            "phase_source",
            "phase_conflict"
        )],
        join_by = c("snp_id", "donor", "zygosity_source"),
        overwrite = TRUE
    )

    # Step 2 recomputes the resolved allele_on_x1 over the *complete* table
    # (every row, not just the ones touched above) so the follow-up merge is a
    # row-complete 1:1 match and cannot null out an untouched row the way a
    # partial-coverage overwrite could.
    resolved <- donor_snp_info(x) %>%
        dplyr::mutate(
            allele_on_x1_em = allele_on_x1,
            allele_on_x1 = dplyr::case_when(
                !is.na(allele_on_x1) ~ allele_on_x1,
                phase_conflict %in% TRUE ~ NA_character_,
                TRUE ~ allele_on_x1_molecule
            )
        ) %>%
        dplyr::select(snp_id, donor, zygosity_source, allele_on_x1_em, allele_on_x1)

    x <- add_donor_snp_metadata(x, resolved, join_by = c("snp_id", "donor", "zygosity_source"), overwrite = TRUE)

    # BAM extraction is the expensive step this function already pays for;
    # attaching the per-molecule calls lets haplotype_expression_by_molecule()
    # reuse them instead of re-extracting from the BAM a second time.
    attr(x, "molecule_calls") <- molecule_calls
    # Strand calibration decides how every ambiguous SNP in a file is resolved,
    # so it is recorded rather than only logged: a file whose orientation came
    # out `NA`, or whose concordance was low enough that
    # `.infer_bam_strand_orientation()` fell back to a bare majority, is the
    # first thing to check when phase blocks look wrong.
    attr(x, "bam_calibration") <- calibration
    x
}


# ==============================================================================
# Orienting read-backed phase blocks to X1/X2
#
# phase_snps() assigns each block an arbitrary local label (H1 = whichever
# allele is REF at that block's first-visited SNP). Orienting H1 to the X1/X2
# labels used elsewhere in the package requires an external reference: a phase
# block's own molecules cannot supply this, because a true escapee's expression
# doesn't track XCI state by definition -- correlating a block against active_x
# would fail on exactly the genes this feature exists to rescue, the same way
# assign_xci()'s own keep_llr filter does. Anchors -- SNPs assign_xci() already
# phased via the EM -- are used instead: whether H1 matches X1 or X2 is read off
# any anchor reachable in the same connected component, and propagated to the
# rest of the block.
#
# Note what this does and does not buy. The *relative* phase within a block is
# genuinely physical: two SNPs seen on one molecule are on one chromosome, and
# that is observed, not inferred. The *absolute* orientation to X1/X2 is not --
# it is inherited wholesale from the EM anchors, which are expression-derived
# (see assign_xci()'s "Phase is inferred from expression, not genotyped"). A
# gene whose EM phase is inverted therefore has its whole read-backed block
# oriented to match that inversion, silently and without conflict, since every
# anchor in the block agrees. Molecules are single transcripts, so blocks never
# span genes and no cross-gene linkage exists to expose it. Read-backed phasing
# refines phase within a gene; it does not replace DNA-based phasing.
# ==============================================================================

#' Orient read-backed phase blocks to X1/X2 using assign_xci()'s EM phase
#'
#' A `phase_snps()` block's H1/H2 labelling is arbitrary and local to that
#' block; this maps it onto the same X1/X2 convention `assign_xci()` and
#' `haplotype_expression()` use, by finding SNPs in the block already phased
#' by the EM ("anchors") and reading off whether H1 agrees with X1 or X2 at
#' each. The chromosome the anchors are estimating is a single physical object,
#' so independently derived anchors within one component should never
#' legitimately disagree; where they do, this is treated as a signal to
#' investigate (most likely a low-power or noisy per-gene EM fit, occasionally a
#' spurious `phase_snps()` edge), not as evidence to average away.
#'
#' Agreement among anchors is therefore evidence of a consistent fit, not of a
#' correct one: the anchors are expression-derived, so a systematically inverted
#' gene yields unanimous anchors pointing the wrong way. See the note above
#' `.orient_phase_blocks()` in the source.
#'
#' @param phase A tibble as returned by `phase_snps()`, with columns
#'   `snp_id`, `block`, `allele_on_h1`, and optionally `block_conflict`, for a
#'   single donor. A SNP in a block `phase_snps()` flagged `block_conflict`
#'   (its edges cannot all be satisfied at once) is reported
#'   `phase_conflict = TRUE` here regardless of how well its anchors agree:
#'   unanimous anchors settle which way round to put a block, not whether the
#'   block's own internal phase is sound.
#' @param anchors A tibble with columns `snp_id` and `allele_on_x1_em` (the
#'   EM-derived phase for that donor, restricted to SNPs where it is not
#'   `NA`).
#' @param donor Character scalar, used only for log messages.
#'
#' @details
#' Per connected component (`block`), every reachable anchor's vote on
#' whether H1 corresponds to X1 or X2 is tallied:
#' \itemize{
#'   \item All anchors agree: oriented confidently.
#'   \item No anchors reachable: left unoriented (`allele_on_x1_molecule`
#'     `NA`); a block's own molecules are not used to orient it, since that
#'     reproduces the EM's own escapee blind spot.
#'   \item >= 3 anchors, exactly one disagrees with an otherwise-unanimous
#'     majority: oriented from the majority; the outlier anchor's own row
#'     is still flagged `phase_conflict = TRUE` for inspection (a repeat
#'     offender across donors would point to a reference-bias locus).
#'   \item Any other split (including exactly 2 anchors disagreeing, or >= 2
#'     anchors dissenting from a majority): not resolved. Every SNP in the
#'     component gets `allele_on_x1_molecule = NA` and `phase_conflict =
#'     TRUE`, since there is no principled way to tell which anchor is
#'     wrong.
#' }
#'
#' An anchor `phase_snps()` never linked to any other SNP (most commonly a
#' gene with only one heterozygous SNP) still carries its own EM-derived
#' phase; it becomes its own block of one (negative `phase_block`, to
#' stay visually distinct from `phase_snps()`'s positive block ids) rather
#' than being left out entirely.
#'
#' @return A tibble with columns `snp_id`, `phase_block`, `allele_on_x1_molecule`
#'   ("REF"/"ALT"/`NA`), `phase_source`, and `phase_conflict` (logical).
#'   Includes both SNPs from `phase` and any unlinked anchor from `anchors`.
#'   `phase_source` records how each SNP's value was arrived at, which differs
#'   in how much molecule evidence stands behind it:
#'   \describe{
#'     \item{`"read_backed_propagated"`}{A non-anchor SNP in a block. Its phase
#'       relative to the block's anchors was observed on molecules spanning
#'       both, so this is the only case with molecule evidence for the SNP
#'       itself.}
#'     \item{`"read_backed_anchor"`}{A SNP in a block that is itself an EM
#'       anchor. The value is the EM's, corroborated by agreeing with its
#'       block, not derived from molecules.}
#'     \item{`"em"`}{An anchor `phase_snps()` never linked to any partner,
#'       returned as its own singleton block with a negative `phase_block`. The
#'       value is the EM's, copied verbatim; no molecule was involved.}
#'   }
#'   The two read-backed values share a prefix so `startsWith(phase_source,
#'   "read_backed")` selects them together.
#'
#' @keywords internal
.orient_phase_blocks <- function(phase, anchors, donor = NA_character_) {
    phase_anchors <- phase %>%
        dplyr::inner_join(anchors, by = "snp_id") %>%
        dplyr::mutate(implies_h1_is_x1 = allele_on_h1 == allele_on_x1_em)

    block_summary <- phase_anchors %>%
        dplyr::summarise(
            n_anchors = dplyr::n(),
            n_x1 = sum(implies_h1_is_x1),
            n_x2 = sum(!implies_h1_is_x1),
            .by = block
        ) %>%
        dplyr::mutate(
            orientation = dplyr::case_when(
                n_x1 == n_anchors ~ "h1_is_x1",
                n_x2 == n_anchors ~ "h1_is_x2",
                n_anchors >= 3 & n_x1 == n_anchors - 1L ~ "h1_is_x1",
                n_anchors >= 3 & n_x2 == n_anchors - 1L ~ "h1_is_x2",
                TRUE ~ "conflict"
            )
        )

    conflicted_blocks <- block_summary$block[block_summary$orientation == "conflict"]
    if (length(conflicted_blocks) > 0) {
        logger::log_warn(
            "[{donor}] {length(conflicted_blocks)} phase block(s) have anchors that ",
            "disagree with no clear majority; leaving unoriented and flagging phase_conflict"
        )
    }

    # An anchor whose own vote disagrees with its block's resolved orientation
    # is a minority-of-one outlier (the only way a block still resolves despite
    # a dissent) -- flagged even though the block orientation itself is trusted.
    anchor_status <- phase_anchors %>%
        dplyr::left_join(dplyr::select(block_summary, block, orientation), by = "block") %>%
        dplyr::mutate(
            own_vote = dplyr::if_else(implies_h1_is_x1, "h1_is_x1", "h1_is_x2"),
            is_outlier_anchor = orientation %in% c("h1_is_x1", "h1_is_x2") & own_vote != orientation
        )

    if (any(anchor_status$is_outlier_anchor)) {
        outliers <- dplyr::filter(anchor_status, is_outlier_anchor)
        logger::log_warn(
            "[{donor}] {nrow(outliers)} anchor SNP(s) disagree with an otherwise-consistent ",
            "phase block, flagging phase_conflict: {paste(outliers$snp_id, collapse = ', ')}"
        )
    }

    # phase_snps() always supplies block_conflict, but the column is optional
    # here so a caller assembling a phase table by hand is not forced to invent
    # one; absent means nothing found a contradiction, which is not the same as
    # an unknown and so is FALSE rather than NA.
    block_conflict_flag <- rep(FALSE, nrow(phase))
    if ("block_conflict" %in% colnames(phase)) {
        block_conflict_flag <- phase[["block_conflict"]]
    }

    linked_snps <- phase %>%
        dplyr::left_join(dplyr::select(block_summary, block, orientation), by = "block") %>%
        dplyr::mutate(
            allele_on_x1_molecule = dplyr::case_when(
                is.na(orientation) ~ NA_character_,
                orientation == "conflict" ~ NA_character_,
                orientation == "h1_is_x1" ~ allele_on_h1,
                orientation == "h1_is_x2" ~ dplyr::if_else(allele_on_h1 == "REF", "ALT", "REF")
            ),
            # An anchor's own value comes from the EM and is merely corroborated
            # by the block it sits in; a non-anchor's is reached by propagating
            # along edges observed on molecules. Both are read-backed in the
            # sense that a block links them, but only the second has any
            # molecule evidence for the SNP itself, so they are named apart.
            phase_source = dplyr::if_else(
                snp_id %in% anchors$snp_id,
                "read_backed_anchor",
                "read_backed_propagated"
            ),
            phase_conflict = !is.na(orientation) & orientation == "conflict"
        ) %>%
        dplyr::left_join(dplyr::select(anchor_status, snp_id, is_outlier_anchor), by = "snp_id") %>%
        # Two independent ways a block can be untrustworthy, folded into one
        # flag because they mean the same thing downstream: do not rely on this
        # SNP's phase. The anchors can disagree about how to orient the block
        # (above), or the block's own edges can contradict each other, in which
        # case phase_snps() already found no orientation satisfies them all --
        # a block that cannot be internally consistent is not made trustworthy
        # by anchors that happen to agree on which way round to put it.
        dplyr::mutate(
            phase_conflict = phase_conflict |
                dplyr::coalesce(is_outlier_anchor, FALSE) |
                dplyr::coalesce(.env$block_conflict_flag, FALSE)
        ) %>%
        dplyr::select(snp_id, phase_block = block, allele_on_x1_molecule, phase_source, phase_conflict)

    # An anchor that phase_snps() never linked to any partner (e.g. a
    # single-heterozygous-SNP gene) still has a perfectly good EM-derived
    # phase of its own -- there is just nothing to pool it with. Rather than
    # leaving it invisible to haplotype_expression_by_molecule(), it becomes
    # its own block of one, mirroring how the original prototype treated an
    # unphased SNP ("a block of one, with H1 defined as REF"). Negative ids
    # keep these visually and numerically distinct from phase_snps()'s
    # (always positive) block ids.
    unlinked_anchors <- anchors %>%
        dplyr::filter(!snp_id %in% phase$snp_id) %>%
        dplyr::transmute(
            snp_id,
            phase_block = -dplyr::row_number(),
            allele_on_x1_molecule = allele_on_x1_em,
            # Copied verbatim from the EM: no molecule ever linked this SNP to
            # another, so nothing here is read-backed and labelling it so would
            # misreport the one case with no molecule evidence at all.
            phase_source = "em",
            phase_conflict = FALSE
        )

    dplyr::bind_rows(linked_snps, unlinked_anchors)
}

# Resolve donor/library/BAM inputs, extract per-molecule allele calls for
# each donor's het SNPs, and read-backed phase them with phase_snps() --
# everything add_molecule_phase() and molecule_haplotype_counts() both need
# before they part ways: the former orients blocks to X1/X2 against
# assign_xci()'s EM anchors, the latter reports them block-local. `snp_ids_fn`
# is called per donor as `snp_ids_fn(donor)` and must return the snp_ids
# to extract for that donor (the het-SNP set, however the caller derived it).
#
# Returns a named list, one entry per donor with phase blocks, each holding
# `phase` (phase_snps() output) and `per_snp` (molecule_snp_alleles() output
# with transcript_strand and donor columns, the per-molecule calls
# downstream functions attach as their "molecule_calls" attribute). A donor
# contributing no rows (no het SNPs, or no phase blocks formed) is absent.
.extract_and_phase_donor_molecules <- function(
    x,
    bam_files,
    snp_ids_fn,
    min_mapq,
    min_baseq,
    min_molecules,
    min_cells,
    error_rate,
    min_llr,
    threads
) {
    barcode_info <- barcode_info(x)
    if (!"donor" %in% colnames(barcode_info)) {
        stop("SNPData object has no donor assignments; phasing is done separately for each donor.")
    }
    if (is.null(bam_files)) {
        bam_files <- .stored_bam_files(x)
    }
    bam_files <- .as_library_bam_list(bam_files)

    # BAM files are keyed by library rather than by donor because that is what
    # they are a property of: one library's BAM holds all of its donors' cells,
    # and the object already records which library each cell came from. A
    # donor's files are therefore looked up, not asked for, which also makes it
    # impossible for the caller to point a donor at another library's reads --
    # where the same barcode names a different cell entirely.
    donor_library <- .donor_library_map(barcode_info)
    if (all(is.na(donor_library$library_id))) {
        # No labels anywhere: the object is a single implicit library, so one
        # entry covers every donor whatever the caller happened to name it.
        if (length(bam_files) != 1) {
            stop(
                "barcode_info$library_id is unset for every cell, so the object is a single library and ",
                "bam_files must have exactly one entry; it has ",
                length(bam_files),
                ". Label each cell's library to phase more than one."
            )
        }
        donor_library$library_id <- names(bam_files)
    }

    unknown_libraries <- setdiff(names(bam_files), donor_library$library_id)
    if (length(unknown_libraries) > 0) {
        stop(
            "bam_files names not found in barcode_info$library_id: ",
            paste(unknown_libraries, collapse = ", ")
        )
    }

    # "doublet"/"unassigned" are not real donors -- a doublet's genotype is a
    # mix of two cells' and an unassigned cell has no confident genotype, so
    # neither has a meaningful het-SNP set to phase against. Dropped here the
    # same way assign_xci() excludes them from its own per-donor EM fit.
    non_donor_labels <- intersect(donor_library$donor, c("doublet", "unassigned"))
    if (length(non_donor_labels) > 0) {
        logger::log_warn(
            "Excluding non-donor label(s) from phasing: {paste(non_donor_labels, collapse = ', ')}"
        )
    }
    donor_library <- donor_library %>%
        dplyr::filter(!donor %in% c("doublet", "unassigned"), library_id %in% names(bam_files))
    if (nrow(donor_library) == 0) {
        logger::log_warn("No real donors have BAM files supplied for their library; nothing to phase.")
        return(list(per_donor = list(), calibration = NULL))
    }

    bam_files <- .check_bam_paths(bam_files[unique(donor_library$library_id)])
    donor_bams <- stats::setNames(bam_files[donor_library$library_id], donor_library$donor)
    calibration <- .calibrate_bam_strands(unique(unlist(bam_files, use.names = FALSE)))

    per_donor <- purrr::map(names(donor_bams), function(donor_id) {
        donor_snp_ids <- snp_ids_fn(donor_id)
        donor_het_snp_info <- snp_info(x) %>%
            dplyr::filter(snp_id %in% donor_snp_ids) %>%
            dplyr::select(snp_id, chrom, pos, ref, alt)
        if (nrow(donor_het_snp_info) == 0) {
            logger::log_warn("[{donor_id}] no het SNPs to phase; skipping")
            return(NULL)
        }

        # A donor's whole barcode set is safe to use against every one of its
        # library's files: `library_id` is what disambiguates a barcode shared
        # with another library, and these files are that library's.
        donor_barcodes <- barcode_info$barcode[barcode_info$donor == donor_id]
        per_file <- purrr::map(donor_bams[[donor_id]], function(bam_file) {
            extracted <- extract_snp_calls(
                bam_file,
                donor_het_snp_info,
                barcodes = donor_barcodes,
                min_mapq = min_mapq,
                min_baseq = min_baseq,
                threads = threads
            )
            # Orientation is a property of the file's pipeline, so it is applied
            # before pooling: two of a donor's files may well be calibrated
            # differently, and pooled reads carry no record of where they came
            # from.
            orientation <- calibration$orientation[calibration$bam_file == bam_file]
            molecule_strand <- molecule_read_strand(extracted$reads)
            aligned_strand <- molecule_strand$strand
            if (is.na(orientation)) {
                molecule_strand$transcript_strand <- rep(NA_character_, length(aligned_strand))
            } else if (orientation == "sense") {
                molecule_strand$transcript_strand <- aligned_strand
            } else {
                molecule_strand$transcript_strand <- ifelse(aligned_strand == "+", "-", "+")
            }
            list(
                tallies = extracted$tallies,
                molecule_strand = dplyr::select(molecule_strand, barcode, umi, transcript_strand)
            )
        })
        pooled <- .pool_donor_calls(per_file)

        per_snp <- molecule_snp_alleles(pooled$tallies) %>%
            dplyr::left_join(pooled$molecule_strand, by = c("barcode", "umi"))
        phase <- phase_snps(
            per_snp,
            min_molecules = min_molecules,
            min_cells = min_cells,
            error_rate = error_rate,
            min_llr = min_llr
        )
        if (nrow(phase) == 0) {
            logger::log_warn("[{donor_id}] no phase blocks formed from {nrow(donor_het_snp_info)} het SNPs")
            return(NULL)
        }
        per_snp$donor <- donor_id
        list(phase = phase, per_snp = per_snp)
    })
    names(per_donor) <- names(donor_bams)
    list(per_donor = purrr::compact(per_donor), calibration = calibration)
}

# The paths recorded against each library at import, as the named list
# `bam_files` would have been given as. Libraries with no stored path are left
# out entirely, so a half-populated object fails the same way an incomplete
# `bam_files` argument would rather than silently phasing only some donors.
.stored_bam_files <- function(x) {
    stored <- library_info(x)
    stored <- stored[lengths(stored$bam_files) > 0, , drop = FALSE]
    if (nrow(stored) == 0) {
        stop(
            "No bam_files given and none recorded on the object. Supply them here, or record them with ",
            "import_cellsnp(..., bam_files = ) or add_library_bams()."
        )
    }
    stats::setNames(stored$bam_files, stored$library_id)
}

# A BAM file listed twice under one library would be extracted twice and its
# tallies summed by `.pool_donor_calls()`, silently doubling every read behind
# that library's molecules. A missing index is just as quiet but costlier:
# `.prefilter_bam()`'s `samtools view -M -L` seeks by index, and without one the
# region-restricted scan degrades to streaming the whole file once per donor.
.check_bam_paths <- function(bam_files) {
    for (library_id in names(bam_files)) {
        paths <- bam_files[[library_id]]
        missing_paths <- paths[!file.exists(paths)]
        if (length(missing_paths) > 0) {
            stop("[", library_id, "] BAM file(s) not found: ", paste(missing_paths, collapse = ", "))
        }
        resolved_paths <- normalizePath(paths)
        if (anyDuplicated(resolved_paths) > 0) {
            stop(
                "[",
                library_id,
                "] the same BAM file is listed more than once: ",
                paste(unique(resolved_paths[duplicated(resolved_paths)]), collapse = ", "),
                ". Repeated files would double every read count they contribute."
            )
        }
        is_indexed <- vapply(paths, .has_bam_index, logical(1), USE.NAMES = FALSE)
        if (!all(is_indexed)) {
            stop(
                "[",
                library_id,
                "] BAM file(s) have no index: ",
                paste(paths[!is_indexed], collapse = ", "),
                ". Index them with samtools index; extraction seeks by index and is far slower without one."
            )
        }
        bam_files[[library_id]] <- resolved_paths
    }
    bam_files
}

.has_bam_index <- function(bam_file) {
    # Both index naming conventions: alongside the full name (reads.bam.bai)
    # and replacing the extension (reads.bai), in either .bai or .csi form.
    index_paths <- c(
        paste0(bam_file, c(".bai", ".csi")),
        sub("\\.bam$", ".bai", bam_file),
        sub("\\.bam$", ".csi", bam_file)
    )
    any(file.exists(index_paths))
}

# Which library each donor's cells came from. This is derived from the object
# rather than asked of the caller: `library_id` is a property of the cell, so
# the object already knows, and a donor whose cells span two libraries breaks
# the assumption every BAM lookup here rests on -- that a donor's reads live in
# exactly one library's files -- and so is an error rather than a guess.
# An object with no library labels at all is treated as one implicit library.
.donor_library_map <- function(barcode_info) {
    n_missing <- sum(is.na(barcode_info$library_id))
    if (n_missing > 0 && n_missing < nrow(barcode_info)) {
        stop(
            "barcode_info$library_id is set for some cells but not others (",
            n_missing,
            " of ",
            nrow(barcode_info),
            " unlabelled). Label every cell's library, or none."
        )
    }
    donor_library <- dplyr::distinct(barcode_info, donor, library_id)
    donor_library <- donor_library[!is.na(donor_library$donor), , drop = FALSE]
    split_donors <- donor_library$donor[duplicated(donor_library$donor)]
    if (length(split_donors) > 0) {
        stop(
            "Donor(s) with cells in more than one library: ",
            paste(unique(split_donors), collapse = ", "),
            ". add_molecule_phase() looks up a donor's BAM files by its library, so each donor must sit in one."
        )
    }
    donor_library
}

# Strand calibration is a property of the BAM's pipeline, not of the donors
# read out of it, so it is computed once per file and up front: eagerly, so a
# file that cannot be calibrated is reported before any expensive extraction is
# paid for, and once, so two donors sharing a library BAM cannot calibrate it
# inconsistently. Unlike extraction this scan starts at the head of the file
# and ignores the index, so repeating it per donor is the one cost that does
# not shrink with the chromosome filter.
.calibrate_bam_strands <- function(bam_paths) {
    per_file_rows <- purrr::map(bam_paths, function(bam_file) {
        calibration <- tryCatch(
            .infer_bam_strand_orientation(bam_file),
            error = function(e) {
                logger::log_warn(
                    "could not infer strand orientation for {bam_file} ({conditionMessage(e)}); ",
                    "strand-ambiguous SNPs will be unresolved for molecules from this file"
                )
                NULL
            }
        )
        if (is.null(calibration)) {
            return(tibble::tibble(
                bam_file = bam_file,
                orientation = NA_character_,
                n_ts_reads = NA_integer_,
                concordance = NA_real_,
                n_scanned = NA_integer_
            ))
        }
        tibble::tibble(
            bam_file = bam_file,
            orientation = calibration$orientation,
            n_ts_reads = as.integer(calibration$n_ts_reads),
            concordance = as.numeric(calibration$concordance),
            n_scanned = as.integer(calibration$n_scanned)
        )
    })
    dplyr::bind_rows(per_file_rows)
}

# Pool one donor's per-file extractions into the shape a single file would have
# produced. The two tables need different reductions because they mean
# different things. Read tallies for a molecule split across files are partial
# counts of one vote, so they must be summed *before* `molecule_snp_alleles()`
# picks a winner: binding alone would instead take the argmax of the per-file
# counts and discard the losing file's reads, and voting per file then binding
# would emit one row per file per molecule, so `phase_snps()` would count a
# single molecule several times towards `min_molecules`. Transcript strand is
# already one call per molecule per file, so it is majority-voted instead, and
# must come out at one row per molecule or the left join onto the calls fans
# out. Molecules seen only in an uncalibrated file are simply absent here and
# pick up `NA` from that join; a genuine "+"/"-" disagreement between files
# resolves to `NA` rather than a guess.
.pool_donor_calls <- function(per_file) {
    tallies <- dplyr::bind_rows(purrr::map(per_file, "tallies")) %>%
        dplyr::summarise(n_calls = sum(n_calls), .by = c(barcode, umi, snp_id, allele))

    molecule_strand <- dplyr::bind_rows(purrr::map(per_file, "molecule_strand")) %>%
        dplyr::filter(!is.na(transcript_strand)) %>%
        dplyr::count(barcode, umi, transcript_strand, name = "n_files") %>%
        dplyr::slice_max(n_files, n = 1, by = c(barcode, umi), with_ties = TRUE) %>%
        dplyr::summarise(
            transcript_strand = dplyr::if_else(dplyr::n() == 1L, transcript_strand[1], NA_character_),
            .by = c(barcode, umi)
        )

    list(tallies = tallies, molecule_strand = molecule_strand)
}

#' Infer whether a BAM's reads are sense or antisense to their transcript
#'
#' Demultiplexing tools such as Flexiplex reorient reads before alignment,
#' and 5' vs 3' protocol data end up flipped in opposite directions (5'
#' reads sense to the transcript, 3' reads its reverse complement), with no
#' record of which happened left in the BAM. A spliced read's `ts:A:+/-` tag
#' (the transcript strand minimap2 calls from the GT-AG splice-junction
#' signal) is independent of that flip, so comparing a read's own alignment
#' strand to its `ts` value reveals the orientation for that read; pooling
#' this over enough `ts`-tagged reads calibrates the whole BAM, since one BAM
#' is assumed to use a single protocol throughout. Reads are scanned in
#' batches and the scan stops as soon as the split is decisive, rather than
#' reading the whole file.
#'
#' @param bam_file Path to a BAM (read sequentially; need not be indexed).
#' @param batch_size Reads scanned per batch. Default 5000.
#' @param min_ts_reads `ts`-tagged reads required before checking for
#'   confidence. Default 200.
#' @param min_concordance Fraction of `ts`-tagged reads that must agree on
#'   sense/antisense to stop early. Default 0.95.
#' @param max_reads Reads scanned before giving up rather than looping over
#'   the whole file. Default 200000.
#'
#' @return A list with `orientation` (`"sense"` or `"antisense"`,
#'   whichever a majority of `ts`-tagged reads support), `n_ts_reads`,
#'   `concordance` (fraction of those agreeing with `orientation`), and
#'   `n_scanned` (total reads read to reach the decision).
#'
#' @keywords internal
.infer_bam_strand_orientation <- function(
    bam_file,
    batch_size = 5000L,
    min_ts_reads = 200L,
    min_concordance = 0.95,
    max_reads = 200000L
) {
    bam_conn <- Rsamtools::BamFile(bam_file, yieldSize = batch_size)
    open(bam_conn)
    on.exit(close(bam_conn), add = TRUE)

    param <- Rsamtools::ScanBamParam(
        tag = "ts",
        flag = Rsamtools::scanBamFlag(
            isSecondaryAlignment = FALSE,
            isSupplementaryAlignment = FALSE,
            isUnmappedQuery = FALSE,
            isDuplicate = FALSE
        )
    )

    n_sense <- 0L
    n_antisense <- 0L
    n_scanned <- 0L

    repeat {
        galn <- GenomicAlignments::readGAlignments(bam_conn, param = param)
        if (length(galn) == 0) {
            break
        }
        n_scanned <- n_scanned + length(galn)

        # A read whose alignment strand matches the transcript strand minimap2
        # called from its splice junctions was sequenced sense to the
        # transcript; a mismatch means the pipeline reverse-complemented it.
        transcript_strand_tag <- S4Vectors::mcols(galn)$ts
        has_ts_tag <- !is.na(transcript_strand_tag)
        if (any(has_ts_tag)) {
            aligned_strand <- as.character(BiocGenerics::strand(galn))[has_ts_tag]
            is_sense_read <- aligned_strand == transcript_strand_tag[has_ts_tag]
            n_sense <- n_sense + sum(is_sense_read)
            n_antisense <- n_antisense + sum(!is_sense_read)
        }

        n_ts_reads <- n_sense + n_antisense
        if (n_ts_reads >= min_ts_reads) {
            concordance <- max(n_sense, n_antisense) / n_ts_reads
            if (concordance >= min_concordance) {
                return(list(
                    orientation = if (n_sense >= n_antisense) "sense" else "antisense",
                    n_ts_reads = n_ts_reads,
                    concordance = concordance,
                    n_scanned = n_scanned
                ))
            }
        }
        if (n_scanned >= max_reads) {
            break
        }
    }

    n_ts_reads <- n_sense + n_antisense
    if (n_ts_reads == 0) {
        stop("No ts-tagged reads found in ", bam_file, " after scanning ", n_scanned, " reads")
    }
    concordance <- max(n_sense, n_antisense) / n_ts_reads
    logger::log_warn(
        "Strand orientation for {bam_file} inconclusive after {n_scanned} reads ",
        "({n_ts_reads} ts-tagged, {round(concordance * 100, 1)}% concordant); using majority"
    )
    list(
        orientation = if (n_sense >= n_antisense) "sense" else "antisense",
        n_ts_reads = n_ts_reads,
        concordance = concordance,
        n_scanned = n_scanned
    )
}

#' Per-gene molecule counts on read-backed haplotype blocks, without XCI
#'
#' Extracts molecule-level allele calls from each donor's BAM files, phases
#' heterozygous SNPs with `phase_snps()`, and counts each gene's molecules
#' once per phase block they fall in -- the general-purpose counterpart to
#' `haplotype_expression_by_molecule()` for genes with no X-inactivation
#' signal to orient blocks against.
#'
#' @details
#' Unlike XCI, an autosomal gene's two haplotypes have no external signal
#' (silencing skew, an EM fit) to say which physical chromosome copy is
#' "haplotype 1" versus "haplotype 2" -- see
#' \code{\link{add_molecule_phase}} and `.orient_phase_blocks()` in the
#' source for how XCI supplies that signal via `assign_xci()`'s anchors.
#' `phase_snps()`'s block-local `H1`/`H2` labels are therefore reported
#' as-is: consistent for every cell of one donor within a single phase
#' block (both haplotypes come from that donor's own genome, so every cell
#' shares them), but arbitrary and \strong{not comparable across blocks or
#' across donors}. A gene whose heterozygous SNPs are not all spanned by a
#' shared molecule splits into more than one block, and each block's `H1`
#' is unrelated to any other block's `H1` -- there is no valid way to sum
#' them into one gene-level count, so this function keeps them as separate
#' rows rather than guessing an orientation.
#'
#' A symmetric imbalance statistic, e.g. \code{pmin(h1_count, h2_count) /
#' (h1_count + h2_count)}, is still meaningful per row and can be compared
#' or pooled across blocks and donors, since it does not depend on which
#' label is `H1`.
#'
#' @param x A SNPData object, required, with a `donor` column in
#'   `barcode_info` and a zygosity source established (Vireo genotypes read
#'   at import, or \code{\link{infer_zygosity}}).
#' @param bam_files A named character vector or list, optional (default
#'   `NULL`, taking the paths recorded in `library_info(x)$bam_files`),
#'   `library_id = path(s)`. See \code{\link{add_molecule_phase}}'s
#'   `bam_files` argument for the full matching rules; the same rules apply
#'   here.
#' @param target_chrom Character vector, optional (default `NULL`, every
#'   chromosome). Canonical chromosome(s) to restrict het-SNP selection to.
#' @param min_mapq,min_baseq,threads Integer (defaults 20, 10, 4). Passed to
#'   `extract_snp_calls()`.
#' @param min_molecules,min_cells,error_rate,min_llr Integer/numeric (defaults
#'   5, 2, 0.05, 3). Passed to `phase_snps()`, which accepts a SNP pair as an
#'   edge on the likelihood ratio between the two haplotype hypotheses, backed
#'   by molecules from at least `min_cells` distinct cells; see its
#'   \sQuote{Accepting an edge} section.
#'
#' @return A tibble with one row per (`donor`, `gene_name`, `phase_block`),
#'   columns `h1_count`, `h2_count` (integer; molecules voting for each
#'   block-local haplotype label, where a tied molecule is dropped as
#'   ambiguous, as in `haplotype_expression_by_molecule()`), `coverage`
#'   (integer; their sum), and `n_molecules` (integer; distinct molecules
#'   backing the block, before the haplotype vote, matching
#'   `haplotype_expression_by_molecule()`'s `dominant_molecules`). A (donor,
#'   gene) pair contributes one row per phase block its heterozygous SNPs
#'   fall in; see Details for why these rows are not pooled into one gene
#'   total.
#'
#' @family molecule-level allele counting functions
#' @export
#'
#' @examples
#' \dontrun{
#' hap_blocks <- molecule_haplotype_counts(snp_data, bam_files = c(lib1 = "lib1.bam"))
#'
#' # A symmetric imbalance statistic, valid to compare across blocks/donors
#' # even though h1_count/h2_count individually are not
#' hap_blocks$maf <- pmin(hap_blocks$h1_count, hap_blocks$h2_count) / hap_blocks$coverage
#' }
molecule_haplotype_counts <- function(
    x,
    bam_files = NULL,
    target_chrom = NULL,
    min_mapq = 20L,
    min_baseq = 10L,
    min_molecules = 5L,
    min_cells = 2L,
    error_rate = 0.05,
    min_llr = 3,
    threads = 4L
) {
    snp_gene_map <- snp_gene_map(x)
    if (nrow(snp_gene_map) == 0) {
        stop(
            "This object carries no SNP-to-gene map. It is built by import_cellsnp() from a gene ",
            "annotation with a strand column; set it on an existing object with ",
            "snp_gene_map(x) <- assign_snp_genes(snp_info(x), gene_anno)."
        )
    }
    if (is.null(bam_files)) {
        bam_files <- .stored_bam_files(x)
    }
    bam_files <- .as_library_bam_list(bam_files)

    x_scope <- x
    if (!is.null(target_chrom)) {
        snp_info <- snp_info(x)
        if (!"chrom_canonical" %in% colnames(snp_info)) {
            stop("No canonical chromosome names available; SNPData must be built with a 'chrom' column.")
        }
        x_scope <- filter_snps(x, chrom_canonical %in% target_chrom)
    }

    het_status <- NULL
    snp_ids_fn <- function(donor_id) {
        if (is.null(het_status)) {
            het_status <<- donor_het_status_df(x_scope) %>% dplyr::filter(zygosity == "het")
        }
        unique(het_status$snp_id[het_status$donor == donor_id])
    }
    extracted <- .extract_and_phase_donor_molecules(
        x_scope,
        bam_files = bam_files,
        snp_ids_fn = snp_ids_fn,
        min_mapq = min_mapq,
        min_baseq = min_baseq,
        min_molecules = min_molecules,
        min_cells = min_cells,
        error_rate = error_rate,
        min_llr = min_llr,
        threads = threads
    )
    per_donor <- extracted$per_donor
    if (length(per_donor) == 0) {
        stop("No phase blocks formed for any donor.")
    }

    # phase_snps() blocks are numbered independently per donor call, so an
    # unlinked SNP would otherwise collide with another donor's block 1; each
    # donor's blocks are namespaced by donor here for that reason (mirroring
    # .orient_phase_blocks()'s negative-id trick for unlinked anchors, but
    # there is no anchor step here to make a block of one for an unlinked SNP,
    # so it is simply absent, same as phase_snps() leaves it).
    phase <- purrr::imap(per_donor, function(donor_result, donor_id) {
        donor_phase <- donor_result$phase
        donor_phase$donor <- donor_id
        donor_phase
    }) %>%
        dplyr::bind_rows() %>%
        dplyr::rename(phase_block = block)
    molecule_calls <- dplyr::bind_rows(purrr::map(per_donor, "per_snp"))

    calls <- .molecule_gene_phase_calls(molecule_calls, snp_gene_map, phase, orientation_col = "allele_on_h1")

    # A molecule votes across every SNP of the block it covers; majority wins,
    # a tie is ambiguous (residual base-calling noise once phase is accounted
    # for) and dropped rather than guessed -- the same rule
    # haplotype_expression_by_molecule() applies per gene, applied here per
    # (gene, block) since blocks are not pooled.
    molecule_votes <- calls %>%
        dplyr::summarise(
            n_h1 = sum(is_oriented_allele),
            n_h2 = sum(!is_oriented_allele),
            .by = c(donor, gene_name, phase_block, barcode, umi)
        ) %>%
        dplyr::mutate(haplotype = dplyr::case_when(n_h1 > n_h2 ~ "H1", n_h2 > n_h1 ~ "H2", TRUE ~ "ambiguous"))

    n_molecules <- molecule_votes %>%
        dplyr::summarise(n_molecules = dplyr::n(), .by = c(donor, gene_name, phase_block))

    molecule_votes %>%
        dplyr::filter(haplotype != "ambiguous") %>%
        dplyr::summarise(
            h1_count = sum(haplotype == "H1"),
            h2_count = sum(haplotype == "H2"),
            .by = c(donor, gene_name, phase_block)
        ) %>%
        dplyr::mutate(coverage = h1_count + h2_count) %>%
        dplyr::left_join(n_molecules, by = c("donor", "gene_name", "phase_block")) %>%
        dplyr::arrange(donor, gene_name, phase_block)
}

# Join per-molecule allele calls to a phase table and the SNP-to-gene map,
# resolving strand-ambiguous SNPs against each molecule's own transcript
# strand. Shared by haplotype_expression_by_molecule() (phase already
# oriented to X1/X2 by .orient_phase_blocks()) and molecule_haplotype_counts()
# (phase left in phase_snps()'s own block-local H1/H2), since the join and the
# strand-resolution rule are identical either way -- only what "oriented
# allele" means differs, named genetically as is_oriented_allele here and
# relabelled by each caller (is_x1, is_h1) for its own convention.
#
# `phase` must have columns snp_id, donor, phase_block, and a column named
# `orientation_col` giving the allele ("REF"/"ALT") this call treats as
# oriented.
.molecule_gene_phase_calls <- function(molecule_calls, snp_gene_map, phase, orientation_col) {
    phase <- dplyr::rename(phase, .oriented_allele = dplyr::all_of(orientation_col))
    calls <- molecule_calls %>%
        dplyr::inner_join(phase, by = c("snp_id", "donor")) %>%
        dplyr::inner_join(snp_gene_map, by = "snp_id", relationship = "many-to-many") %>%
        dplyr::filter(!ambiguous | (!is.na(transcript_strand) & transcript_strand == gene_strand)) %>%
        dplyr::mutate(is_oriented_allele = allele == .oriented_allele) %>%
        dplyr::select(-.oriented_allele)

    if (nrow(calls) == 0) {
        stop("No molecule calls could be matched to a phased, singly-mapped-gene SNP.")
    }
    calls
}

# Per-(donor, gene, phase_block) molecule bookkeeping shared by
# haplotype_expression_by_molecule() and molecule_haplotype_counts(): which
# block has the most distinct molecules backing it, how many molecules sit in
# any other block, and which blocks are included in `counted` once
# `pool_blocks` is applied. Pooling blocks into one gene-level count is only
# valid when every SNP's `is_oriented_allele` means the same thing across
# blocks (an externally, globally oriented phase, as XCI's EM anchors give);
# `calls` here is agnostic to that and only counts molecules, leaving the
# caller to decide whether pooling the returned blocks is sound for its phase
# source.
.molecule_gene_block_counts <- function(calls, pool_blocks) {
    blocks <- calls %>%
        dplyr::distinct(donor, gene_name, phase_block, barcode, umi) %>%
        dplyr::count(donor, gene_name, phase_block, name = "molecules")
    best_block <- blocks %>%
        dplyr::slice_max(molecules, n = 1, by = c(donor, gene_name), with_ties = FALSE) %>%
        dplyr::select(donor, gene_name, phase_block, dominant_molecules = molecules)
    stranded <- blocks %>%
        dplyr::anti_join(best_block, by = c("donor", "gene_name", "phase_block")) %>%
        dplyr::summarise(n_stranded_molecules = sum(molecules), .by = c(donor, gene_name))

    counted <- calls
    if (!pool_blocks) {
        counted <- dplyr::semi_join(counted, best_block, by = c("donor", "gene_name", "phase_block"))
    }
    blocks_counted <- counted %>%
        dplyr::distinct(donor, gene_name, phase_block) %>%
        dplyr::summarise(n_blocks_pooled = dplyr::n(), .by = c(donor, gene_name))

    list(best_block = best_block, stranded = stranded, blocks_counted = blocks_counted, counted = counted)
}

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

    per_chunk <- purrr::map(chunks, function(ii) {
        sub <- galn[ii]
        calls <- .map_bases_to_reads(sub, snp_gr, snp_info$snp_id) %>%
            dplyr::filter(base_quality >= min_baseq) %>%
            dplyr::mutate(
                barcode = S4Vectors::mcols(sub)$CB[read_idx],
                umi = S4Vectors::mcols(sub)$UB[read_idx],
                qname = S4Vectors::mcols(sub)$qname[read_idx],
                strand = as.character(BiocGenerics::strand(sub))[read_idx],
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

    results <- purrr::map(unique(snp_chrom), function(chrom) {
        snp_on_chrom <- which(snp_chrom == chrom)
        read_on_chrom <- galn_chrom[read_of_range] == chrom

        hits <- IRanges::findOverlaps(
            IRanges::IRanges(snp_pos[snp_on_chrom], width = 1L),
            unlisted_ref[read_on_chrom]
        )
        range_idx <- which(read_on_chrom)[S4Vectors::subjectHits(hits)]
        snp_idx <- snp_on_chrom[S4Vectors::queryHits(hits)]
        read_idx <- read_of_range[range_idx]
        query_pos <- snp_pos[snp_idx] - query_shift[range_idx]

        qual_chars <- as.character(
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
            base_quality = as.integer(charToRaw(paste0(qual_chars, collapse = ""))) - 33L
        )
    })

    dplyr::bind_rows(results)
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
    args <- c("view", "-@", threads, "-M", "-L", bed_file, "-q", min_mapq, "-F", "0xD04")
    if (is.null(barcodes)) {
        args <- c(args, "-e", shQuote("[CB] && [UB]"))
    } else {
        barcode_file <- tempfile(fileext = ".txt")
        writeLines(barcodes, barcode_file)
        on.exit(unlink(barcode_file), add = TRUE)
        # samtools permits only one -d/-D tag filter, so a second tag requirement
        # has to go through the expression filter -e instead.
        args <- c(args, "-D", paste0("CB:", barcode_file), "-e", shQuote("[UB]"))
    }

    status <- system2("samtools", c(args, "-u", "-o", out_bam, bam_file))
    if (status != 0) {
        unlink(out_bam)
        warning("samtools pre-filter failed; falling back to reading the full BAM")
        return(NULL)
    }
    out_bam
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

        ts <- S4Vectors::mcols(galn)$ts
        has_ts <- !is.na(ts)
        if (any(has_ts)) {
            read_strand <- as.character(BiocGenerics::strand(galn))[has_ts]
            concordant <- read_strand == ts[has_ts]
            n_sense <- n_sense + sum(concordant)
            n_antisense <- n_antisense + sum(!concordant)
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

#' Resolve the allele call for each (molecule, SNP)
#'
#' Duplicate reads of one molecule vote on the allele at each SNP. Only REF
#' and ALT are retained; OTH is a sequencing error at a known biallelic site
#' and carries no haplotype information.
#'
#' @param tallies A tibble, required, as returned by
#'   `extract_snp_calls()$tallies`, with columns `barcode`, `umi`, `snp_id`,
#'   `allele`, and `n_calls`.
#'
#' @return A tibble with one row per (`barcode`, `umi`, `snp_id`) and columns
#'   `barcode`, `umi`, `snp_id`, `allele` (the majority call, "REF" or "ALT"),
#'   and `n_calls`.
#'
#' @family molecule-level allele counting functions
#' @export
molecule_snp_alleles <- function(tallies) {
    tallies %>%
        dplyr::slice_max(n_calls, n = 1, by = c(barcode, umi, snp_id), with_ties = FALSE) %>%
        dplyr::filter(allele %in% c("REF", "ALT"))
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
#' @param reads A tibble, required, as returned by `extract_snp_calls()$reads`,
#'   with columns `barcode`, `umi`, `qname`, `strand`.
#'
#' @return A tibble with one row per (`barcode`, `umi`) and columns
#'   `barcode`, `umi`, `strand` (the majority call, `"+"` or `"-"`).
#'
#' @family molecule-level allele counting functions
#' @export
molecule_read_strand <- function(reads) {
    reads %>%
        dplyr::count(barcode, umi, strand, name = "n_reads") %>%
        dplyr::slice_max(n_reads, n = 1, by = c(barcode, umi), with_ties = FALSE) %>%
        dplyr::select(barcode, umi, strand)
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
#' @param per_snp A tibble, required, as returned by `molecule_snp_alleles()`,
#'   with columns `barcode`, `umi`, `snp_id`, `allele`.
#' @param min_molecules Integer (default 5). Molecules required to accept a
#'   SNP pair as an edge.
#' @param min_consistency Numeric, in `[0, 1]` (default 0.9). Fraction of
#'   those molecules that must agree on the same/opposite relation for the
#'   edge to be accepted.
#'
#' @return A tibble with columns `snp_id`, `block` (integer phase-block id,
#'   unique within this call), and `allele_on_h1` ("REF" or "ALT", the
#'   allele carried by haplotype 1 at that SNP; H1 is an arbitrary label local
#'   to each block, not oriented to X1/X2). SNPs that could not be linked to
#'   any other SNP by an accepted edge are absent from the result.
#'
#' @family molecule-level allele counting functions
#' @export
phase_snps <- function(per_snp, min_molecules = 5L, min_consistency = 0.9) {
    molecule_snps <- per_snp %>%
        dplyr::arrange(barcode, umi, snp_id) %>%
        dplyr::summarise(snps = list(snp_id), alleles = list(allele), .by = c(barcode, umi)) %>%
        dplyr::filter(lengths(snps) >= 2)

    if (nrow(molecule_snps) == 0) {
        return(tibble::tibble(snp_id = character(), block = integer(), allele_on_h1 = character()))
    }

    pair_rows <- purrr::map(seq_len(nrow(molecule_snps)), function(i) {
        s <- molecule_snps$snps[[i]]
        a <- molecule_snps$alleles[[i]]
        idx <- utils::combn(length(s), 2)
        tibble::tibble(snp_a = s[idx[1, ]], snp_b = s[idx[2, ]], same = a[idx[1, ]] == a[idx[2, ]])
    })

    edges <- dplyr::bind_rows(pair_rows) %>%
        dplyr::summarise(n = dplyr::n(), n_same = sum(same), .by = c(snp_a, snp_b)) %>%
        dplyr::mutate(
            relation = dplyr::if_else(n_same >= n - n_same, "same", "opposite"),
            consistency = pmax(n_same, n - n_same) / n
        ) %>%
        dplyr::filter(n >= min_molecules, consistency >= min_consistency)

    if (nrow(edges) == 0) {
        return(tibble::tibble(snp_id = character(), block = integer(), allele_on_h1 = character()))
    }

    snp_ids <- sort(unique(c(edges$snp_a, edges$snp_b)))
    orientation <- stats::setNames(rep(NA_integer_, length(snp_ids)), snp_ids)
    block <- stats::setNames(rep(NA_integer_, length(snp_ids)), snp_ids)

    # Undirected adjacency: every edge is walkable from both ends.
    adj <- split(
        rbind(
            data.frame(to = edges$snp_b, flip = edges$relation == "opposite"),
            data.frame(to = edges$snp_a, flip = edges$relation == "opposite")
        ),
        c(edges$snp_a, edges$snp_b)
    )

    block_id <- 0L
    for (seed in snp_ids) {
        if (!is.na(orientation[seed])) {
            next
        }
        block_id <- block_id + 1L
        orientation[seed] <- 0L
        block[seed] <- block_id
        queue <- seed
        while (length(queue) > 0) {
            node <- queue[1]
            queue <- queue[-1]
            nbrs <- adj[[node]]
            if (is.null(nbrs)) {
                next
            }
            for (k in seq_len(nrow(nbrs))) {
                to <- nbrs$to[k]
                want <- as.integer(xor(orientation[[node]] == 1L, nbrs$flip[k]))
                if (is.na(orientation[to])) {
                    orientation[to] <- want
                    block[to] <- block_id
                    queue <- c(queue, to)
                }
            }
        }
    }

    tibble::tibble(
        snp_id = snp_ids,
        block = unname(block[snp_ids]),
        allele_on_h1 = dplyr::if_else(unname(orientation[snp_ids]) == 0L, "REF", "ALT")
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
#'   `snp_id`, `block`, `allele_on_h1`, for a single donor.
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
#'   ("REF"/"ALT"/`NA`), `phase_source` (always `"read_backed"`), and
#'   `phase_conflict` (logical). Includes both SNPs from `phase` and any
#'   unlinked anchor from `anchors`.
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

    linked <- phase %>%
        dplyr::left_join(dplyr::select(block_summary, block, orientation), by = "block") %>%
        dplyr::mutate(
            allele_on_x1_molecule = dplyr::case_when(
                is.na(orientation) ~ NA_character_,
                orientation == "conflict" ~ NA_character_,
                orientation == "h1_is_x1" ~ allele_on_h1,
                orientation == "h1_is_x2" ~ dplyr::if_else(allele_on_h1 == "REF", "ALT", "REF")
            ),
            phase_source = "read_backed",
            phase_conflict = !is.na(orientation) & orientation == "conflict"
        ) %>%
        dplyr::left_join(dplyr::select(anchor_status, snp_id, is_outlier_anchor), by = "snp_id") %>%
        dplyr::mutate(phase_conflict = phase_conflict | dplyr::coalesce(is_outlier_anchor, FALSE)) %>%
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
            phase_source = "read_backed",
            phase_conflict = FALSE
        )

    dplyr::bind_rows(linked, unlinked_anchors)
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

# Normalise `bam_files` into a named list, library_id = character vector of
# paths. Shape only: membership against the object and the state of the files
# themselves are checked later, so that an unrecognised library is reported
# before a missing file is.
.normalise_bam_files <- function(bam_files) {
    if (length(bam_files) == 0) {
        stop("bam_files is empty; supply at least one library's BAM file(s).")
    }
    nms <- names(bam_files)
    if (is.null(nms) || anyNA(nms) || any(!nzchar(nms))) {
        stop("bam_files must be named, library_id = path(s), matching barcode_info(x)$library_id.")
    }
    if (anyDuplicated(nms) > 0) {
        stop(
            "bam_files has repeated library_id name(s): ",
            paste(unique(nms[duplicated(nms)]), collapse = ", "),
            ". Give each library one entry listing all of its BAM files."
        )
    }
    bam_files <- as.list(bam_files)
    for (lib in nms) {
        paths <- bam_files[[lib]]
        if (!is.character(paths) || length(paths) == 0 || anyNA(paths)) {
            stop("bam_files[['", lib, "']] must be a non-empty character vector of BAM paths.")
        }
    }
    bam_files
}

# A BAM file listed twice under one library would be extracted twice and its
# tallies summed by `.pool_donor_calls()`, silently doubling every read behind
# that library's molecules. A missing index is just as quiet but costlier:
# `.prefilter_bam()`'s `samtools view -M -L` seeks by index, and without one the
# region-restricted scan degrades to streaming the whole file once per donor.
.check_bam_paths <- function(bam_files) {
    for (lib in names(bam_files)) {
        paths <- bam_files[[lib]]
        missing <- paths[!file.exists(paths)]
        if (length(missing) > 0) {
            stop("[", lib, "] BAM file(s) not found: ", paste(missing, collapse = ", "))
        }
        resolved <- normalizePath(paths)
        if (anyDuplicated(resolved) > 0) {
            stop(
                "[",
                lib,
                "] the same BAM file is listed more than once: ",
                paste(unique(resolved[duplicated(resolved)]), collapse = ", "),
                ". Repeated files would double every read count they contribute."
            )
        }
        indexed <- vapply(paths, .has_bam_index, logical(1), USE.NAMES = FALSE)
        if (!all(indexed)) {
            stop(
                "[",
                lib,
                "] BAM file(s) have no index: ",
                paste(paths[!indexed], collapse = ", "),
                ". Index them with samtools index; extraction seeks by index and is far slower without one."
            )
        }
        bam_files[[lib]] <- resolved
    }
    bam_files
}

.has_bam_index <- function(bam_file) {
    candidates <- c(
        paste0(bam_file, c(".bai", ".csi")),
        sub("\\.bam$", ".bai", bam_file),
        sub("\\.bam$", ".csi", bam_file)
    )
    any(file.exists(candidates))
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
    map <- dplyr::distinct(barcode_info, donor, library_id)
    map <- map[!is.na(map$donor), , drop = FALSE]
    split_donors <- map$donor[duplicated(map$donor)]
    if (length(split_donors) > 0) {
        stop(
            "Donor(s) with cells in more than one library: ",
            paste(unique(split_donors), collapse = ", "),
            ". add_molecule_phase() looks up a donor's BAM files by its library, so each donor must sit in one."
        )
    }
    map
}

# Strand calibration is a property of the BAM's pipeline, not of the donors
# read out of it, so it is computed once per file and up front: eagerly, so a
# file that cannot be calibrated is reported before any expensive extraction is
# paid for, and once, so two donors sharing a library BAM cannot calibrate it
# inconsistently. Unlike extraction this scan starts at the head of the file
# and ignores the index, so repeating it per donor is the one cost that does
# not shrink with the chromosome filter.
.calibrate_bam_strands <- function(bam_paths) {
    rows <- purrr::map(bam_paths, function(f) {
        calibration <- tryCatch(
            .infer_bam_strand_orientation(f),
            error = function(e) {
                logger::log_warn(
                    "could not infer strand orientation for {f} ({conditionMessage(e)}); ",
                    "strand-ambiguous SNPs will be unresolved for molecules from this file"
                )
                NULL
            }
        )
        if (is.null(calibration)) {
            return(tibble::tibble(
                bam_file = f,
                orientation = NA_character_,
                n_ts_reads = NA_integer_,
                concordance = NA_real_,
                n_scanned = NA_integer_
            ))
        }
        tibble::tibble(
            bam_file = f,
            orientation = calibration$orientation,
            n_ts_reads = as.integer(calibration$n_ts_reads),
            concordance = as.numeric(calibration$concordance),
            n_scanned = as.integer(calibration$n_scanned)
        )
    })
    dplyr::bind_rows(rows)
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
#' @param min_molecules,min_consistency Integer/numeric (defaults 5, 0.9).
#'   Passed to `phase_snps()`.
#'
#' @return A SNPData object with `donor_snp_info` gaining the columns
#'   `allele_on_x1_em` (the EM's own phase, unchanged, given its own name for
#'   symmetry with `allele_on_x1_molecule`), `allele_on_x1_molecule`,
#'   `phase_block`, `phase_source`, `phase_conflict`, and a re-derived
#'   `allele_on_x1`: the EM's value where present (so nothing downstream that
#'   already reads `allele_on_x1` changes behaviour), else the molecule value,
#'   else `NA` where `phase_conflict` is `TRUE`. Also carries a
#'   `"molecule_calls"` attribute (a tibble with columns `donor`, `barcode`,
#'   `umi`, `snp_id`, `allele`, `transcript_strand` (`"+"`/`"-"`/`NA`, the
#'   molecule's inferred transcript strand, see
#'   `.infer_bam_strand_orientation()` and `molecule_read_strand()`, used by
#'   `haplotype_expression_by_molecule()` to resolve SNPs `assign_snp_genes()`
#'   flagged `ambiguous`), the per-donor `molecule_snp_alleles()` output
#'   already computed here) so `haplotype_expression_by_molecule()` does not
#'   need to re-extract from the BAM. A second attribute, `"bam_calibration"`,
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
    min_consistency = 0.9,
    threads = 4L
) {
    if (!.has_xci_diagnostics(x)) {
        stop("No stored XCI diagnostics found. Run assign_xci(x) first.")
    }
    barcode_info <- barcode_info(x)
    if (!"donor" %in% colnames(barcode_info)) {
        stop("SNPData object has no donor assignments; add_molecule_phase() phases each donor separately.")
    }
    bam_files <- bam_files %||% .stored_bam_files(x)
    bam_files <- .normalise_bam_files(bam_files)

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
        logger::log_warn("No real donors have BAM files supplied for their library; SNPData unchanged.")
        return(x)
    }

    bam_files <- .check_bam_paths(bam_files[unique(donor_library$library_id)])
    donor_bams <- stats::setNames(bam_files[donor_library$library_id], donor_library$donor)
    calibration <- .calibrate_bam_strands(unique(unlist(bam_files, use.names = FALSE)))

    snp_info <- snp_info(x)
    if (!"chrom_canonical" %in% colnames(snp_info)) {
        stop("No canonical chromosome names available; SNPData must be built with a 'chrom' column.")
    }
    x_chrom <- filter_snps(x, chrom_canonical == target_chrom)
    het_status <- donor_het_status_df(x_chrom) %>%
        dplyr::filter(zygosity == "het", donor %in% names(donor_bams))

    donor_snp_info <- donor_snp_info(x)

    per_donor_phase <- purrr::map(names(donor_bams), function(d) {
        donor_snp_ids <- unique(het_status$snp_id[het_status$donor == d])
        this_snp_info <- snp_info(x_chrom) %>%
            dplyr::filter(snp_id %in% donor_snp_ids) %>%
            dplyr::select(snp_id, chrom, pos, ref, alt)
        if (nrow(this_snp_info) == 0) {
            logger::log_warn("[{d}] no het SNPs on {target_chrom}; skipping")
            return(NULL)
        }

        # A donor's whole barcode set is safe to use against every one of its
        # library's files: `library_id` is what disambiguates a barcode shared
        # with another library, and these files are that library's.
        donor_barcodes <- barcode_info$barcode[barcode_info$donor == d]
        per_file <- purrr::map(donor_bams[[d]], function(bam_file) {
            extracted <- extract_snp_calls(
                bam_file,
                this_snp_info,
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
        phase <- phase_snps(per_snp, min_molecules = min_molecules, min_consistency = min_consistency)
        if (nrow(phase) == 0) {
            logger::log_warn("[{d}] no phase blocks formed from {nrow(this_snp_info)} het SNPs")
            return(NULL)
        }

        anchors <- donor_snp_info %>%
            dplyr::filter(donor == d, !is.na(allele_on_x1)) %>%
            dplyr::transmute(snp_id, allele_on_x1_em = allele_on_x1)

        oriented <- .orient_phase_blocks(phase, anchors, donor = d)
        oriented$donor <- d
        per_snp$donor <- d
        list(oriented = oriented, per_snp = per_snp)
    })
    per_donor_phase <- purrr::compact(per_donor_phase)

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

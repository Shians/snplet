# ==============================================================================
# Molecule phasing: gene attribution and per-donor orchestration
#
# The entry points here sit on top of molecule_extraction.R (per-molecule
# allele calls from a BAM) and molecule_phase_graph.R (turning those calls
# into phased blocks): they resolve which BAM file(s) belong to each donor,
# extract and phase that donor's molecules, and write the result back onto a
# SNPData object or summarise it into gene-level counts.
#
# add_molecule_phase() orients blocks to X1/X2 against assign_xci()'s EM
# phase and stores read-backed phase alongside it. molecule_haplotype_counts()
# instead reports blocks in their own block-local H1/H2 labelling, for callers
# that need molecule-level counts without committing to an absolute
# orientation. assign_snp_genes() supplies the gene attribution both rely on
# for SNPs that overlap more than one gene body.
# ==============================================================================

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

    # Two candidates sharing a strand stay indistinguishable even once a
    # molecule's strand is known, so they are dropped; a strand with exactly
    # one gene survives.
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

    # Deferred and cached: a caller whose only donors are "doublet"/"unassigned"
    # should see that short-circuit rather than an unrelated zygosity error.
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

    # Step 1: new columns only, so this cannot clobber any pre-existing allele_on_x1.
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

    # Step 2: recomputed over the *complete* table so the merge is 1:1 and
    # cannot null out an untouched row.
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

    # Lets haplotype_expression_by_molecule() reuse this instead of re-extracting.
    attr(x, "molecule_calls") <- molecule_calls
    # Recorded, not just logged, so a bad orientation call is inspectable later.
    attr(x, "bam_calibration") <- calibration
    x
}


# Resolves donor/library/BAM inputs, extracts per-molecule allele calls for
# each donor's het SNPs, and read-backed phases them -- everything
# add_molecule_phase() and molecule_haplotype_counts() both need before they
# part ways (the former orients blocks to X1/X2, the latter reports them
# block-local). `snp_ids_fn(donor_id)` returns the het-SNP set to extract.
#
# Returns a named list, one entry per donor with phase blocks formed, each
# holding `phase` (phase_snps() output) and `per_snp` (molecule_snp_alleles()
# output plus transcript_strand and donor columns).
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

    # BAM files are keyed by library, not donor: that is what they are a
    # property of, and it stops a caller pointing a donor at another
    # library's reads.
    donor_library <- .donor_library_map(barcode_info)
    if (all(is.na(donor_library$library_id))) {
        # No labels anywhere: treat the object as one implicit library.
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

    # Neither has a genuine genotype to phase against; excluded the same way
    # assign_xci() excludes them from its own per-donor EM fit.
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
            # Applied before pooling: two of a donor's files may be calibrated
            # differently, and pooled reads carry no record of origin.
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

# Paths recorded against each library at import, in the shape a `bam_files`
# argument would take. Libraries with no stored path are left out, so a
# half-populated object fails rather than silently phasing only some donors.
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

# A duplicate path would double-count every read via .pool_donor_calls(); a
# missing index is quieter but costlier, since .prefilter_bam() seeks by index.
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

# Which library each donor's cells came from, derived from the object rather
# than asked of the caller. A donor spanning two libraries breaks the
# assumption every BAM lookup here rests on, so is an error, not a guess.
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

# Computed once per file, up front: strand calibration is a property of the
# BAM's pipeline, not of the donors read from it, so two donors sharing a file
# must not calibrate it inconsistently, and a failure surfaces before
# extraction is paid for.
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

# Pool one donor's per-file extractions into the shape a single file would
# have produced. Read tallies are partial counts of one vote and must be
# summed *before* molecule_snp_alleles() picks a winner, or a molecule split
# across files gets double-counted (or its losing file's reads discarded).
# Transcript strand is already one call per molecule per file, so it is
# majority-voted down to one row per molecule instead, to avoid fanning out
# the later join; a genuine disagreement between files resolves to NA.
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

        # A match between alignment strand and the ts tag means sense; a mismatch
        # means the pipeline reverse-complemented the read.
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

    # phase_snps() numbers blocks independently per call, so they are
    # namespaced by donor here to avoid collisions (e.g. two donors' block 1).
    phase <- purrr::imap(per_donor, function(donor_result, donor_id) {
        donor_phase <- donor_result$phase
        donor_phase$donor <- donor_id
        donor_phase
    }) %>%
        dplyr::bind_rows() %>%
        dplyr::rename(phase_block = block)
    molecule_calls <- dplyr::bind_rows(purrr::map(per_donor, "per_snp"))

    calls <- .molecule_gene_phase_calls(molecule_calls, snp_gene_map, phase, orientation_col = "allele_on_h1")

    # Majority wins per (gene, block); a tie is dropped, not guessed -- the
    # same rule haplotype_expression_by_molecule() applies per gene.
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

# Joins per-molecule allele calls to a phase table and the SNP-to-gene map,
# resolving strand-ambiguous SNPs against each molecule's transcript strand.
# Shared by haplotype_expression_by_molecule() (X1/X2-oriented phase) and
# molecule_haplotype_counts() (block-local H1/H2 phase): the join and
# strand-resolution rule are identical, only what "oriented allele" means
# differs, via the `orientation_col` argument.
#
# `phase` must have columns snp_id, donor, phase_block, and `orientation_col`
# giving the allele ("REF"/"ALT") this call treats as oriented.
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
# haplotype_expression_by_molecule() and molecule_haplotype_counts(): the
# dominant block, molecules stranded in other blocks, and which blocks
# `pool_blocks` includes in `counted`. Pooling is only valid when
# `is_oriented_allele` means the same thing across blocks (a globally oriented
# phase); this function only counts molecules and leaves that call to the caller.
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

# ==============================================================================
# Read-backed phasing graph
#
# molecule_snp_alleles() (see molecule_extraction.R) gives one allele call per
# (molecule, SNP). The functions here turn those calls into phase. phase_snps()
# builds a graph of SNP pairs linked by molecules that span both, accepts an
# edge where the evidence for "same haplotype" or "opposite" is strong enough,
# and reads off each connected component's relative phase from a spanning
# tree, arbitrarily labelled H1/H2 within that component. .orient_phase_blocks()
# then anchors H1/H2 to the X1/X2 labelling assign_xci() already uses
# elsewhere in the package, via SNPs the EM has independently phased.
# ==============================================================================

#' Phase heterozygous SNPs directly from the molecules that span them
#'
#' Two het SNPs observed on the same molecule sit on the same physical
#' haplotype, so their alleles co-occur non-randomly: mostly REF/REF and
#' ALT/ALT if the two REF alleles share a haplotype ("same"), or REF/ALT and
#' ALT/REF if they do not ("opposite"). This is read-backed phasing: long
#' reads supply it directly, with no statistical fit needed.
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
#' additive per molecule: each agreeing molecule adds that weight, each
#' disagreeing one subtracts it, and an edge is accepted once the total
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
#' a barcode: they flip that cell's molecules together rather than one at a
#' time. Molecules from a single cell therefore inflate \code{n} without
#' adding proportionate evidence, so `min_cells` requires an edge to be
#' corroborated across barcodes before it is accepted. It is a floor on
#' independence, not a second likelihood: the LLR itself still counts
#' molecules.
#'
#' The count is taken over the cells backing the relation the edge is accepted
#' on, not over every cell spanning the pair. The two differ precisely when a
#' cell dissents: four molecules from one cell reading "same" alongside one
#' from another reading "opposite" spans two barcodes but rests, as far as
#' "same" is concerned, on a single cell. Counting the dissenter towards the
#' relation it contradicts would let it vouch for exactly the single-cell edge
#' `min_cells` exists to reject.
#'
#' \code{e} is treated as independent across molecules, which sequencing
#' error is and contamination is not: ambient RNA is correlated within a cell
#' and a doublet flips many molecules together. Where contamination rather
#' than base-calling error dominates, the LLR is optimistic, so set \code{e}
#' from observed discordance rather than from a sequencing-error prior, and
#' treat `min_llr` as a threshold to calibrate rather than a p-value.
#'
#' @section Internally inconsistent blocks:
#' Orientations are fixed by a spanning tree of each component, so any edge
#' that closes a cycle is not needed to phase its endpoints, but it is still
#' an independent prediction of the relation between them. Where such an edge
#' contradicts the orientations already assigned, no assignment of alleles to
#' two haplotypes can satisfy every accepted edge at once. That is impossible
#' for a diploid genome, so at least one edge in the block is wrong: most
#' often a spurious link from ambient RNA, a mismapped paralogue, or an
#' undetected doublet whose two genotypes are read as one.
#'
#' Such blocks are still returned, oriented by the spanning tree as before:
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
    # (0, 0.5): at 0 the weight is infinite, at 0.5 it is 0, beyond it the test reads backwards.
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

    # One vote per pair of SNPs a molecule spans, "same" if the two alleles it
    # read agree; snp_pairs' two combn() rows index the molecule's SNPs and
    # alleles in lockstep.
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
            # Carried through for the `min_cells` count below.
            barcode = molecule_snps$barcode[molecule_idx]
        )
    })

    # LLR and per-relation cell counts implement @section Accepting an edge.
    weight_per_molecule <- log((1 - error_rate) / error_rate)
    edges <- dplyr::bind_rows(pair_votes) %>%
        dplyr::summarise(
            n = dplyr::n(),
            n_same = sum(same),
            # n_cells_same/opposite, not one count over the whole pair, since a
            # cell is independent evidence only for the relation it voted for.
            n_cells_same = dplyr::n_distinct(barcode[same]),
            n_cells_opposite = dplyr::n_distinct(barcode[!same]),
            .by = c(snp_a, snp_b)
        ) %>%
        dplyr::mutate(
            relation = dplyr::if_else(n_same >= n - n_same, "same", "opposite"),
            consistency = pmax(n_same, n - n_same) / n,
            n_cells = dplyr::if_else(relation == "same", n_cells_same, n_cells_opposite),
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
    # orientation[snp]: 0 = REF on H1, 1 = REF on H2; NA doubles as "unvisited".
    orientation <- stats::setNames(rep(NA_integer_, length(snp_ids)), snp_ids)
    block <- stats::setNames(rep(NA_integer_, length(snp_ids)), snp_ids)

    # Undirected adjacency: each edge appears twice (once per end), and
    # edge_idx traces a half-edge back to its row in `edges` for conflict
    # reporting below.
    neighbours_of_snp <- split(
        rbind(
            data.frame(to = edges$snp_b, flip = edges$relation == "opposite", edge_idx = seq_len(nrow(edges))),
            data.frame(to = edges$snp_a, flip = edges$relation == "opposite", edge_idx = seq_len(nrow(edges)))
        ),
        c(edges$snp_a, edges$snp_b)
    )

    # An edge reaching an already-oriented SNP closes a cycle; see @section
    # Internally inconsistent blocks for what a disagreement there means.
    is_edge_inconsistent <- rep(FALSE, nrow(edges))
    conflicts_per_block <- integer(0)

    # BFS: each unvisited SNP seeds a block, arbitrarily orienting its REF
    # allele to H1 (0); every SNP reached from it inherits an orientation
    # forced by the relations along the way.
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
                # "opposite" flips the orientation across the edge, "same"
                # carries it through: XOR the current SNP's orientation with flip.
                implied_orientation <- as.integer(
                    xor(orientation[[current_snp]] == 1L, neighbours$flip[neighbour_idx])
                )
                if (is.na(orientation[neighbour_snp])) {
                    orientation[neighbour_snp] <- implied_orientation
                    block[neighbour_snp] <- block_id
                    snp_queue <- c(snp_queue, neighbour_snp)
                } else if (orientation[[neighbour_snp]] != implied_orientation) {
                    # Walked from both ends, so record via edge_idx once, not twice.
                    is_edge_inconsistent[neighbours$edge_idx[neighbour_idx]] <- TRUE
                }
            }
        }
    }

    conflicting_edges <- edges[is_edge_inconsistent, , drop = FALSE]
    if (nrow(conflicting_edges) > 0) {
        # Either endpoint identifies the block: both are in it by construction.
        conflicts_per_block <- table(block[conflicting_edges$snp_a])
        logger::log_warn(
            "{nrow(conflicting_edges)} edge(s) contradict the phase they were assigned, in ",
            "{length(conflicts_per_block)} block(s); flagging those SNPs block_conflict = TRUE"
        )
    }

    conflicted_block_ids <- as.integer(names(conflicts_per_block))

    # Resolved outside tibble(): a `block = ` column defined inside it would
    # mask this `block` lookup vector for every argument after it.
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


# ==============================================================================
# Orienting read-backed phase blocks to X1/X2
#
# phase_snps() assigns each block an arbitrary local label (H1 = whichever
# allele is REF at that block's first-visited SNP). Orienting H1 to the X1/X2
# labels used elsewhere in the package requires an external reference: a phase
# block's own molecules cannot supply this, because a true escapee's expression
# doesn't track XCI state by definition. Correlating a block against active_x
# would fail on exactly the genes this feature exists to rescue, the same way
# assign_xci()'s own keep_llr filter does. Anchors (SNPs assign_xci() already
# phased via the EM) are used instead: whether H1 matches X1 or X2 is read off
# any anchor reachable in the same connected component, and propagated to the
# rest of the block.
#
# What this does and does not buy: the *relative* phase within a block is
# genuinely physical, since two SNPs seen on one molecule are on one
# chromosome, and that is observed, not inferred. The *absolute* orientation
# to X1/X2 is not: it is inherited wholesale from the EM anchors, which are
# expression-derived (see assign_xci()'s "Phase is inferred from expression,
# not genotyped"). A gene whose EM phase is inverted therefore has its whole
# read-backed block oriented to match that inversion, silently and without
# conflict, since every anchor in the block agrees. Molecules are single
# transcripts, so blocks never span genes and no cross-gene linkage exists to
# expose it. Read-backed phasing refines phase within a gene; it does not
# replace DNA-based phasing.
# ==============================================================================

#' Orient read-backed phase blocks to X1/X2 using assign_xci()'s EM phase
#'
#' A `phase_snps()` block's H1/H2 labelling is arbitrary and local to that
#' block. This maps it onto the same X1/X2 convention `assign_xci()` and
#' `haplotype_expression()` use, by finding SNPs in the block already phased
#' by the EM ("anchors") and reading off whether H1 agrees with X1 or X2 at
#' each. The chromosome the anchors are estimating is a single physical
#' object, so independently derived anchors within one component should never
#' legitimately disagree. Where they do, this is treated as a signal to
#' investigate (most likely a low-power or noisy per-gene EM fit, occasionally
#' a spurious `phase_snps()` edge), not as evidence to average away.
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

    # A minority-of-one anchor is flagged even though the block orientation
    # itself is trusted (the majority resolved it).
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

    # block_conflict is optional here (phase_snps() always supplies it, a
    # hand-built phase table need not); absent means FALSE, not unknown.
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
            # See @return: only the propagated case has molecule evidence for
            # the SNP itself; an anchor's value is the EM's, merely corroborated.
            phase_source = dplyr::if_else(
                snp_id %in% anchors$snp_id,
                "read_backed_anchor",
                "read_backed_propagated"
            ),
            phase_conflict = !is.na(orientation) & orientation == "conflict"
        ) %>%
        dplyr::left_join(dplyr::select(anchor_status, snp_id, is_outlier_anchor), by = "snp_id") %>%
        # Anchor disagreement and phase_snps()'s own edge conflicts are two
        # independent reasons not to trust a SNP's phase, folded into one flag.
        dplyr::mutate(
            phase_conflict = phase_conflict |
                dplyr::coalesce(is_outlier_anchor, FALSE) |
                dplyr::coalesce(.env$block_conflict_flag, FALSE)
        ) %>%
        dplyr::select(snp_id, phase_block = block, allele_on_x1_molecule, phase_source, phase_conflict)

    # See @details: an anchor phase_snps() never linked becomes its own
    # negative-id block of one rather than being dropped.
    unlinked_anchors <- anchors %>%
        dplyr::filter(!snp_id %in% phase$snp_id) %>%
        dplyr::transmute(
            snp_id,
            phase_block = -dplyr::row_number(),
            allele_on_x1_molecule = allele_on_x1_em,
            phase_source = "em",
            phase_conflict = FALSE
        )

    dplyr::bind_rows(linked_snps, unlinked_anchors)
}

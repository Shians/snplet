#' Phased active/inactive haplotype expression per gene
#'
#' Splits phased het SNPs' ALT/REF counts into active- and inactive-haplotype
#' reads using the phase \code{\link{assign_xci}} fit against its frozen
#' cell-state calls, computed separately within the X1-active and X2-active cell
#' groups, then reports one elected representative SNP per donor and gene (see
#' \sQuote{Gene-level election}). \code{by_snp = TRUE} reports every SNP instead,
#' which is the surface to inspect when a gene is missing or suspect.
#'
#' Unlike a per-group \code{pmin}/\code{pmax} on raw counts, the
#' allele-to-haplotype assignment comes from the stored per-donor phase
#' (\code{donor_snp_info$allele_on_x1}) and is applied consistently to both
#' groups, so the same physical allele can never be silently relabelled as
#' "active" in both — that contradiction is instead surfaced explicitly.
#'
#' Both grains span every donor-genotyped het SNP with a stored phase, including
#' SNPs \code{assign_xci} excluded from active-X calling as uninformative or
#' escaping, unless \code{xci_informative_only = TRUE}.
#'
#' @details
#' For each SNP and donor, the stored phase records which allele (REF or ALT)
#' the model assigned to the X1 cluster; the other is assigned to X2. That
#' assignment is inferred from expression, never genotyped -- see
#' \sQuote{Phase is inferred from expression, not genotyped}. The active
#' haplotype is the one named by \code{active_x}: X1 in the X1-active group, X2
#' in the X2-active group; the other X is taken to be silenced. Counts are
#' assigned to \code{active_count} / \code{inactive_count} accordingly.
#'
#' Under canonical XCI, the active haplotype dominates in both groups, so the
#' \emph{physical} allele that dominates should flip between the groups (X1's
#' allele in X1-active cells, X2's allele in X2-active cells). When the same
#' physical allele dominates in both groups — biologically impossible for a
#' cleanly inactivated SNP — \code{same_allele_dominant} is \code{TRUE},
#' identifying escaping or XCI-independent allelic imbalance rather than
#' flattening it. The \code{escape_fraction} column quantifies the magnitude:
#' the fraction of reads coming from the haplotype that should be silenced.
#'
#' \strong{\code{escape_fraction} is effectively bounded above by 0.5.} From
#' expression alone there is no way to tell a gene transcribed mostly from the
#' inactive X from a gene transcribed mostly from the active X whose phase is
#' recorded the other way round: the two produce identical read counts.
#' \code{\link{assign_xci}} resolves that symmetry by construction -- its phase
#' step picks whichever orientation makes the silenced allele the minority one,
#' and \code{pi_g} is clamped below 0.5 -- so the stored phase always names the
#' majority allele as active. The scale therefore runs from 0 (strict
#' inactivation) to 0.5 (complete escape, both haplotypes equal), and 0.5 is the
#' maximum escape that is measurable rather than merely the midpoint. Set
#' \code{escape_threshold}, and any null passed to \code{\link{test_escape}}, on
#' that scale.
#'
#' A consequence worth internalising: a gene expressing predominantly from the
#' inactive X is not reported as escaping at all. It is stored with inverted
#' phase and reads as an unusually \emph{clean} XCI gene -- see \sQuote{Genes
#' masked by phase inversion}.
#'
#' Because the pooled fit is held below 0.5, an individual group can still exceed
#' it by sampling -- most readily when the gene's true escape is near 0.5, where
#' each group lands on a random side. Such a group is marked
#' \code{phase_contradiction} and is \emph{reported, not dropped}: the same
#' signature arises from a quantification artefact (reference-mapping bias,
#' ambient RNA, low coverage) and from a complete escapee, and excluding it would
#' let sampling decide whether a real escapee appears at all.
#'
#' \code{phase_contradiction} and \code{same_allele_dominant} are related but not
#' equivalent. With both groups covered, \code{same_allele_dominant} is
#' \code{TRUE} exactly when the two groups disagree about which side of 0.5 they
#' sit on -- one flagged and the other not. Two groups that both exceed 0.5 would
#' still flip, since each is dominated by a different physical allele, but that
#' state does not arise from a fit produced by \code{\link{assign_xci}}: the
#' phase step would have inverted such a SNP.
#'
#' @inheritSection assign_xci Phase is inferred from expression, not genotyped
#'
#' @param x A SNPData object that had XCI diagnostics stored by
#'   \code{\link{assign_xci}} or
#'   \code{\link{assign_xci_by_clonotype}}.
#' @param escape_threshold Inactive-haplotype fraction at or above which a
#'   group is flagged as escaping in \code{escapes}. Default 0.1.
#' @param xci_informative_only Logical. If \code{TRUE}, restrict to SNPs that
#'   drove active-X calling in \code{assign_xci} (\code{donor_snp_info$xci_informative}),
#'   excluding escaping or otherwise uninformative SNPs. Default \code{FALSE},
#'   reporting every SNP with a stored phase.
#' @param by_snp Logical. If \code{TRUE}, report every phased SNP rather than
#'   electing one representative per gene, changing the output grain (see Value).
#'   Default \code{FALSE}, i.e. gene-level, which requires gene annotation.
#' @param inverted_phase_genes Character vector of gene names known to be
#'   transcribed predominantly from the \emph{inactive} X, whose stored phase is
#'   therefore expected to be inverted (see \sQuote{Genes masked by phase
#'   inversion}). Matching rows are marked \code{phase_likely_inverted}; nothing
#'   is dropped or corrected. Default \code{"XIST"}. Pass
#'   \code{character(0)} to disable, or extend it with genes you have
#'   independent evidence for -- this is curated prior biology, not a
#'   measurement, and the default is deliberately minimal because \code{XIST} is
#'   the only such gene well established in human somatic cells.
#'
#' @section Gene-level election:
#' Summing a gene's per-SNP counts is not a valid gene total: a read (or a
#' long-read molecule) spanning several of a gene's het SNPs is counted once
#' per SNP, so the sum inflates well-covered multi-SNP genes relative to
#' single-SNP ones. The default therefore elects a single representative SNP per
#' (\code{donor}, \code{gene_name}) and reports only its two groups, so the
#' returned counts remain a genuine read count.
#'
#' The representative is the highest-coverage SNP (summed over both active-X
#' groups) \emph{whose dominant physical allele differs between the two
#' groups}. That flip is what canonical XCI requires -- the expressed allele
#' must change when the active X does -- so requiring it selects a SNP whose
#' data are consistent with the stored phase, rather than letting the
#' best-covered SNP represent the gene regardless of whether it contradicts
#' that phase.
#'
#' Note this test is stricter than \code{!same_allele_dominant}, which is also
#' \code{FALSE} for a SNP with zero coverage in one group: such a SNP
#' demonstrates no flip either way and cannot serve as a representative.
#'
#' It is not, however, a test on the \emph{magnitude} of escape. A SNP whose two
#' groups both exceed 0.5 still flips and remains electable; only a SNP whose
#' groups fall on opposite sides of 0.5 -- i.e. one group flagged
#' \code{phase_contradiction} and the other not -- fails. The criterion is
#' mutual consistency between the two cell populations, not low escape.
#'
#' A (donor, gene) pair where no SNP qualifies has no trustworthy
#' representative and is dropped entirely, with the number dropped reported via
#' \code{logger}; it is not silently backed off to the best-covered SNP.
#'
#' Both \code{same_allele_dominant} and \code{phase_contradiction} are dropped
#' from the elected output: each is uniformly \code{FALSE} by construction,
#' and would read as evidence rather than as the selection criterion it is.
#'
#' @section Genes masked by phase inversion:
#' Because escape above 0.5 is unidentifiable (see Details), a gene transcribed
#' mainly from the inactive X is stored with its phase inverted: the allele it
#' actually expresses is recorded as sitting on the \emph{active} X. Such a gene
#' does not appear as escaping. It appears as one of the cleanest inactivated
#' genes in the output, with \code{escape_fraction} near zero, and its
#' \code{active_count} and \code{inactive_count} are swapped relative to the
#' truth. Verified on \code{XIST}: given its real expression pattern it surfaces
#' here at an \code{escape_fraction} of roughly 0.05.
#'
#' The masking applies to any inactive-X-biased or imprinted gene, not only
#' \code{XIST}, and no measurement in the data distinguishes it -- an inverted
#' gene and a strictly inactivated one produce identical counts. Read-backed
#' phase does not help either: \code{\link{phase_snps}} links SNPs co-observed on
#' one molecule, and a molecule is a single transcript, so phase blocks never
#' span genes and cannot expose a \emph{gene-level} inversion. Resolving it needs
#' external chromosome-level phase, or prior knowledge of the gene.
#'
#' \code{inverted_phase_genes} is that prior knowledge, applied as a label only.
#' Rows whose \code{gene_name} matches are marked \code{phase_likely_inverted};
#' counts, \code{escape_fraction} and the election are all left untouched, since
#' the flag asserts biology rather than anything measured here. For a flagged
#' gene, read \code{active_count} and \code{inactive_count} as reversed and the
#' true escape as roughly \code{1 - escape_fraction}.
#'
#' @section What the default grain omits:
#' Partial escape is \emph{not} excluded: a gene escaping at 0.2, 0.35 or 0.45
#' flips normally, is elected, and has its \code{escape_fraction} reported
#' intact. The election removes only a SNP whose two groups disagree about which
#' side of 0.5 they sit on.
#'
#' The exposure is therefore narrow, but it lands in the worst possible place --
#' on \emph{complete} escapees. Since \code{escape_fraction} saturates at 0.5
#' (see Details), a fully escaping gene sits exactly at the line, both groups
#' land on a random side of it, and the SNP is excluded with probability
#' approaching one half \emph{however deep the coverage}. Away from the line the
#' risk falls off with sampling noise, as \code{1/sqrt(coverage)}, and is
#' negligible by 0.4 at typical depths.
#'
#' The practical failure mode is a maximally escaping gene reported for some
#' donors and silently missing for others, purely by sampling -- which is exactly
#' the pattern that would be misread as a biological difference between donors or
#' experimental groups. When escape is the quantity of interest, always check the
#' \code{logger} count of dropped pairs, inspect them under \code{by_snp = TRUE},
#' and prefer an independent measurement -- read-backed phase via
#' \code{\link{haplotype_expression_by_molecule}} -- to settle whether such a SNP
#' is an escapee or an artefact.
#'
#' @return A tibble with one row per donor, gene, and active-X group -- two rows
#'   per gene -- with columns \code{donor} (when the object carries donor
#'   assignments), \code{snp_id} (the elected representative),
#'   \code{gene_name}, \code{active_x} (the expressed X, "X1" or "X2"),
#'   \code{n_cells} (cells in the group), \code{active_count},
#'   \code{inactive_count}, \code{coverage}, \code{dominant_allele} (the physical
#'   allele, "REF"/"ALT", with more reads in this group; \code{NA} when
#'   uncovered), \code{escape_fraction} (\code{inactive_count / coverage}), and
#'   \code{escapes} (\code{escape_fraction >= escape_threshold}), and
#'   \code{phase_likely_inverted} (\code{TRUE} when \code{gene_name} is in
#'   \code{inverted_phase_genes}). SNPs with no gene annotation (\code{NA}
#'   \code{gene_name}) are excluded, as are genes with no electable SNP.
#'
#'   With \code{by_snp = TRUE} the grain is instead one row per donor, phased
#'   SNP, and active-X group, and two further columns are reported:
#'   \code{phase_contradiction} (per group; \code{TRUE} when
#'   \code{escape_fraction > 0.5}, i.e. the "inactive" haplotype outnumbers the
#'   "active" one -- see Details) and \code{same_allele_dominant} (per SNP;
#'   \code{TRUE} when both covered groups favour the same physical allele). No
#'   group is dropped for contradicting the phase; both flags are reported for
#'   the caller to act on. When the object carries no gene annotation,
#'   \code{gene_name} and \code{phase_likely_inverted} are both omitted.
#'
#' @family X-chromosome inactivation functions
#' @export
#'
#' @examples
#' \dontrun{
#' snp_data <- assign_xci(snp_data)
#'
#' # One elected representative SNP per gene, safe to carry into test_escape()
#' hap <- haplotype_expression(snp_data)
#'
#' # Every SNP, including those no gene could elect -- where escape candidates
#' # and phase artefacts both live
#' hap_by_snp <- haplotype_expression(snp_data, by_snp = TRUE)
#'
#' # SNPs whose dominant allele failed to flip between the active-X groups
#' dplyr::filter(hap_by_snp, same_allele_dominant)
#'
#' # Restrict to the narrower set that drove active-X calling
#' hap_informative <- haplotype_expression(snp_data, xci_informative_only = TRUE)
#' }
setGeneric(
    "haplotype_expression",
    function(
        x,
        escape_threshold = 0.1,
        xci_informative_only = FALSE,
        by_snp = FALSE,
        inverted_phase_genes = "XIST"
    ) {
        standardGeneric("haplotype_expression")
    }
)

#' @rdname haplotype_expression
#' @include SNPData-class.R
setMethod(
    "haplotype_expression",
    signature(x = "SNPData"),
    function(
        x,
        escape_threshold = 0.1,
        xci_informative_only = FALSE,
        by_snp = FALSE,
        inverted_phase_genes = "XIST"
    ) {
        barcode_info <- barcode_info(x)
        snp_info <- snp_info(x)

        if (!"active_x" %in% colnames(barcode_info) || !.has_xci_diagnostics(x)) {
            stop("No stored XCI diagnostics found. Run assign_xci(x) first.")
        }
        # Electing a representative per gene is meaningless without knowing which
        # gene each SNP belongs to. Erroring beats silently falling back to
        # per-SNP output, which would return a different grain than the caller
        # asked for -- the message names the argument that does that on purpose.
        if (!by_snp && !"gene_name" %in% colnames(snp_info)) {
            stop(
                "Gene-level output requires gene annotation; run add_snp_gene_names(x) first, ",
                "or pass by_snp = TRUE for per-SNP output."
            )
        }

        donor_snp_info <- donor_snp_info(x)
        if (xci_informative_only) {
            donor_snp_info <- dplyr::filter(donor_snp_info, xci_informative)
        }
        informative_snp_info <- donor_snp_info %>%
            dplyr::filter(!is.na(allele_on_x1)) %>%
            dplyr::select(snp_id, donor, allele_on_x1)

        # Informative SNPs carry a phase; cells carry a confident active-X call.
        snp_keep <- snp_info$snp_id %in% informative_snp_info$snp_id
        cell_keep <- barcode_info$active_x %in% c("X1", "X2")
        if (!any(snp_keep) || !any(cell_keep)) {
            stop("No phased SNPs or confidently assigned cells available.")
        }

        ref <- ref_count(x)[snp_keep, cell_keep, drop = FALSE]
        alt <- alt_count(x)[snp_keep, cell_keep, drop = FALSE]
        group <- barcode_info$active_x[cell_keep]
        snp_keep_info <- snp_info[snp_keep, ]

        # The XCI fit is per-donor (X1/X2 labels and phase are donor-specific), so
        # counts must be pooled within a donor, never across donors. Absent a donor
        # column, treat the whole object as a single donor.
        has_donor <- "donor" %in% colnames(barcode_info)
        donor <- if (has_donor) as.character(barcode_info$donor[cell_keep]) else rep("all", ncol(ref))

        # Per-donor phase: a SNP retained in several donors carries a distinct
        # X1-allele in each, looked up from donor_snp_info (one row per SNP x
        # donor). Without a donor column the object is single-donor, so each SNP's
        # sole stored phase applies regardless of which donor label it was stored
        # under.
        single_donor_phase <- dplyr::distinct(informative_snp_info, snp_id, .keep_all = TRUE)
        phase_for_donor <- function(d) {
            if (!has_donor) {
                return(single_donor_phase$allele_on_x1[match(snp_keep_info$snp_id, single_donor_phase$snp_id)])
            }
            # allele_on_x1 for one donor, per SNP (NA where the SNP was not
            # informative in that donor, so its counts cannot be phased).
            donor_phase <- dplyr::filter(informative_snp_info, donor == d)
            donor_phase$allele_on_x1[match(snp_keep_info$snp_id, donor_phase$snp_id)]
        }

        # Summarise one donor x active-X group: pool counts across its cells, then
        # split each SNP into active- vs inactive-haplotype reads using the stored phase.
        summarise_group <- function(d, g) {
            cols <- donor == d & group == g
            ref_g <- Matrix::rowSums(ref[, cols, drop = FALSE])
            alt_g <- Matrix::rowSums(alt[, cols, drop = FALSE])

            # Map physical alleles onto haplotypes via this donor's phase.
            allele_on_x1 <- phase_for_donor(d)
            x1_count <- ifelse(allele_on_x1 == "REF", ref_g, alt_g)
            x2_count <- ifelse(allele_on_x1 == "REF", alt_g, ref_g)

            # g names the active haplotype in this group; the other X is silenced.
            active_count <- if (g == "X1") x1_count else x2_count
            inactive_count <- if (g == "X1") x2_count else x1_count
            cov <- active_count + inactive_count

            tibble::tibble(
                donor = d,
                snp_id = snp_keep_info$snp_id,
                gene_name = snp_keep_info$gene_name %||% NA_character_,
                active_x = g,
                n_cells = sum(cols),
                active_count = active_count,
                inactive_count = inactive_count,
                coverage = cov,
                dominant_allele = dplyr::case_when(
                    cov == 0 ~ NA_character_,
                    ref_g >= alt_g ~ "REF",
                    TRUE ~ "ALT"
                )
            )
        }

        combos <- expand.grid(
            donor = sort(unique(donor)),
            active_x = c("X1", "X2"),
            stringsAsFactors = FALSE
        )
        result <- purrr::map2(combos$donor, combos$active_x, summarise_group) %>%
            dplyr::bind_rows() %>%
            # A SNP not informative in this donor has no phase there; its counts
            # cannot be split into active/inactive, so drop those rows.
            dplyr::filter(!is.na(active_count)) %>%
            dplyr::mutate(
                escape_fraction = dplyr::if_else(coverage > 0, inactive_count / coverage, NA_real_),
                escapes = escape_fraction >= escape_threshold
            ) %>%
            # A group where the "inactive" haplotype outnumbers the "active" one
            # contradicts the stored phase rather than showing partial leakage.
            # It is flagged, not dropped: the same signature is produced by a
            # quantification artefact (reference-mapping bias, ambient RNA) and
            # by a near-complete escapee, whose two groups straddle 0.5 by noise
            # and whose side of the line depends on the EM's arbitrary X1/X2
            # labelling. Dropping would make that coin-flip decide whether a real
            # escapee is reported at all.
            #
            # Purely mechanical, with no per-gene exemption. An earlier XIST
            # carve-out was removed: assign_xci() cannot store a phase that makes
            # XIST read above 0.5 in the first place (see phase_likely_inverted
            # below), so the exemption never fired for the case it was written
            # for, and in the one case it could fire -- a near-0.5 straddle -- it
            # suppressed the flag for XIST alone while changing nothing else.
            dplyr::mutate(
                phase_contradiction = !is.na(escape_fraction) & escape_fraction > 0.5
            ) %>%
            # Prior biology, not a measurement: a gene transcribed mainly from
            # the inactive X is stored by assign_xci() with inverted phase and so
            # reports a near-zero escape_fraction, indistinguishable from strict
            # inactivation. Nothing in expression data separates the two, so the
            # only available signal is the gene's identity -- hence a curated
            # list rather than a derived flag.
            dplyr::mutate(phase_likely_inverted = gene_name %in% inverted_phase_genes) %>%
            # The flip must be assessed within a donor: X1/X2 labels are not
            # comparable across donors.
            dplyr::group_by(donor, snp_id) %>%
            # Both covered groups favouring the same physical allele is the XCI
            # contradiction: the active allele failed to flip between them.
            dplyr::mutate(
                same_allele_dominant = dplyr::n_distinct(dominant_allele[coverage > 0]) == 1L &
                    sum(coverage > 0) == 2L
            ) %>%
            dplyr::ungroup()

        # Elect before the donor/gene_name columns are conditionally dropped
        # below: the election is per (donor, gene) and needs both.
        if (!by_snp) {
            result <- .elect_gene_representative_snps(result)
        }

        if (!has_donor) {
            result <- dplyr::select(result, -donor)
        }

        # Drop the placeholder gene_name column when the object carries no
        # annotation, and with it the curated flag, which is keyed on gene name
        # and would otherwise be a uniformly FALSE column asserting that no gene
        # is inverted when in truth none could be checked.
        if (!"gene_name" %in% colnames(snp_info)) {
            result <- dplyr::select(result, -gene_name, -phase_likely_inverted)
        }
        result
    }
)

#' Elect one representative SNP per donor and gene
#'
#' Collapses \code{\link{haplotype_expression}}'s per-SNP rows to the
#' highest-coverage SNP whose dominant physical allele flips between the two
#' active-X groups. Summing a gene's SNPs instead would count a read spanning
#' several of them once per SNP; electing one keeps the output a genuine read
#' count. See the \sQuote{Collapsing to the gene level} section of
#' \code{\link{haplotype_expression}} for why the flip is required.
#'
#' @param result The per-SNP tibble assembled by \code{\link{haplotype_expression}},
#'   still carrying its \code{donor}, \code{gene_name} and
#'   \code{same_allele_dominant} columns.
#'
#' @return \code{result} restricted to the elected SNPs' rows, without
#'   \code{same_allele_dominant}, ordered by donor, gene and active X.
#'
#' @keywords internal
.elect_gene_representative_snps <- function(result) {
    annotated <- dplyr::filter(result, !is.na(gene_name))

    # Both groups covered AND favouring different physical alleles. Deliberately
    # not !same_allele_dominant: that is also FALSE for a SNP covered in only one
    # group, which demonstrates no flip and so cannot represent the gene. Total
    # coverage is summed over both groups, so the ranking is not decided by
    # whichever group happens to hold more cells.
    #
    # This is deliberately not a test on the magnitude of escape. If both groups
    # exceed 0.5, each is dominated by the other X's allele -- two different
    # physical alleles -- so the pair still flips and stays electable. Only groups
    # on opposite sides of 0.5 fail, which is mutual inconsistency between the two
    # cell populations rather than high escape. Symmetric escape at any level
    # therefore survives; see the "What the default grain omits" section.
    candidates <- annotated %>%
        dplyr::group_by(donor, snp_id) %>%
        dplyr::filter(
            sum(coverage > 0) == 2L,
            dplyr::n_distinct(dominant_allele[coverage > 0]) == 2L
        ) %>%
        dplyr::mutate(snp_coverage = sum(coverage)) %>%
        dplyr::ungroup()

    # Ties on coverage are broken by snp_id so repeat runs elect the same SNP;
    # slice_max(with_ties = FALSE) alone would depend on incoming row order.
    elected <- candidates %>%
        dplyr::distinct(donor, gene_name, snp_id, snp_coverage) %>%
        dplyr::arrange(donor, gene_name, dplyr::desc(snp_coverage), snp_id) %>%
        dplyr::slice_head(n = 1, by = c(donor, gene_name))

    n_before <- nrow(dplyr::distinct(annotated, donor, gene_name))
    n_dropped <- n_before - nrow(elected)
    if (n_dropped > 0) {
        logger::log_warn(
            "Dropped {n_dropped}/{n_before} donor x gene pair(s) from gene-level output: ",
            "no SNP had its dominant allele flip between the active-X groups"
        )
    }

    # same_allele_dominant and phase_contradiction are both uniformly FALSE for
    # every elected SNP, by the argument above; keeping them would read as
    # evidence rather than as the selection criteria they are.
    candidates %>%
        dplyr::semi_join(elected, by = c("donor", "gene_name", "snp_id")) %>%
        dplyr::select(-snp_coverage, -same_allele_dominant, -phase_contradiction) %>%
        dplyr::arrange(donor, gene_name, active_x)
}

#' Phased active/inactive haplotype expression per gene, pooling molecules across SNPs
#'
#' Counts each molecule once per gene regardless of how many of the gene's
#' heterozygous SNPs it covers, unlike \code{\link{haplotype_expression}}
#' (which reports one row per SNP): summing that function's per-SNP counts
#' back up to a gene total would recreate the multi-SNP-per-molecule
#' over-count read-backed phasing exists to remove, since a molecule
#' spanning several of a gene's SNPs would cast one independent vote per row
#' again. Pooling instead at the molecule level, within the gene's
#' best-supported phase block, is what actually increases usable evidence
#' per gene -- the motivation being genes like escapees with many SNPs
#' (e.g. PRKX, 292 SNPs) that \code{assign_xci}'s EM reduces to a single
#' representative SNP.
#'
#' @details
#' A gene's heterozygous SNPs can fall into more than one read-backed phase
#' block (\code{\link{phase_snps}}) when no single molecule spans all of
#' them -- most commonly a spliced/unspliced boundary, since a mature mRNA
#' and its unspliced precursor rarely share a molecule. For each
#' (\code{gene_name}, \code{donor}), only the block backed by the most
#' distinct molecules is pooled; molecules whose SNPs fall in any other
#' block for the same gene are real evidence but cannot be pooled with the
#' dominant block's haplotype labels (their phase has no established
#' relationship to it), so they are reported as \code{n_stranded_molecules}
#' rather than dropped silently.
#'
#' A molecule's haplotype is called by majority vote across the SNPs it
#' covers within the dominant block only (a tie is base-calling noise, once
#' phase is accounted for, and is dropped as \code{"ambiguous"}). Only genes
#' with a stored, resolved \code{allele_on_x1} and \code{phase_block} (i.e.
#' processed by \code{\link{add_molecule_phase}}) contribute. A gene whose
#' only heterozygous SNP could not be linked to any other by
#' \code{\link{phase_snps}} still contributes here if that SNP already has
#' an EM-derived phase from \code{\link{assign_xci}}: \code{\link{add_molecule_phase}}
#' gives it its own singleton block, so a single-SNP gene benefits from
#' correct molecule-level (rather than read-level) counting too. Only a
#' single-SNP gene with \emph{no} EM-derived phase at all (never selected by
#' \code{assign_xci}'s per-gene SNP pick) has no way to be oriented and is
#' left to \code{haplotype_expression}, which already handles a single SNP
#' correctly.
#'
#' @inheritSection assign_xci Phase is inferred from expression, not genotyped
#'
#' @param x A SNPData object that has had XCI diagnostics stored by
#'   \code{\link{assign_xci}} (or \code{\link{assign_xci_by_clonotype}}) and
#'   subsequently had read-backed phase added by
#'   \code{\link{add_molecule_phase}}.
#' @param molecule_calls A tibble with columns \code{donor}, \code{barcode},
#'   \code{umi}, \code{snp_id}, \code{allele}, \code{transcript_strand} --
#'   the union, across every donor, of \code{\link{add_molecule_phase}}'s
#'   \code{"molecule_calls"} attribute (or \code{\link{molecule_snp_alleles}}
#'   plus a \code{transcript_strand} column from
#'   \code{\link{molecule_read_strand}}, bind rows and add a \code{donor}
#'   column).
#' @param snp_gene_map A tibble with columns \code{snp_id}, \code{gene_name},
#'   \code{gene_strand}, \code{ambiguous}, as returned by
#'   \code{\link{assign_snp_genes}}. A SNP with \code{ambiguous = TRUE}
#'   overlaps more than one gene and is only counted towards a given
#'   candidate's \code{gene_name} for molecules whose own
#'   \code{transcript_strand} matches that candidate's \code{gene_strand};
#'   molecules with no strand information (\code{NA}) or the wrong strand are
#'   dropped for that SNP rather than guessed. \code{ambiguous = FALSE} rows
#'   apply regardless of strand.
#' @param escape_threshold Inactive-haplotype fraction at or above which a
#'   group is flagged as escaping. Default 0.1.
#'
#' @return A tibble with one row per (\code{gene_name}, \code{donor},
#'   \code{active_x}) triple, with columns \code{donor}, \code{gene_name},
#'   \code{active_x} (the expressed X, "X1" or "X2"), \code{active_count},
#'   \code{inactive_count}, \code{coverage}, \code{escape_fraction}
#'   (\code{inactive_count / coverage}), \code{escapes}
#'   (\code{escape_fraction >= escape_threshold}), \code{phase_block_used}
#'   (the phase block whose molecules were pooled), \code{dominant_molecules}
#'   (molecules backing the pooled block), and \code{n_stranded_molecules}
#'   (molecules of the same gene in a different, unpooled block).
#'
#' @family X-chromosome inactivation functions
#' @export
#'
#' @examples
#' \dontrun{
#' snp_data <- assign_xci(snp_data)
#' snp_data <- add_molecule_phase(snp_data, bam_files = c(donor0 = "donor0.bam"))
#' hap <- haplotype_expression_by_molecule(snp_data, molecule_calls, snp_gene_map)
#' }
setGeneric(
    "haplotype_expression_by_molecule",
    function(x, molecule_calls, snp_gene_map, escape_threshold = 0.1) {
        standardGeneric("haplotype_expression_by_molecule")
    }
)

#' @rdname haplotype_expression_by_molecule
#' @include SNPData-class.R
setMethod(
    "haplotype_expression_by_molecule",
    signature(x = "SNPData"),
    function(x, molecule_calls, snp_gene_map, escape_threshold = 0.1) {
        barcode_info <- barcode_info(x)
        donor_snp_info <- donor_snp_info(x)

        if (!"active_x" %in% colnames(barcode_info) || !.has_xci_diagnostics(x)) {
            stop("No stored XCI diagnostics found. Run assign_xci(x) first.")
        }
        if (!all(c("phase_block", "allele_on_x1") %in% colnames(donor_snp_info))) {
            stop("No stored molecule phase found. Run add_molecule_phase(x) first.")
        }
        required_call_cols <- c("donor", "barcode", "umi", "snp_id", "allele", "transcript_strand")
        missing_call_cols <- setdiff(required_call_cols, colnames(molecule_calls))
        if (length(missing_call_cols) > 0) {
            stop("molecule_calls is missing required column(s): ", paste(missing_call_cols, collapse = ", "))
        }
        required_gene_cols <- c("snp_id", "gene_name", "gene_strand", "ambiguous")
        missing_gene_cols <- setdiff(required_gene_cols, colnames(snp_gene_map))
        if (length(missing_gene_cols) > 0) {
            stop("snp_gene_map is missing required column(s): ", paste(missing_gene_cols, collapse = ", "))
        }

        phase <- donor_snp_info %>%
            dplyr::filter(!is.na(allele_on_x1), !is.na(phase_block)) %>%
            dplyr::select(snp_id, donor, allele_on_x1, phase_block)

        # is_x1: whether this molecule's allele at this SNP is the one the
        # resolved phase names as sitting on X1 -- the gene-level analogue of
        # haplotype_expression()'s per-SNP active/inactive split. A SNP
        # assign_snp_genes() flagged ambiguous only counts towards a candidate
        # gene for molecules whose own transcript strand matches that
        # candidate's strand; strand-unknown or mismatched molecules are
        # dropped for that SNP rather than guessed. Joining an ambiguous SNP's
        # multiple molecules against its multiple gene candidates is a genuine
        # many-to-many fan-out, narrowed back down by the strand filter below.
        calls <- molecule_calls %>%
            dplyr::inner_join(phase, by = c("snp_id", "donor")) %>%
            dplyr::inner_join(snp_gene_map, by = "snp_id", relationship = "many-to-many") %>%
            dplyr::filter(!ambiguous | (!is.na(transcript_strand) & transcript_strand == gene_strand)) %>%
            dplyr::mutate(is_x1 = allele == allele_on_x1)

        if (nrow(calls) == 0) {
            stop("No molecule calls could be matched to a phased, singly-mapped-gene SNP.")
        }

        blocks <- calls %>%
            dplyr::distinct(donor, gene_name, phase_block, barcode, umi) %>%
            dplyr::count(donor, gene_name, phase_block, name = "molecules")
        best_block <- blocks %>%
            dplyr::slice_max(molecules, n = 1, by = c(donor, gene_name), with_ties = FALSE) %>%
            dplyr::select(donor, gene_name, phase_block, dominant_molecules = molecules)
        stranded <- blocks %>%
            dplyr::anti_join(best_block, by = c("donor", "gene_name", "phase_block")) %>%
            dplyr::summarise(n_stranded_molecules = sum(molecules), .by = c(donor, gene_name))

        # A molecule votes across every SNP it covers within the dominant block
        # only; majority wins, a tie is ambiguous (residual base-calling noise
        # once phase is accounted for) and dropped rather than guessed.
        molecules <- calls %>%
            dplyr::semi_join(best_block, by = c("donor", "gene_name", "phase_block")) %>%
            dplyr::summarise(n_x1 = sum(is_x1), n_x2 = sum(!is_x1), .by = c(donor, gene_name, barcode, umi)) %>%
            dplyr::mutate(
                haplotype = dplyr::case_when(n_x1 > n_x2 ~ "X1", n_x2 > n_x1 ~ "X2", TRUE ~ "ambiguous")
            ) %>%
            dplyr::filter(haplotype != "ambiguous")

        cell_groups <- barcode_info %>%
            dplyr::filter(active_x %in% c("X1", "X2")) %>%
            dplyr::select(barcode, donor, active_x)

        # Pool molecules per (donor, gene, active-X group): active_count is
        # molecules landing on the haplotype the group's own active_x names,
        # inactive_count the other -- the same convention as
        # haplotype_expression()'s summarise_group().
        counts <- molecules %>%
            dplyr::inner_join(cell_groups, by = c("donor", "barcode")) %>%
            dplyr::summarise(
                n_x1 = sum(haplotype == "X1"),
                n_x2 = sum(haplotype == "X2"),
                .by = c(donor, gene_name, active_x)
            ) %>%
            dplyr::mutate(
                active_count = dplyr::if_else(active_x == "X1", n_x1, n_x2),
                inactive_count = dplyr::if_else(active_x == "X1", n_x2, n_x1),
                coverage = active_count + inactive_count
            ) %>%
            dplyr::select(donor, gene_name, active_x, active_count, inactive_count, coverage)

        # Every (donor, gene) gets both active_x groups reported, even one with
        # zero molecules, so a completely silenced haplotype reads as coverage
        # zero rather than being silently absent from the output.
        tidyr::expand_grid(
            dplyr::distinct(counts, donor, gene_name),
            active_x = c("X1", "X2")
        ) %>%
            dplyr::left_join(counts, by = c("donor", "gene_name", "active_x")) %>%
            dplyr::mutate(
                active_count = dplyr::coalesce(active_count, 0L),
                inactive_count = dplyr::coalesce(inactive_count, 0L),
                coverage = dplyr::coalesce(coverage, 0L),
                escape_fraction = dplyr::if_else(coverage > 0, inactive_count / coverage, NA_real_),
                escapes = escape_fraction >= escape_threshold
            ) %>%
            dplyr::left_join(best_block, by = c("donor", "gene_name")) %>%
            dplyr::left_join(stranded, by = c("donor", "gene_name")) %>%
            dplyr::mutate(n_stranded_molecules = dplyr::coalesce(n_stranded_molecules, 0L)) %>%
            dplyr::rename(phase_block_used = phase_block) %>%
            dplyr::arrange(donor, gene_name, active_x)
    }
)

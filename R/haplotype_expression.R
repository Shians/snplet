#' Phased active/inactive haplotype expression per gene
#'
#' Splits phased het SNPs' ALT/REF counts into active- and inactive-haplotype
#' reads using the per-SNP phase that \code{\link{assign_xci}} assigned,
#' computed separately within the X1-active and X2-active cell groups it
#' called, then reports one selected representative SNP per donor and gene (see
#' \sQuote{Gene-level representative selection}), with the two groups pooled
#' into a single row. \code{by_snp = TRUE} reports every SNP instead, useful
#' for working out why a gene is missing from the gene-level output, or for
#' double-checking one whose result looks wrong; \code{by_active_x = TRUE} keeps
#' the two active-X groups on separate rows.
#'
#' Unlike a per-group \code{pmin}/\code{pmax} on raw counts, the
#' allele-to-haplotype assignment comes from the stored per-donor phase
#' (\code{donor_snp_info$allele_on_x1}) and is applied consistently to both
#' groups, so the same physical allele is never relabelled as "active" in
#' both groups; that contradiction is instead flagged explicitly.
#'
#' Both the gene-level and per-SNP output span every donor-genotyped het SNP with a stored phase, including
#' SNPs \code{assign_xci} excluded from active-X calling as uninformative or
#' escaping, unless \code{xci_informative_only = TRUE}.
#'
#' @details
#' For each SNP and donor, the stored phase assigns REF or ALT to X1 (see
#' \sQuote{Phase is inferred from expression, not genotyped}); the other
#' allele is X2. Within each active-X group the expected-active haplotype (X1
#' in the X1-active group, X2 in the X2-active group) determines how counts are
#' split into \code{active_count}/\code{inactive_count}. The two groups are
#' then summed into one row per gene unless \code{by_active_x = TRUE}: they are
#' disjoint sets of cells, so the sum double-counts nothing, and the pooled
#' \code{escape_fraction} is better powered than either group's own. Pooling
#' also matches the grain \code{\link{test_escape}} tests at and that
#' \code{donor_info(x)$xci_rho} is fitted at, so the default output can be
#' tested directly.
#'
#' Columns meaningful only within a single group are dropped when pooling:
#' \code{active_x}, \code{dominant_allele} (which flips between the groups by
#' construction) and \code{phase_contradiction}. Pass \code{by_active_x = TRUE}
#' to inspect them, or to compare the two cell populations against each other.
#'
#' \code{same_allele_dominant} summarises both groups rather than describing
#' one, so it survives pooling and is reported under \code{by_snp = TRUE}
#' whichever way \code{by_active_x} is set. It is \code{TRUE} when the same
#' physical allele dominates in both groups, flagging escape or
#' XCI-independent imbalance
#' rather than hiding it. Under canonical XCI the dominant allele should flip
#' between groups (X1's allele in X1-active cells, X2's in X2-active cells),
#' so a failure to flip is what triggers the flag. \code{escape_fraction}
#' gives the magnitude, as the fraction of reads from the haplotype that
#' should be silenced.
#'
#' \strong{\code{escape_fraction} is bounded between 0 and 0.5.} Expression
#' alone cannot distinguish a gene transcribed mostly from the inactive X
#' from one transcribed mostly from the active X with inverted phase, so
#' \code{\link{assign_xci}} always names the majority allele active. The
#' scale therefore runs from 0 (strict inactivation) to 0.5 (complete
#' escape), not to 1. Set \code{escape_threshold}, and any null passed to
#' \code{\link{test_escape}}, on this 0-0.5 scale.
#'
#' One consequence: a gene expressed mainly from the inactive X is stored
#' with inverted phase and reads as unusually \emph{clean} XCI rather than as
#' escaping; see \sQuote{Genes masked by phase inversion}.
#'
#' \code{phase_contradiction} is \code{TRUE} for a group whose own
#' \code{escape_fraction} exceeds 0.5, even though \code{\link{assign_xci}}'s
#' donor-pooled estimate for the gene is capped at 0.5: each active-X group
#' here draws on fewer cells than that pooled fit, so its own estimate is
#' noisier and can cross 0.5 by chance alone, especially when the gene's true
#' escape is near 0.5. Such rows are kept rather than dropped, since the same
#' pattern can come from a quantification artefact or from a genuine complete
#' escapee, and dropping it would let chance decide whether real escapees are
#' reported.
#'
#' \code{same_allele_dominant} is \code{TRUE} exactly when the two groups'
#' \code{phase_contradiction} status disagrees (with both groups covered):
#' the two flags are related but not equivalent.
#'
#' @inheritSection assign_xci Phase is inferred from expression, not genotyped
#'
#' @param x A SNPData object, required, that had XCI diagnostics stored by
#'   \code{\link{assign_xci}} or \code{\link{assign_xci_by_clonotype}}.
#' @param escape_threshold Numeric, in \code{[0, 1]} (default 0.1).
#'   Inactive-haplotype fraction at or above which a row is flagged as
#'   escaping in \code{escapes}, applied to the pooled \code{escape_fraction}
#'   unless \code{by_active_x = TRUE}.
#' @param xci_informative_only Logical (default \code{FALSE}, reporting every
#'   SNP with a stored phase). If \code{TRUE}, restrict to SNPs that drove
#'   active-X calling in \code{assign_xci}
#'   (\code{donor_snp_info$xci_informative}), excluding escaping or otherwise
#'   uninformative SNPs.
#' @param by_snp Logical (default \code{FALSE}, i.e. gene-level, which
#'   requires gene annotation). If \code{TRUE}, report every phased SNP
#'   rather than selecting one representative per gene, changing what each
#'   output row represents (see Value).
#' @param by_active_x Logical (default \code{FALSE}, i.e. the X1-active and
#'   X2-active cell groups pooled into one row). If \code{TRUE}, report the
#'   two groups on separate rows, adding \code{active_x},
#'   \code{dominant_allele} and \code{phase_contradiction}. Independent of
#'   \code{by_snp}: it changes how the cells are split, not whether rows are
#'   SNPs or genes.
#' @param inverted_phase_genes Character vector (default \code{"XIST"}).
#'   Gene names known to be transcribed predominantly from the
#'   \emph{inactive} X, whose stored phase is therefore expected to be
#'   inverted (see \sQuote{Genes masked by phase inversion}). Matching rows
#'   are marked \code{phase_likely_inverted}; nothing is dropped or
#'   corrected. Pass \code{character(0)} to disable, or extend it with genes
#'   you have independent evidence for; this is curated prior biology, not a
#'   measurement, and the default is deliberately minimal because \code{XIST}
#'   is the only such gene well established in human somatic cells.
#'
#' @section Gene-level representative selection:
#' The default selects a single representative SNP per (\code{donor},
#' \code{gene_name}) rather than summing per-SNP counts, since summing is not
#' a valid gene total: a read (or long-read molecule) spanning several of a
#' gene's het SNPs would be counted once per SNP, inflating well-covered
#' multi-SNP genes relative to single-SNP ones. Only the representative's two
#' groups are reported, keeping the returned counts a genuine read count.
#'
#' The representative is the highest-coverage SNP (summed over both active-X
#' groups) whose dominant physical allele differs between the two groups,
#' since that flip is what canonical XCI requires: the expressed allele must
#' change when the active X does. This selects a SNP consistent with the
#' stored phase, rather than simply the best-covered SNP regardless of
#' whether it contradicts that phase.
#'
#' This is stricter than \code{!same_allele_dominant} (also \code{FALSE} for
#' zero coverage in one group, which demonstrates no flip either way) and is
#' not a test on escape \emph{magnitude}: a SNP whose two groups both exceed
#' 0.5 still flips and remains a valid representative, failing only when the
#' groups fall on opposite sides of 0.5 (one flagged
#' \code{phase_contradiction}, the other not). The criterion is mutual
#' consistency between the two cell populations, not low escape.
#'
#' A (donor, gene) pair with no qualifying SNP is dropped entirely, with the
#' number dropped reported via \code{logger}, rather than falling back to the
#' best-covered SNP.
#'
#' \code{same_allele_dominant} and \code{phase_contradiction} are dropped
#' from the selected output, since both are uniformly \code{FALSE} for every
#' retained row and would read as evidence rather than as the selection
#' criterion they are.
#'
#' @section Genes masked by phase inversion:
#' A gene transcribed mainly from the inactive X is stored with its phase
#' inverted (the allele it actually expresses is recorded as sitting on the
#' \emph{active} X), because escape above 0.5 is unidentifiable (see
#' Details). It does not appear as escaping: it reads as one of the cleanest
#' inactivated genes in the output, with \code{escape_fraction} near zero and
#' \code{active_count}/\code{inactive_count} swapped relative to the truth.
#'
#' The masking applies to any inactive-X-biased or imprinted gene, not only
#' \code{XIST}: no measurement distinguishes an inverted gene from a strictly
#' inactivated one, since both produce identical counts. Read-backed phase
#' cannot help either, since \code{\link{phase_snps}} links SNPs co-observed
#' on one molecule (a single transcript), so phase blocks never span genes
#' and cannot expose a \emph{gene-level} inversion; resolving it needs
#' external chromosome-level phase, or prior knowledge of the gene.
#'
#' \code{inverted_phase_genes} supplies that prior knowledge, applied as a
#' label only: matching rows are marked \code{phase_likely_inverted}, with
#' counts, \code{escape_fraction} and the representative selection left
#' untouched, since the flag asserts biology rather than anything measured
#' here. For a flagged gene, read \code{active_count} and
#' \code{inactive_count} as reversed and the true escape as roughly
#' \code{1 - escape_fraction}.
#'
#' @section What the default output omits:
#' Partial escape is \emph{not} excluded: a gene escaping at 0.2, 0.35 or 0.45
#' flips normally, is selected as a representative, and has its
#' \code{escape_fraction} reported intact. The selection step removes only a
#' SNP whose two groups disagree about which side of 0.5 they sit on.
#'
#' The exposure lands specifically on \emph{complete} escapees, since
#' \code{escape_fraction} saturates at 0.5 (see Details): a fully escaping
#' gene sits exactly on that line, so each SNP's two groups land on a random
#' side of it and the SNP is excluded with probability approaching one half
#' \emph{however deep the coverage}. Away from the line the risk falls off
#' with sampling noise, as \code{1/sqrt(coverage)}, and is negligible by 0.4
#' at typical depths.
#'
#' The practical failure mode is a maximally escaping gene reported for some
#' donors and missing for others without any flag, purely by sampling, which
#' is exactly the pattern that would be misread as a biological difference
#' between donors or experimental groups. When escape is the quantity of
#' interest, always check the \code{logger} count of dropped pairs, inspect
#' them under \code{by_snp = TRUE}, and prefer an independent measurement,
#' read-backed phase via \code{\link{haplotype_expression_by_molecule}}, to
#' settle whether such a SNP is an escapee or an artefact.
#'
#' @return A tibble with one row per donor and gene, with columns \code{donor}
#'   (when the object carries donor assignments), \code{snp_id} (the selected
#'   representative), \code{gene_name}, \code{n_cells} (cells contributing),
#'   \code{active_count}, \code{inactive_count}, \code{coverage},
#'   \code{escape_fraction} (\code{inactive_count / coverage}), \code{escapes}
#'   (\code{escape_fraction >= escape_threshold}) and
#'   \code{phase_likely_inverted} (\code{TRUE} when \code{gene_name} is in
#'   \code{inverted_phase_genes}). SNPs with no gene annotation (\code{NA}
#'   \code{gene_name}) are excluded, as are genes with no qualifying SNP. This
#'   is the grain \code{\link{test_escape}} expects, so the result can be
#'   passed to it directly.
#'
#'   With \code{by_snp = TRUE} each row instead represents one donor and phased
#'   SNP, with \code{same_allele_dominant} added (per SNP; \code{TRUE} when
#'   both covered groups favour the same physical allele).
#'
#'   With \code{by_active_x = TRUE} each row is split in two, one per active-X
#'   group, and three further columns are reported: \code{active_x} (the
#'   expressed X, "X1" or "X2"), \code{dominant_allele} (the physical allele,
#'   "REF"/"ALT", with more reads in this group; \code{NA} when uncovered) and
#'   \code{phase_contradiction} (per group; \code{TRUE} when
#'   \code{escape_fraction > 0.5}, i.e. the "inactive" haplotype outnumbers the
#'   "active" one, see Details). No group is dropped for contradicting the
#'   phase; the flags are reported for you to act on.
#'
#'   When the object carries no gene annotation, \code{gene_name} and
#'   \code{phase_likely_inverted} are both omitted.
#'
#' @family X-chromosome inactivation functions
#' @export
#'
#' @examples
#' \dontrun{
#' # Requires a snp_data object already processed by assign_xci() (see its
#' # examples), which in turn needs enough X-linked heterozygous SNP coverage
#' # to fit the model.
#' snp_data <- assign_xci(snp_data)
#'
#' # One representative SNP per gene, safe to carry into test_escape()
#' hap <- haplotype_expression(snp_data)
#'
#' # Every SNP, including those no gene could select as its representative,
#' # where escape candidates and phase artefacts both live
#' hap_by_snp <- haplotype_expression(snp_data, by_snp = TRUE)
#'
#' # SNPs whose dominant allele failed to flip between the active-X groups
#' dplyr::filter(hap_by_snp, same_allele_dominant)
#'
#' # The X1-active and X2-active cell populations kept apart, for comparing
#' # them against each other rather than testing the gene
#' hap_by_active_x <- haplotype_expression(snp_data, by_active_x = TRUE)
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
        by_active_x = FALSE,
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
        by_active_x = FALSE,
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

        # Both the election above and the flip diagnostics need the two groups
        # apart, so pooling is the last step rather than the first.
        if (!by_active_x) {
            result <- .pool_active_x_groups(result, escape_threshold)
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

#' Select one representative SNP per donor and gene
#'
#' Collapses \code{\link{haplotype_expression}}'s per-SNP rows to the
#' highest-coverage SNP whose dominant physical allele flips between the two
#' active-X groups. Summing a gene's SNPs instead would count a read spanning
#' several of them once per SNP; selecting one keeps the output a genuine read
#' count. See the \sQuote{Gene-level representative selection} section of
#' \code{\link{haplotype_expression}} for why the flip is required.
#'
#' @param result The per-SNP tibble assembled by \code{\link{haplotype_expression}},
#'   still carrying its \code{donor}, \code{gene_name} and
#'   \code{same_allele_dominant} columns.
#'
#' @return \code{result} restricted to the selected SNPs' rows, without
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
    # therefore survives; see the "What the default output omits" section.
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

#' Pool the two active-X groups into one row
#'
#' Sums \code{active_count}/\code{inactive_count} across the X1-active and
#' X2-active groups, giving one row per gene (or per SNP under
#' \code{by_snp = TRUE}). The two groups are disjoint sets of cells, so the sum
#' double-counts nothing; it is the same pooling
#' \code{\link{.fit_pooled_rho_by_donor}} fits \code{xci_rho} against.
#'
#' Columns that only mean something within a single group are dropped:
#' \code{active_x} itself, \code{dominant_allele} (which flips between the
#' groups by construction, so a pooled value would be an artefact of whichever
#' group had more reads) and \code{phase_contradiction} (defined per group).
#' \code{same_allele_dominant} is already a per-SNP summary of both groups and
#' survives.
#'
#' @param result A per-group tibble assembled by
#'   \code{\link{haplotype_expression}}, after any representative election.
#' @param escape_threshold Numeric, in \code{[0, 1]}, required. Applied to the
#'   pooled \code{escape_fraction}, which is a better-powered estimate than
#'   either group's own.
#'
#' @return \code{result} with one row per donor and SNP, ordered by donor and
#'   gene where those columns are present.
#'
#' @keywords internal
.pool_active_x_groups <- function(result, escape_threshold) {
    id_cols <- intersect(
        c("donor", "snp_id", "gene_name", "phase_likely_inverted", "same_allele_dominant"),
        colnames(result)
    )

    pooled <- result %>%
        dplyr::summarise(
            n_cells = sum(n_cells),
            active_count = sum(active_count),
            inactive_count = sum(inactive_count),
            .by = dplyr::all_of(id_cols)
        ) %>%
        dplyr::mutate(
            coverage = active_count + inactive_count,
            escape_fraction = dplyr::if_else(coverage > 0, inactive_count / coverage, NA_real_),
            escapes = escape_fraction >= escape_threshold
        ) %>%
        # The two flags are per-SNP metadata rather than measurements, so they
        # trail the counts as they do in the per-group output.
        dplyr::relocate(dplyr::any_of(c("phase_likely_inverted", "same_allele_dominant")), .after = dplyr::last_col())

    sort_cols <- intersect(c("donor", "gene_name", "snp_id"), colnames(pooled))
    dplyr::arrange(pooled, dplyr::pick(dplyr::all_of(sort_cols)))
}

#' Phased active/inactive haplotype expression per gene, pooling molecules across SNPs
#'
#' Counts each molecule once per gene regardless of how many of the gene's
#' heterozygous SNPs it covers, unlike \code{\link{haplotype_expression}}
#' (one row per SNP). Summing that function's per-SNP counts to a gene total
#' would recreate the multi-SNP-per-molecule over-count read-backed phasing
#' exists to remove, since a molecule spanning several SNPs would cast one
#' vote per row. Pooling instead at the molecule level, within the gene's
#' best-supported phase block, increases usable evidence per gene,
#' particularly for escapees with many SNPs (e.g. PRKX, 292 SNPs) that
#' \code{assign_xci}'s EM reduces to a single representative SNP.
#'
#' @details
#' A gene's heterozygous SNPs can fall into more than one read-backed phase
#' block (\code{\link{phase_snps}}) when no single molecule spans all of
#' them, most commonly at a spliced/unspliced boundary (a mature mRNA and its
#' unspliced precursor rarely share a molecule). For each (\code{gene_name},
#' \code{donor}), only the block backed by the most distinct molecules is
#' pooled; molecules whose SNPs fall in any other block are real evidence but
#' cannot be pooled with the dominant block's haplotype labels, so they are
#' reported as \code{n_stranded_molecules} rather than dropped without record.
#'
#' A molecule's haplotype is called by majority vote across the SNPs it
#' covers within the dominant block only; a tie is dropped as
#' \code{"ambiguous"} (residual base-calling noise once phase is accounted
#' for). Only genes with a stored, resolved \code{allele_on_x1} and
#' \code{phase_block} (i.e. processed by \code{\link{add_molecule_phase}})
#' contribute. A single-SNP gene contributes too if that SNP already has an
#' EM-derived phase from \code{\link{assign_xci}}: \code{add_molecule_phase()}
#' gives it its own block of one, extending correct molecule-level counting
#' to single-SNP genes. Only a single-SNP gene with \emph{no} EM-derived
#' phase at all has no way to be oriented, and is left to
#' \code{haplotype_expression}, which already handles it correctly.
#'
#' @inheritSection assign_xci Phase is inferred from expression, not genotyped
#'
#' @param x A SNPData object, required, that has had XCI diagnostics stored
#'   by \code{\link{assign_xci}} (or \code{\link{assign_xci_by_clonotype}})
#'   and subsequently had read-backed phase added by
#'   \code{\link{add_molecule_phase}}, which is also where the per-molecule
#'   allele calls come from (see \sQuote{Where the molecule calls come from}).
#' @param escape_threshold Numeric, in \code{[0, 1]} (default 0.1).
#'   Inactive-haplotype fraction at or above which a row is flagged as
#'   escaping.
#' @param by_active_x Logical (default \code{FALSE}, i.e. the X1-active and
#'   X2-active cell groups pooled into one row per gene, matching
#'   \code{\link{haplotype_expression}} and the grain
#'   \code{\link{test_escape}} tests at). If \code{TRUE}, report the two
#'   groups on separate rows, adding \code{active_x}.
#' @param pool_blocks Logical (default \code{TRUE}). If \code{TRUE}, count
#'   molecules from every one of a gene's phase blocks. If \code{FALSE}, count
#'   only the gene's largest block, leaving the rest reported but uncounted in
#'   \code{n_stranded_molecules}. See \sQuote{Pooling a gene's phase blocks}.
#'
#' @section Pooling a gene's phase blocks:
#' \code{\link{phase_snps}} cannot link two SNPs no single molecule spans, so
#' one gene routinely ends up split across several phase blocks. Every SNP
#' counted here nonetheless carries a globally oriented \code{allele_on_x1} --
#' either the EM's own call, or an orientation \code{\link{add_molecule_phase}}
#' propagated to the block from an EM anchor -- so a molecule's haplotype means
#' the same thing in every block, and pooling blocks is arithmetically sound
#' rather than a mixing of incompatible labels.
#'
#' What pooling stakes is that each block's orientation is individually
#' correct. Within a block, relative phase is physical: two SNPs seen on one
#' molecule sit on one chromosome, which is observed rather than inferred.
#' Across blocks it is inherited from the EM anchors, which are
#' expression-derived, so a block oriented backwards contributes its molecules
#' to the wrong haplotype and inflates \code{escape_fraction}. Counting only
#' the largest block (\code{pool_blocks = FALSE}) makes that impossible, at the
#' cost of discarding the smaller blocks entirely.
#'
#' Because \code{escape_fraction} is bounded above by 0.5 by construction (see
#' \code{\link{haplotype_expression}}), a block whose own inactive fraction
#' exceeds 0.5 is not an extreme escapee but a block oriented the wrong way
#' round. The molecules in such blocks are counted, but reported separately as
#' \code{discordant_block_molecules} so the one failure mode pooling introduces
#' is visible rather than silently averaged in.
#'
#' @section Where the inputs come from:
#' Beyond \code{x}, there is nothing to supply. Both inputs this function needs
#' are already on the object:
#'
#' \describe{
#'   \item{The per-molecule allele calls}{Read from the
#'     \code{"molecule_calls"} attribute \code{\link{add_molecule_phase}} left
#'     there when it extracted them from the BAM files. They are an attribute
#'     rather than a slot because they are BAM-derived working data keyed by
#'     molecule, not part of the object's SNP-by-cell counts, and so are lost
#'     by operations that rebuild the object: pass the object
#'     \code{add_molecule_phase()} returned, and subset \emph{before} that call
#'     rather than after.}
#'   \item{The SNP-to-gene map}{Read from \code{\link{snp_gene_map}}, built by
#'     \code{\link{import_cellsnp}} from the gene annotation it was given, and
#'     carried through subsetting and merging with its SNPs. It is distinct
#'     from \code{snp_info$gene_name}, which comma-joins overlapping genes into
#'     one label: attributing a molecule at a SNP overlapping two genes needs
#'     each candidate as its own row with its strand. An annotation without a
#'     \code{strand} column cannot supply that, so an object imported with one
#'     has an empty map; \code{snp_gene_map(x) <- assign_snp_genes(snp_info(x),
#'     gene_anno)} fills it in without re-importing.}
#' }
#'
#' A missing input is an error naming the step that produces it, rather than a
#' silently empty result.
#'
#' @section How the map resolves multi-gene SNPs:
#' A SNP that \code{\link{assign_snp_genes}} flagged \code{ambiguous} overlaps
#' more than one gene, and is only counted towards a given candidate's
#' \code{gene_name} for molecules whose own \code{transcript_strand} matches
#' that candidate's \code{gene_strand}; molecules with no strand information
#' (\code{NA}) or the wrong strand are dropped for that SNP rather than
#' guessed. Unambiguous rows apply regardless of strand.
#'
#' @return A tibble with one row per (\code{donor}, \code{gene_name}) pair,
#'   with columns \code{donor}, \code{gene_name}, \code{active_count},
#'   \code{inactive_count}, \code{coverage}, \code{escape_fraction}
#'   (\code{inactive_count / coverage}), \code{escapes}
#'   (\code{escape_fraction >= escape_threshold}), \code{phase_block_used}
#'   (the gene's largest phase block), \code{dominant_molecules}
#'   (molecules backing that block), \code{n_stranded_molecules}
#'   (molecules of the same gene in any other block -- counted when
#'   \code{pool_blocks} is \code{TRUE}, reported but uncounted when it is
#'   \code{FALSE}), \code{n_blocks_pooled} (blocks actually counted, always 1
#'   when \code{pool_blocks} is \code{FALSE}), and
#'   \code{discordant_block_molecules} (counted molecules sitting in a block
#'   whose own inactive fraction exceeds 0.5, i.e. one that looks oriented
#'   backwards; see \sQuote{Pooling a gene's phase blocks}).
#'
#'   With \code{by_active_x = TRUE} each row is split into two, one per active-X
#'   group, with \code{active_x} (the expressed X, "X1" or "X2") added. A group
#'   with no molecules is still reported, as coverage zero rather than being
#'   absent.
#'
#' @family X-chromosome inactivation functions
#' @export
#'
#' @examples
#' \dontrun{
#' snp_data <- assign_xci(snp_data)
#' # the molecule calls ride along on snp_data from here on
#' snp_data <- add_molecule_phase(snp_data, bam_files = c(lib1 = "lib1.bam"))
#'
#' hap <- haplotype_expression_by_molecule(snp_data)
#'
#' # Only needed if the annotation given to import_cellsnp() had no strand
#' # column, leaving snp_gene_map(snp_data) empty
#' snp_gene_map(snp_data) <- assign_snp_genes(snp_info(snp_data), gene_anno)
#' }
setGeneric(
    "haplotype_expression_by_molecule",
    function(x, escape_threshold = 0.1, by_active_x = FALSE, pool_blocks = TRUE) {
        standardGeneric("haplotype_expression_by_molecule")
    }
)

#' @rdname haplotype_expression_by_molecule
#' @include SNPData-class.R
setMethod(
    "haplotype_expression_by_molecule",
    signature(x = "SNPData"),
    function(x, escape_threshold = 0.1, by_active_x = FALSE, pool_blocks = TRUE) {
        barcode_info <- barcode_info(x)
        donor_snp_info <- donor_snp_info(x)

        if (!"active_x" %in% colnames(barcode_info) || !.has_xci_diagnostics(x)) {
            stop("No stored XCI diagnostics found. Run assign_xci(x) first.")
        }
        if (!all(c("phase_block", "allele_on_x1") %in% colnames(donor_snp_info))) {
            stop("No stored molecule phase found. Run add_molecule_phase(x) first.")
        }
        # The molecule calls are BAM-derived working data that add_molecule_phase()
        # already extracted, so they are taken from the object rather than asked
        # for: there is no other way to derive them that would agree with the
        # phase blocks stored alongside. Being an attribute, they do not survive
        # operations that rebuild the object, hence a named error rather than an
        # empty result.
        molecule_calls <- attr(x, "molecule_calls")
        if (is.null(molecule_calls)) {
            stop(
                "This object carries no molecule calls; they are attached by add_molecule_phase(x). ",
                "Attributes are lost by operations that rebuild the object, so re-run it, ",
                "or subset before it rather than after."
            )
        }
        required_call_cols <- c("donor", "barcode", "umi", "snp_id", "allele", "transcript_strand")
        missing_call_cols <- setdiff(required_call_cols, colnames(molecule_calls))
        if (length(missing_call_cols) > 0) {
            stop("molecule_calls is missing required column(s): ", paste(missing_call_cols, collapse = ", "))
        }
        # Built at import from the gene annotation, where strand is available;
        # snp_info$gene_name cannot stand in for it, being a comma-joined label
        # with no strand and no per-candidate rows.
        snp_gene_map <- snp_gene_map(x)
        if (nrow(snp_gene_map) == 0) {
            stop(
                "This object carries no SNP-to-gene map. It is built by import_cellsnp() from a gene ",
                "annotation with a strand column; set it on an existing object with ",
                "snp_gene_map(x) <- assign_snp_genes(snp_info(x), gene_anno)."
            )
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

        # Pooling is the default because every SNP reaching this point already
        # carries a globally oriented allele_on_x1, so the blocks of one gene
        # share a scale; counting only the largest is the conservative option
        # rather than the correct one (see 'Pooling a gene's phase blocks').
        counted <- calls
        if (!pool_blocks) {
            counted <- dplyr::semi_join(counted, best_block, by = c("donor", "gene_name", "phase_block"))
        }
        blocks_counted <- counted %>%
            dplyr::distinct(donor, gene_name, phase_block) %>%
            dplyr::summarise(n_blocks_pooled = dplyr::n(), .by = c(donor, gene_name))

        cell_groups <- barcode_info %>%
            dplyr::filter(active_x %in% c("X1", "X2")) %>%
            dplyr::select(barcode, donor, active_x)

        # A molecule votes across every SNP it covers in the blocks being
        # counted; majority wins, a tie is ambiguous (residual base-calling
        # noise once phase is accounted for) and dropped rather than guessed.
        # Under pooling a molecule spanning two of the gene's blocks votes with
        # all of its SNPs at once: the same rule over a wider set, not a new one.
        molecules <- counted %>%
            dplyr::summarise(n_x1 = sum(is_x1), n_x2 = sum(!is_x1), .by = c(donor, gene_name, barcode, umi)) %>%
            dplyr::mutate(
                haplotype = dplyr::case_when(n_x1 > n_x2 ~ "X1", n_x2 > n_x1 ~ "X2", TRUE ~ "ambiguous")
            ) %>%
            dplyr::filter(haplotype != "ambiguous")

        # The same vote taken per block, purely to audit orientation:
        # escape_fraction cannot exceed 0.5 by construction, so a block whose
        # own inactive fraction does is oriented backwards rather than
        # escaping. Its molecules are still counted -- flipping them here would
        # be guessing which of the two orientations is the wrong one -- but the
        # count is surfaced so the gene can be re-examined or excluded.
        discordant <- counted %>%
            dplyr::summarise(
                n_x1 = sum(is_x1),
                n_x2 = sum(!is_x1),
                .by = c(donor, gene_name, phase_block, barcode, umi)
            ) %>%
            dplyr::mutate(
                haplotype = dplyr::case_when(n_x1 > n_x2 ~ "X1", n_x2 > n_x1 ~ "X2", TRUE ~ "ambiguous")
            ) %>%
            dplyr::filter(haplotype != "ambiguous") %>%
            dplyr::inner_join(cell_groups, by = c("donor", "barcode")) %>%
            dplyr::summarise(
                block_molecules = dplyr::n(),
                block_inactive = sum(haplotype != active_x),
                .by = c(donor, gene_name, phase_block)
            ) %>%
            dplyr::summarise(
                discordant_block_molecules = sum(block_molecules[block_inactive * 2 > block_molecules]),
                .by = c(donor, gene_name)
            )

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
        by_gene <- tidyr::expand_grid(
            dplyr::distinct(counts, donor, gene_name),
            active_x = c("X1", "X2")
        ) %>%
            dplyr::left_join(counts, by = c("donor", "gene_name", "active_x")) %>%
            dplyr::mutate(
                active_count = dplyr::coalesce(active_count, 0L),
                inactive_count = dplyr::coalesce(inactive_count, 0L),
                coverage = dplyr::coalesce(coverage, 0L)
            )

        # Pooling the disjoint groups is the default here for the same reason as
        # in haplotype_expression(): one row per (donor, gene) is the grain
        # test_escape() tests at.
        if (!by_active_x) {
            by_gene <- by_gene %>%
                dplyr::summarise(
                    active_count = sum(active_count),
                    inactive_count = sum(inactive_count),
                    coverage = sum(coverage),
                    .by = c(donor, gene_name)
                )
        }

        by_gene %>%
            dplyr::mutate(
                escape_fraction = dplyr::if_else(coverage > 0, inactive_count / coverage, NA_real_),
                escapes = escape_fraction >= escape_threshold
            ) %>%
            dplyr::left_join(best_block, by = c("donor", "gene_name")) %>%
            dplyr::left_join(stranded, by = c("donor", "gene_name")) %>%
            dplyr::left_join(blocks_counted, by = c("donor", "gene_name")) %>%
            dplyr::left_join(discordant, by = c("donor", "gene_name")) %>%
            dplyr::mutate(
                n_stranded_molecules = dplyr::coalesce(n_stranded_molecules, 0L),
                discordant_block_molecules = dplyr::coalesce(discordant_block_molecules, 0L)
            ) %>%
            dplyr::rename(phase_block_used = phase_block) %>%
            dplyr::arrange(dplyr::pick(dplyr::any_of(c("donor", "gene_name", "active_x"))))
    }
)

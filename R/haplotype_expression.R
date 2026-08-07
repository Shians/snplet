#' Phased active/inactive haplotype expression per SNP
#'
#' Splits each phased SNP's ALT/REF counts into active- and
#' inactive-haplotype reads using the phase \code{\link{assign_xci}} fit
#' against its frozen cell-state calls, computed separately within the
#' X1-active and X2-active cell groups. By default this spans every
#' donor-genotyped het SNP with a stored phase, including SNPs
#' \code{assign_xci} excluded from active-X calling as uninformative or
#' escaping — exactly the SNPs worth inspecting for escape. Unlike a
#' per-group \code{pmin}/\code{pmax} on raw counts, the allele-to-haplotype
#' assignment comes from the stored per-donor phase
#' (\code{donor_snp_info$allele_on_x1}) and is applied consistently to both
#' groups, so the same physical allele can never be silently relabelled as
#' "active" in both — that contradiction is instead surfaced explicitly.
#'
#' @details
#' For each SNP and donor, the stored phase records which physical allele (REF
#' or ALT) is carried by the X1 haplotype; the other allele is carried by X2. The
#' active haplotype is the one named by \code{active_x}: X1 in the X1-active
#' group, X2 in the X2-active group; the other X is silenced. Counts are
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
#' A group whose \code{escape_fraction} exceeds 0.5 -- the "inactive"
#' haplotype outnumbering the "active" one -- is a phase contradiction
#' rather than partial leakage, and is far more likely to reflect a
#' quantification artefact (reference-mapping bias, ambient RNA, low-coverage
#' noise) than genuine escape or imprinting. Such groups are excluded from
#' the returned tibble for now, rather than reported as escaping -- with one
#' exception: \code{XIST} is transcribed almost exclusively from the
#' \emph{inactive} X by design, so both its groups legitimately read with
#' \code{escape_fraction} near 1 and are kept regardless.
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
#'
#' @return A tibble with one row per donor, phased SNP, and active-X
#'   group, with columns \code{donor} (when the object carries donor
#'   assignments), \code{snp_id}, \code{gene_name} (when the object carries gene
#'   annotation), \code{active_x} (the expressed X, "X1" or "X2"),
#'   \code{n_cells} (cells in the group), \code{active_count},
#'   \code{inactive_count}, \code{coverage}, \code{dominant_allele} (the physical
#'   allele, "REF"/"ALT", with more reads in this group; \code{NA} when
#'   uncovered), \code{escape_fraction} (\code{inactive_count / coverage}),
#'   \code{escapes} (\code{escape_fraction >= escape_threshold}), and
#'   \code{same_allele_dominant} (per SNP; \code{TRUE} when both covered groups
#'   favour the same physical allele). Groups with \code{escape_fraction > 0.5}
#'   (a phase contradiction, see Details) are dropped before this tibble is
#'   returned, except for \code{XIST}.
#'
#' @family X-chromosome inactivation functions
#' @export
#'
#' @examples
#' \dontrun{
#' snp_data <- assign_xci(snp_data)
#' hap <- haplotype_expression(snp_data)
#'
#' # SNPs that escape inactivation surface directly
#' dplyr::filter(hap, same_allele_dominant)
#'
#' # Restrict to the narrower set that drove active-X calling
#' hap_informative <- haplotype_expression(snp_data, xci_informative_only = TRUE)
#' }
setGeneric("haplotype_expression", function(x, escape_threshold = 0.1, xci_informative_only = FALSE) {
    standardGeneric("haplotype_expression")
})

#' @rdname haplotype_expression
#' @include SNPData-class.R
setMethod(
    "haplotype_expression",
    signature(x = "SNPData"),
    function(x, escape_threshold = 0.1, xci_informative_only = FALSE) {
        barcode_info <- barcode_info(x)
        snp_info <- snp_info(x)

        if (!"active_x" %in% colnames(barcode_info) || !.has_xci_diagnostics(x)) {
            stop("No stored XCI diagnostics found. Run assign_xci(x) first.")
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
            # contradicts the stored phase outright rather than showing partial
            # leakage, and is more likely a quantification artefact than genuine
            # escape -- exclude it for now rather than report it as escaping.
            # XIST is exempt: it is transcribed almost exclusively from the
            # inactive X by design, so escape_fraction near 1 in both groups is
            # expected, not a contradiction.
            dplyr::filter(
                dplyr::coalesce(gene_name == "XIST", FALSE) | is.na(escape_fraction) | escape_fraction <= 0.5
            ) %>%
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

        if (!has_donor) {
            result <- dplyr::select(result, -donor)
        }

        # Drop the placeholder gene_name column when the object carries no annotation.
        if (!"gene_name" %in% colnames(snp_info)) {
            result <- dplyr::select(result, -gene_name)
        }
        result
    }
)

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

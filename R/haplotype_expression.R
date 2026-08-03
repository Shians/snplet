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
#'   favour the same physical allele).
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
        barcode_info <- get_barcode_info(x)
        snp_info <- get_snp_info(x)

        if (!"active_x" %in% colnames(barcode_info) || !.has_xci_diagnostics(x)) {
            stop("No stored XCI diagnostics found. Run assign_xci(x) first.")
        }

        donor_snp_info <- get_donor_snp_info(x)
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

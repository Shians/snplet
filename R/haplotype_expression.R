#' Phased active/inactive haplotype expression per SNP
#'
#' Splits each informative SNP's ALT/REF counts into active- and
#' inactive-haplotype reads using the phase inferred by
#' \code{\link{assign_inactive_x}}, computed separately within the X1-inactive
#' and X2-inactive cell groups. Unlike a per-group \code{pmin}/\code{pmax} on
#' raw counts, the allele-to-haplotype assignment comes from the stored
#' per-donor phase (\code{xci_allele_on_x1_by_donor}) and is applied
#' consistently to both groups, so the same physical allele can never be
#' silently relabelled as "active" in both — that contradiction is instead
#' surfaced explicitly.
#'
#' @details
#' For each SNP and donor, the stored phase records which physical allele (REF
#' or ALT) is carried by the X1 haplotype; the other allele is carried by X2. The
#' active haplotype is the one that is \emph{not} inactivated: X2 in the
#' X1-inactive group, X1 in the X2-inactive group. Counts are assigned to
#' \code{active_count} / \code{inactive_count} accordingly.
#'
#' Under canonical XCI, the active haplotype dominates in both groups, so the
#' \emph{physical} allele that dominates should flip between the groups (X2's
#' allele in X1-inactive cells, X1's allele in X2-inactive cells). When the same
#' physical allele dominates in both groups — biologically impossible for a
#' cleanly inactivated SNP — \code{same_allele_dominant} is \code{TRUE},
#' identifying escaping or XCI-independent allelic imbalance rather than
#' flattening it. The \code{escape_fraction} column quantifies the magnitude:
#' the fraction of reads coming from the haplotype that should be silenced.
#'
#' @param x A SNPData object that had XCI diagnostics stored by
#'   \code{\link{assign_inactive_x}} or
#'   \code{\link{assign_inactive_x_by_clonotype}}.
#' @param escape_threshold Inactive-haplotype fraction at or above which a
#'   group is flagged as escaping in \code{escapes}. Default 0.1.
#'
#' @return A tibble with one row per donor, informative SNP, and inactive-X
#'   group, with columns \code{donor} (when the object carries donor
#'   assignments), \code{snp_id}, \code{gene_name} (when the object carries gene
#'   annotation), \code{inactive_x} (the silenced X, "X1" or "X2"),
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
#' snp_data <- assign_inactive_x(snp_data)
#' hap <- haplotype_expression(snp_data)
#'
#' # SNPs that escape inactivation surface directly
#' dplyr::filter(hap, same_allele_dominant)
#' }
setGeneric("haplotype_expression", function(x, escape_threshold = 0.1) {
    standardGeneric("haplotype_expression")
})

#' @rdname haplotype_expression
#' @include SNPData-class.R
setMethod("haplotype_expression", signature(x = "SNPData"), function(x, escape_threshold = 0.1) {
    barcode_info <- get_barcode_info(x)
    snp_info <- get_snp_info(x)

    if (!"inactive_x" %in% colnames(barcode_info) || !"xci_allele_on_x1_by_donor" %in% colnames(snp_info)) {
        stop("No stored XCI diagnostics found. Run assign_inactive_x(x) first.")
    }

    # Informative SNPs carry a phase; cells carry a confident inactive-X call.
    snp_keep <- !is.na(snp_info$xci_allele_on_x1_by_donor)
    cell_keep <- barcode_info$inactive_x %in% c("X1", "X2")
    if (!any(snp_keep) || !any(cell_keep)) {
        stop("No phased SNPs or confidently assigned cells available.")
    }

    ref <- ref_count(x)[snp_keep, cell_keep, drop = FALSE]
    alt <- alt_count(x)[snp_keep, cell_keep, drop = FALSE]
    group <- barcode_info$inactive_x[cell_keep]
    snp_keep_info <- snp_info[snp_keep, ]

    # The XCI fit is per-donor (X1/X2 labels and phase are donor-specific), so
    # counts must be pooled within a donor, never across donors. Absent a donor
    # column, treat the whole object as a single donor.
    has_donor <- "donor" %in% colnames(barcode_info)
    donor <- if (has_donor) as.character(barcode_info$donor[cell_keep]) else rep("all", ncol(ref))

    # Per-donor phase: a SNP retained in several donors carries a distinct
    # X1-allele in each. The *_by_donor column aligns those alleles to
    # xci_informative_donor; parse them so each donor gets its own phase rather
    # than a flattened value.
    donor_ids <- strsplit(snp_keep_info$xci_informative_donor, ",")
    donor_alleles <- strsplit(snp_keep_info$xci_allele_on_x1_by_donor, ",")
    # allele_on_x1 for one donor, per SNP (NA where the SNP was not informative
    # in that donor, so its counts cannot be phased). Without a donor column the
    # object is single-donor, so each SNP's sole stored phase applies.
    phase_for_donor <- function(d) {
        if (!has_donor) {
            return(vapply(donor_alleles, `[`, character(1), 1L))
        }
        mapply(
            function(ids, alleles) alleles[match(d, ids)],
            donor_ids,
            donor_alleles,
            USE.NAMES = FALSE
        )
    }

    # Summarise one donor x inactive-X group: pool counts across its cells, then
    # split each SNP into active- vs inactive-haplotype reads using the stored phase.
    summarise_group <- function(d, g) {
        cols <- donor == d & group == g
        ref_g <- Matrix::rowSums(ref[, cols, drop = FALSE])
        alt_g <- Matrix::rowSums(alt[, cols, drop = FALSE])

        # Map physical alleles onto haplotypes via this donor's phase.
        allele_on_x1 <- phase_for_donor(d)
        x1_count <- ifelse(allele_on_x1 == "REF", ref_g, alt_g)
        x2_count <- ifelse(allele_on_x1 == "REF", alt_g, ref_g)

        # Active haplotype is the one not silenced in this group.
        active_count <- if (g == "X1") x2_count else x1_count
        inactive_count <- if (g == "X1") x1_count else x2_count
        cov <- active_count + inactive_count

        tibble::tibble(
            donor = d,
            snp_id = snp_keep_info$snp_id,
            gene_name = snp_keep_info$gene_name %||% NA_character_,
            inactive_x = g,
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
        inactive_x = c("X1", "X2"),
        stringsAsFactors = FALSE
    )
    result <- purrr::map2(combos$donor, combos$inactive_x, summarise_group) %>%
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
})

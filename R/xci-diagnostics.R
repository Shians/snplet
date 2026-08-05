#' Extract per-cell active-X assignments
#'
#' Pulls the stored X-chromosome inactivation call and posterior for every cell
#' out of a SNPData object's barcode metadata. Requires that
#' \code{\link{assign_xci}} or \code{\link{assign_xci_by_clonotype}}
#' has been run first; errors otherwise.
#'
#' @param x A SNPData object that had XCI diagnostics stored by
#'   \code{\link{assign_xci}} or
#'   \code{\link{assign_xci_by_clonotype}}.
#'
#' @return A tibble with one row per cell and columns \code{cell_id}
#'   (character), \code{donor} (character; present only when the object carries
#'   donor information), \code{active_x} (character "X1"/"X2", or \code{NA}
#'   for cells below the confidence threshold), and \code{xci_post_X1_active}
#'   (numeric posterior probability that X1 is the active chromosome).
#'
#' @family X-chromosome inactivation functions
#' @export
setGeneric("xci_assignments", function(x) standardGeneric("xci_assignments"))

#' @rdname xci_assignments
#' @include SNPData-class.R
setMethod("xci_assignments", signature(x = "SNPData"), function(x) {
    barcode_info <- barcode_info(x)
    if (!"active_x" %in% colnames(barcode_info)) {
        stop("No stored XCI diagnostics found. Run assign_xci(x) first.")
    }
    barcode_info %>%
        dplyr::select(cell_id, dplyr::any_of("donor"), active_x, xci_post_X1_active)
})

#' Extract SNP haplotypes from a SNPData object with stored XCI diagnostics
#'
#' Returns the inferred phase and escape fraction for each SNP in the final
#' gene set used by the EM model. Phase is donor-specific — a SNP retained in
#' several donors carries a distinct X1-allele in each — so the result has one
#' row per SNP and donor rather than flattening to a single value.
#'
#' @param x A SNPData object that had XCI diagnostics stored by
#'   \code{\link{assign_xci}} or
#'   \code{\link{assign_xci_by_clonotype}}.
#'
#' @return A tibble with one row per informative SNP and donor, with columns
#'   \code{snp_id} (character), \code{gene_name} (character; present only when
#'   the object carries gene annotation), \code{donor} (character),
#'   \code{allele_on_x1} (character, "REF" or "ALT" — the allele carried by the
#'   X1 haplotype in that donor's model), and \code{escape_fraction} (numeric
#'   estimated fraction of reads from the inactive allele in that donor).
#'
#' @family X-chromosome inactivation functions
#' @importFrom tidyr separate_longer_delim
#' @export
setGeneric("xci_haplotypes", function(x) standardGeneric("xci_haplotypes"))

#' @rdname xci_haplotypes
#' @include SNPData-class.R
setMethod("xci_haplotypes", signature(x = "SNPData"), function(x) {
    if (!.has_xci_diagnostics(x)) {
        stop("No stored XCI diagnostics found. Run assign_xci(x) first.")
    }
    snp_info <- snp_info(x)
    donor_snp_info(x) %>%
        dplyr::filter(xci_informative) %>%
        dplyr::left_join(dplyr::select(snp_info, snp_id, dplyr::any_of("gene_name")), by = "snp_id") %>%
        dplyr::select(
            snp_id,
            dplyr::any_of("gene_name"),
            donor,
            allele_on_x1,
            escape_fraction = xci_escape_fraction
        )
})

#' Write XCI diagnostics from a fit into a SNPData object
#'
#' Promotes the fit diagnostics into the object's indexable metadata slots so
#' they survive subsetting: barcode metadata gains \code{active_x} and
#' \code{xci_post_X1_active}; \code{donor_snp_info} gains one row per phased
#' SNP x donor with \code{xci_informative} (whether the gene drove active-X
#' calling), \code{allele_on_x1} and \code{xci_escape_fraction}. For a
#' clonotype-level fit the per-cell projection is used for barcode
#' annotation.
#'
#' @keywords internal
.store_xci_fit <- function(x, fit) {
    donor_fits <- purrr::keep(fit, ~ !is.null(.x))

    barcode_diag <- purrr::map(donor_fits, function(f) {
        assignments <- if (!is.null(f$cell_assignments)) f$cell_assignments else f$assignments
        assignments %>%
            dplyr::transmute(
                cell_id,
                # assignment already names the active X; "unassigned" means no call.
                active_x = dplyr::na_if(assignment, "unassigned"),
                xci_post_X1_active = post_X1_active,
                # Record the unit the model was fit on so downstream plotting can
                # aggregate cells back to clonotypes when appropriate.
                xci_fit_unit = f$unit %||% "cell"
            )
    }) %>%
        dplyr::bind_rows()

    snp_diag <- purrr::map(donor_fits, function(f) {
        f$haplotypes %>%
            dplyr::transmute(
                snp_id,
                donor = f$donor,
                xci_informative,
                allele_on_x1,
                xci_escape_fraction = escape_fraction
            )
    }) %>%
        dplyr::bind_rows()

    if (nrow(barcode_diag) > 0) {
        x <- add_barcode_metadata(x, barcode_diag, join_by = "cell_id", overwrite = TRUE)
    }
    if (nrow(snp_diag) > 0) {
        x <- add_donor_snp_metadata(x, snp_diag, join_by = c("snp_id", "donor"), overwrite = TRUE)
    }
    x
}

#' Whether a SNPData object carries stored XCI fit diagnostics
#'
#' \code{\link{.store_xci_fit}} is the only writer of \code{xci_informative},
#' so at least one \code{TRUE} in \code{donor_snp_info} is a reliable signal
#' that \code{assign_xci}/\code{assign_xci_by_clonotype} has run, even though
#' the column itself may already exist (as \code{NA}) from zygosity calling.
#'
#' @keywords internal
.has_xci_diagnostics <- function(x) {
    donor_snp_info <- donor_snp_info(x)
    "xci_informative" %in% colnames(donor_snp_info) && any(donor_snp_info$xci_informative, na.rm = TRUE)
}

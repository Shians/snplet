#' Assign inactive X chromosome to cells
#'
#' Identifies which X chromosome is inactive in female cells based on allelic
#' imbalance at heterozygous SNPs using an Expectation-Maximisation (EM)
#' algorithm with a beta-binomial likelihood.
#'
#' @details
#' X-chromosome inactivation (XCI) is a dosage compensation mechanism in female
#' mammals where one of the two X chromosomes is randomly silenced in each cell.
#' This function infers which X is inactive by fitting an EM model to allelic
#' read counts at heterozygous SNPs on the X chromosome.
#'
#' The algorithm works as follows:
#' \enumerate{
#'   \item Filters to heterozygous SNPs on the X chromosome for each donor
#'   \item Selects the SNP with highest coverage per gene to avoid redundancy
#'   \item Removes outlier genes with atypical allelic skew
#'   \item Runs a beta-binomial EM algorithm with multiple random initialisations
#'   \item Assigns cells to X1 or X2 based on posterior probability
#' }
#'
#' Cells that do not meet \code{confidence_threshold} receive \code{NA} in
#' the \code{inactive_x} column.
#'
#' @param x SNPData object containing X chromosome SNP data with donor
#'   assignments and heterozygosity information
#' @param n_inits Number of random initialisations for the EM algorithm.
#'   The run with the highest log-likelihood is returned. Default 10.
#' @param confidence_threshold Posterior probability threshold for hard
#'   assignment. Cells with \code{post_X1 >= confidence_threshold} are
#'   assigned "X1", those with \code{post_X1 <= 1 - confidence_threshold}
#'   are assigned "X2", and the remainder receive \code{NA}. Default 0.95.
#' @param refit_after_filter Logical; if TRUE, re-run the EM algorithm after
#'   filtering genes with inconsistent allelic patterns. Provides sharper
#'   posteriors on the cleaned gene set. Default FALSE.
#'
#' @return SNPData object with an additional \code{inactive_x} column in
#'   barcode metadata, with values "X1" or "X2" indicating the inferred
#'   inactive X chromosome state. Cells that do not meet the confidence
#'   threshold receive \code{NA}.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Assign inactive X chromosome to cells
#' snp_data <- assign_inactive_x(snp_data)
#'
#' # View results
#' get_barcode_info(snp_data) %>%
#'   count(donor, inactive_x)
#' }
setGeneric("assign_inactive_x", function(x, n_inits = 10, confidence_threshold = 0.95, refit_after_filter = FALSE) {
    standardGeneric("assign_inactive_x")
})

#' @rdname assign_inactive_x
#' @include SNPData-class.R
setMethod(
    "assign_inactive_x",
    signature(x = "SNPData"),
    function(x, n_inits = 10, confidence_threshold = 0.95, refit_after_filter = FALSE) {
        unique_donors <- sort(unique(get_barcode_info(x)$donor))
        result <- furrr::future_map(
            unique_donors,
            function(d) {
                tryCatch(
                    .assign_inactive_x_single_donor(
                        filter_samples(x, donor == d),
                        n_inits,
                        confidence_threshold,
                        refit_after_filter
                    ),
                    error = function(e) {
                        logger::log_warn("Failed to assign inactive X for donor {d}: {conditionMessage(e)}")
                        tibble::tibble(cell_id = character(), inactive_x = character())
                    }
                )
            },
            .options = furrr::furrr_options(packages = "snplet", seed = TRUE)
        ) %>%
            dplyr::bind_rows()
        add_barcode_metadata(x, result)
    }
)

#' Assign inactive X chromosome to cells by clonotype
#'
#' Identifies which X chromosome is inactive in female cells based on allelic
#' imbalance at heterozygous SNPs, aggregating counts to the clonotype level
#' before running the EM model and projecting assignments back to individual cells.
#'
#' @details
#' This function is similar to \code{\link{assign_inactive_x}} but runs the
#' beta-binomial EM on clonotype-aggregated counts rather than per-cell counts.
#' This approach:
#' \itemize{
#'   \item Increases statistical power by aggregating counts across clonally related cells
#'   \item Reduces noise from sparse cell-level data
#'   \item Is biologically motivated since clonally related cells should share X-inactivation state
#' }
#'
#' The algorithm works as follows:
#' \enumerate{
#'   \item Filters to heterozygous SNPs on the X chromosome for each donor
#'   \item Selects the SNP with highest coverage per gene to avoid redundancy
#'   \item Aggregates ALT and REF counts by clonotype
#'   \item Removes outlier genes with atypical allelic skew
#'   \item Runs a beta-binomial EM algorithm with multiple random initialisations
#'   \item Assigns clonotypes to X1 or X2 based on posterior probability
#'   \item Projects clonotype assignments back to individual cells
#' }
#'
#' @param x SNPData object containing X chromosome SNP data with donor
#'   assignments, clonotype information, and heterozygosity information
#' @param n_inits Number of random initialisations for the EM algorithm.
#'   The run with the highest log-likelihood is returned. Default 10.
#' @param confidence_threshold Posterior probability threshold for hard
#'   assignment. Clonotypes with \code{post_X1 >= confidence_threshold} are
#'   assigned "X1", those with \code{post_X1 <= 1 - confidence_threshold}
#'   are assigned "X2", and the remainder receive \code{NA}. Default 0.95.
#' @param refit_after_filter Logical; if TRUE, re-run the EM algorithm after
#'   filtering genes with inconsistent allelic patterns. Provides sharper
#'   posteriors on the cleaned gene set. Default FALSE.
#'
#' @return SNPData object with an additional \code{inactive_x} column in
#'   barcode metadata, with values "X1" or "X2" indicating the inferred
#'   inactive X chromosome state. Cells from clonotypes that do not meet the
#'   confidence threshold receive \code{NA}.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Assign inactive X chromosome to cells by clonotype
#' snp_data <- assign_inactive_x_by_clonotype(snp_data)
#'
#' # View results
#' get_barcode_info(snp_data) %>%
#'   count(donor, clonotype, inactive_x)
#' }
setGeneric(
    "assign_inactive_x_by_clonotype",
    function(x, n_inits = 10, confidence_threshold = 0.95, refit_after_filter = FALSE) {
        standardGeneric("assign_inactive_x_by_clonotype")
    }
)

#' @rdname assign_inactive_x_by_clonotype
#' @include SNPData-class.R
setMethod(
    "assign_inactive_x_by_clonotype",
    signature(x = "SNPData"),
    function(x, n_inits = 10, confidence_threshold = 0.95, refit_after_filter = FALSE) {
        .check_clonotype_available(x)
        unique_donors <- sort(unique(get_barcode_info(x)$donor))
        result <- furrr::future_map(
            unique_donors,
            function(d) {
                tryCatch(
                    .assign_inactive_x_single_donor_by_clonotype(
                        filter_samples(x, donor == d),
                        n_inits,
                        confidence_threshold,
                        refit_after_filter
                    ),
                    error = function(e) {
                        logger::log_warn("Failed to assign inactive X for donor {d}: {conditionMessage(e)}")
                        tibble::tibble(cell_id = character(), inactive_x = character())
                    }
                )
            },
            .options = furrr::furrr_options(packages = "snplet", seed = TRUE)
        ) %>%
            dplyr::bind_rows()
        add_barcode_metadata(x, result)
    }
)

#' Fit X-chromosome inactivation model
#'
#' Runs the beta-binomial EM algorithm for each donor and returns a full fit
#' object containing per-unit assignments, per-gene haplotypes, and the final
#' allele count matrices used in the model. The modelling unit is either the
#' cell (\code{by = "cell"}) or the clonotype (\code{by = "clonotype"}). Use
#' this when diagnostic outputs are needed; for barcode annotation only, use
#' \code{\link{assign_inactive_x}}.
#'
#' @param x SNPData object containing X chromosome SNP data with donor
#'   assignments and heterozygosity information. For \code{by = "clonotype"},
#'   clonotype information must also be present.
#' @param n_inits Number of random initialisations for the EM algorithm.
#'   Default 10.
#' @param confidence_threshold Posterior probability threshold for hard
#'   assignment. Default 0.95.
#' @param refit_after_filter Logical; if TRUE, re-run the EM algorithm after
#'   filtering genes with inconsistent allelic patterns. Provides sharper
#'   posteriors on the cleaned gene set. Default FALSE.
#' @param by Modelling unit. \code{"cell"} runs the EM on per-cell counts;
#'   \code{"clonotype"} aggregates ALT/REF counts by clonotype first and
#'   projects assignments back to cells. Default \code{"cell"}.
#' @param store Logical; if TRUE, return the input SNPData object with the
#'   diagnostics written into its metadata slots rather than returning an
#'   \code{xci_fit} object. Barcode metadata gains \code{inactive_x} and
#'   \code{xci_post_X1}; SNP metadata gains \code{xci_informative},
#'   \code{xci_allele_on_x1}, and \code{xci_escape_fraction}. These columns survive
#'   subsetting because they live in the object's indexable slots. Default
#'   FALSE.
#'
#' @return If \code{store = FALSE}, an \code{xci_fit} object (named list by
#'   donor); use \code{\link{xci_assignments}}, \code{\link{xci_haplotypes}},
#'   and \code{\link{plot_inactive_x_assignment_heatmap}} to extract results.
#'   If \code{store = TRUE}, the input SNPData object with diagnostics written
#'   into its metadata slots; the same accessors and heatmap function accept
#'   the returned SNPData object.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' fit <- fit_inactive_x(snp_data)
#' xci_assignments(fit)
#' xci_haplotypes(fit)
#' plot_inactive_x_assignment_heatmap(fit, donor = "donor1")
#'
#' # Clonotype-level fit, diagnostics stored back into the object
#' snp_data <- fit_inactive_x(snp_data, by = "clonotype", store = TRUE)
#' xci_assignments(snp_data)
#' plot_inactive_x_assignment_heatmap(snp_data, donor = "donor1")
#' }
setGeneric(
    "fit_inactive_x",
    function(
        x,
        n_inits = 10,
        confidence_threshold = 0.95,
        refit_after_filter = FALSE,
        by = c("cell", "clonotype"),
        store = FALSE
    ) {
        standardGeneric("fit_inactive_x")
    }
)

#' @rdname fit_inactive_x
#' @include SNPData-class.R
setMethod(
    "fit_inactive_x",
    signature(x = "SNPData"),
    function(
        x,
        n_inits = 10,
        confidence_threshold = 0.95,
        refit_after_filter = FALSE,
        by = c("cell", "clonotype"),
        store = FALSE
    ) {
        by <- match.arg(by)
        if (by == "clonotype") {
            .check_clonotype_available(x)
        }

        donor_fitter <- if (by == "clonotype") .fit_xci_donor_by_clonotype else .fit_xci_donor
        unique_donors <- sort(unique(get_barcode_info(x)$donor))

        result <- furrr::future_map(
            unique_donors,
            function(d) {
                tryCatch(
                    donor_fitter(filter_samples(x, donor == d), n_inits, confidence_threshold, refit_after_filter),
                    error = function(e) {
                        logger::log_warn("Failed to fit XCI for donor {d}: {conditionMessage(e)}")
                        NULL
                    }
                )
            },
            .options = furrr::furrr_options(packages = "snplet", seed = TRUE)
        ) %>%
            magrittr::set_names(unique_donors)

        fit <- structure(result, class = "xci_fit")

        if (store) {
            return(.store_xci_fit(x, fit))
        }
        fit
    }
)

#' @keywords internal
.check_clonotype_available <- function(x) {
    barcode_info <- get_barcode_info(x)
    if (!"clonotype" %in% colnames(barcode_info)) {
        stop(
            "Clonotype information not available. Add clonotype data using add_barcode_metadata() or import_cellsnp() with vdj_file parameter."
        )
    }
    if (all(is.na(barcode_info$clonotype))) {
        stop(
            "All clonotype values are NA. Cannot perform clonotype-level X-inactivation assignment. Add clonotype data using add_barcode_metadata() or import_cellsnp() with vdj_file parameter."
        )
    }
    invisible(TRUE)
}

# Register the S3 `xci_fit` class so it can be used in S4 method signatures.
methods::setOldClass("xci_fit")

#' Extract assignments from an XCI fit or a stored SNPData object
#'
#' @param x An \code{xci_fit} object from \code{\link{fit_inactive_x}}, or a
#'   SNPData object that had diagnostics stored via
#'   \code{fit_inactive_x(..., store = TRUE)}.
#'
#' @return A tibble of assignments. For an \code{xci_fit} the columns are
#'   \code{donor}, \code{cell_id}, \code{post_X1}, \code{post_X2}, and
#'   \code{assignment} (for a clonotype-level fit, \code{cell_id} holds the
#'   clonotype ID). For a SNPData object the columns are \code{cell_id},
#'   \code{donor}, \code{inactive_x}, and \code{xci_post_X1}.
#'
#' @export
setGeneric("xci_assignments", function(x) standardGeneric("xci_assignments"))

#' @rdname xci_assignments
setMethod("xci_assignments", signature(x = "xci_fit"), function(x) {
    purrr::keep(x, ~ !is.null(.x)) %>%
        purrr::map(~ .x$assignments) %>%
        dplyr::bind_rows(.id = "donor")
})

#' @rdname xci_assignments
#' @include SNPData-class.R
setMethod("xci_assignments", signature(x = "SNPData"), function(x) {
    barcode_info <- get_barcode_info(x)
    if (!"inactive_x" %in% colnames(barcode_info)) {
        stop("No stored XCI diagnostics found. Run fit_inactive_x(x, store = TRUE) first.")
    }
    barcode_info %>%
        dplyr::select(cell_id, dplyr::any_of("donor"), inactive_x, xci_post_X1)
})

#' Extract SNP haplotypes from an XCI fit or a stored SNPData object
#'
#' Returns the inferred phase and escape fraction for each SNP in the final
#' gene set used by the EM model.
#'
#' @param x An \code{xci_fit} object from \code{\link{fit_inactive_x}}, or a
#'   SNPData object that had diagnostics stored via
#'   \code{fit_inactive_x(..., store = TRUE)}.
#'
#' @return A tibble with the inferred phase (\code{allele_on_x1}, "REF" or
#'   "ALT") and estimated escape fraction (\code{escape_fraction}) for each informative
#'   SNP. An \code{xci_fit} additionally carries \code{donor} and
#'   \code{gene_name}.
#'
#' @export
setGeneric("xci_haplotypes", function(x) standardGeneric("xci_haplotypes"))

#' @rdname xci_haplotypes
setMethod("xci_haplotypes", signature(x = "xci_fit"), function(x) {
    purrr::keep(x, ~ !is.null(.x)) %>%
        purrr::map(~ .x$haplotypes) %>%
        dplyr::bind_rows(.id = "donor")
})

#' @rdname xci_haplotypes
#' @include SNPData-class.R
setMethod("xci_haplotypes", signature(x = "SNPData"), function(x) {
    snp_info <- get_snp_info(x)
    if (!"xci_informative" %in% colnames(snp_info)) {
        stop("No stored XCI diagnostics found. Run fit_inactive_x(x, store = TRUE) first.")
    }
    snp_info %>%
        dplyr::filter(xci_informative) %>%
        dplyr::select(
            snp_id,
            dplyr::any_of("gene_name"),
            allele_on_x1 = xci_allele_on_x1,
            escape_fraction = xci_escape_fraction
        )
})

#' Plot assignment heatmap from an XCI fit or a stored SNPData object
#'
#' Visualizes the REF allele fraction at the final EM-selected genes for a
#' single donor. Columns (cells or clonotypes) are ordered by posterior and
#' annotated at the top. Rows are the genes that survived both the outlier
#' filter and the post-convergence escapee filter.
#'
#' @param x An \code{xci_fit} object from \code{\link{fit_inactive_x}}, or a
#'   SNPData object that had diagnostics stored via
#'   \code{fit_inactive_x(..., store = TRUE)}.
#' @param donor Character string specifying which donor to visualize
#'
#' @return A ComplexHeatmap object
#'
#' @importFrom circlize colorRamp2
#' @export
setGeneric("plot_inactive_x_assignment_heatmap", function(x, donor) {
    standardGeneric("plot_inactive_x_assignment_heatmap")
})

#' @rdname plot_inactive_x_assignment_heatmap
setMethod("plot_inactive_x_assignment_heatmap", signature(x = "xci_fit"), function(x, donor) {
    donor_fit <- x[[donor]]
    if (is.null(donor_fit)) {
        stop(paste0("No fit available for donor: ", donor))
    }
    gene_name_map <- stats::setNames(donor_fit$haplotypes$gene_name, donor_fit$haplotypes$snp_id)
    .plot_xci_heatmap_from_parts(
        ref_mat = donor_fit$ref_mat,
        alt_mat = donor_fit$alt_mat,
        assignment = donor_fit$assignments$assignment,
        post_X1 = donor_fit$assignments$post_X1,
        unit_ids = donor_fit$assignments$cell_id,
        gene_name_map = gene_name_map
    )
})

#' @rdname plot_inactive_x_assignment_heatmap
#' @include SNPData-class.R
setMethod("plot_inactive_x_assignment_heatmap", signature(x = "SNPData"), function(x, donor) {
    barcode_info <- get_barcode_info(x)
    snp_info <- get_snp_info(x)
    if (!"inactive_x" %in% colnames(barcode_info) || !"xci_informative" %in% colnames(snp_info)) {
        stop("No stored XCI diagnostics found. Run fit_inactive_x(x, store = TRUE) first.")
    }

    donor_data <- filter_samples(x, donor == !!donor)
    donor_data <- filter_snps(donor_data, xci_informative)

    donor_snp_info <- get_snp_info(donor_data)
    donor_barcode_info <- get_barcode_info(donor_data)

    # Cells this donor's model actually assigned (drop NA inactive_x)
    assigned <- donor_barcode_info %>%
        dplyr::mutate(assignment = ifelse(is.na(inactive_x), "unassigned", inactive_x))

    ref_mat <- ref_count(donor_data)
    alt_mat <- alt_count(donor_data)
    rownames(ref_mat) <- donor_snp_info$snp_id
    rownames(alt_mat) <- donor_snp_info$snp_id

    gene_name_map <- stats::setNames(donor_snp_info$gene_name, donor_snp_info$snp_id)
    .plot_xci_heatmap_from_parts(
        ref_mat = ref_mat,
        alt_mat = alt_mat,
        assignment = assigned$assignment,
        post_X1 = assigned$xci_post_X1,
        unit_ids = assigned$cell_id,
        gene_name_map = gene_name_map
    )
})

#' @keywords internal
.plot_xci_heatmap_from_parts <- function(ref_mat, alt_mat, assignment, post_X1, unit_ids, gene_name_map) {
    # REF allele fraction; NA where uncovered
    cov_mat <- ref_mat + alt_mat
    frac_mat <- ref_mat / cov_mat
    frac_mat[cov_mat == 0] <- NA

    # Replace SNP ID row names with gene names
    rownames(frac_mat) <- gene_name_map[rownames(frac_mat)]

    # Drop rows that are entirely NA — they carry no signal and break hclust
    row_has_data <- rowSums(!is.na(frac_mat)) > 0
    frac_mat <- frac_mat[row_has_data, , drop = FALSE]

    # Order units: X1 → X2 → unassigned
    idx <- order(post_X1, decreasing = TRUE)
    frac_mat <- frac_mat[, unit_ids[idx], drop = FALSE]

    col_ann <- ComplexHeatmap::HeatmapAnnotation(
        assignment = assignment[idx],
        posterior_X1 = post_X1[idx],
        col = list(
            assignment = c(X1 = "#4e79a7", X2 = "#f28e2b", unassigned = "grey70"),
            posterior_X1 = circlize::colorRamp2(c(0, 0.1, 0.9, 1), c("#e9a3c9", "white", "white", "#a1d76a"))
        ),
        show_annotation_name = TRUE
    )

    # Impute NAs with 0.5 (balanced) only for row clustering so that sparse
    # coverage does not produce NaN distances; the display matrix keeps NAs.
    frac_for_cluster <- frac_mat
    frac_for_cluster[is.na(frac_for_cluster)] <- 0.5
    row_dend <- stats::hclust(stats::dist(as.matrix(frac_for_cluster)), method = "ward.D2")

    ComplexHeatmap::Heatmap(
        as.matrix(frac_mat),
        top_annotation = col_ann,
        cluster_columns = FALSE,
        cluster_rows = row_dend,
        show_row_names = TRUE,
        show_column_names = FALSE,
        name = "REF fraction",
        na_col = "grey95"
    )
}

#' @keywords internal
.fit_xci_donor <- function(snp_data, n_inits = 10, confidence_threshold = 0.95, refit_after_filter = FALSE) {
    donor <- unique(get_barcode_info(snp_data)$donor)
    logger::log_info("[{donor}] Fitting XCI model")

    snp_data <- .filter_to_informative_het_snps(snp_data, donor)

    ref_mat <- ref_count(snp_data)
    alt_mat <- alt_count(snp_data)
    cell_ids <- colnames(ref_mat)
    snp_info <- get_snp_info(snp_data)

    if (nrow(ref_mat) == 0 || ncol(ref_mat) == 0) {
        logger::log_warn("Insufficient data for donor {donor}")
        return(NULL)
    }

    xci_result <- .infer_xci(
        ref_mat,
        alt_mat,
        n_inits = n_inits,
        confidence_threshold = confidence_threshold,
        refit_after_filter = refit_after_filter,
        donor = donor
    )

    .assemble_xci_fit(
        xci_result = xci_result,
        donor = donor,
        ref_mat = ref_mat,
        alt_mat = alt_mat,
        snp_info = snp_info,
        unit_ids = cell_ids,
        unit_label = "cells"
    )
}

#' @keywords internal
.fit_xci_donor_by_clonotype <- function(
    snp_data,
    n_inits = 10,
    confidence_threshold = 0.95,
    refit_after_filter = FALSE
) {
    donor <- unique(get_barcode_info(snp_data)$donor)
    logger::log_info("[{donor}] Fitting XCI model at clonotype level")

    snp_data <- .filter_to_informative_het_snps(snp_data, donor)

    barcode_info <- get_barcode_info(snp_data)
    has_clonotype <- !is.na(barcode_info$clonotype)
    if (!any(has_clonotype)) {
        stop(glue::glue("No cells with non-NA clonotype values for donor {donor}"))
    }
    snp_data <- snp_data[, has_clonotype]

    cell_to_clonotype <- get_barcode_info(snp_data) %>%
        dplyr::select(cell_id, clonotype)
    clonotypes <- cell_to_clonotype$clonotype

    ref_mat <- groupedRowSums(ref_count(snp_data), clonotypes)
    alt_mat <- groupedRowSums(alt_count(snp_data), clonotypes)
    clonotype_ids <- colnames(ref_mat)
    snp_info <- get_snp_info(snp_data)

    if (nrow(ref_mat) == 0 || ncol(ref_mat) == 0) {
        logger::log_warn("Insufficient data for donor {donor}")
        return(NULL)
    }

    xci_result <- .infer_xci(
        ref_mat,
        alt_mat,
        n_inits = n_inits,
        confidence_threshold = confidence_threshold,
        refit_after_filter = refit_after_filter,
        donor = donor
    )

    fit <- .assemble_xci_fit(
        xci_result = xci_result,
        donor = donor,
        ref_mat = ref_mat,
        alt_mat = alt_mat,
        snp_info = snp_info,
        unit_ids = clonotype_ids,
        unit_label = "clonotypes"
    )

    # Project clonotype-level posteriors down to individual cells so downstream
    # cell annotation works. Each cell inherits its clonotype's posterior.
    fit$cell_assignments <- fit$assignments %>%
        dplyr::rename(clonotype = cell_id) %>%
        dplyr::inner_join(cell_to_clonotype, by = "clonotype") %>%
        dplyr::select(cell_id, post_X1, post_X2, assignment)

    fit$unit <- "clonotype"
    fit
}

#' Assemble a per-donor XCI fit from an EM result
#'
#' Shared post-EM assembly used by both the cell- and clonotype-level donor
#' fitters. \code{unit_ids} label the columns of the count matrices (cell IDs
#' or clonotype IDs); the resulting \code{assignments} tibble keys on
#' \code{cell_id} regardless so accessors and the heatmap treat the modelling
#' unit uniformly.
#'
#' @keywords internal
.assemble_xci_fit <- function(xci_result, donor, ref_mat, alt_mat, snp_info, unit_ids, unit_label) {
    assignments <- xci_result$post %>%
        dplyr::mutate(cell_id = unit_ids[cell], post_X2 = 1 - post_X1) %>%
        dplyr::select(cell_id, post_X1, post_X2, assignment)

    snp_info_filtered <- snp_info[xci_result$gene_keep, ]
    haplotypes <- tibble::tibble(
        snp_id = snp_info_filtered$snp_id,
        gene_name = snp_info_filtered$gene_name,
        allele_on_x1 = ifelse(xci_result$h_g == 0, "REF", "ALT"),
        escape_fraction = xci_result$pi_g
    )

    counts <- table(factor(assignments$assignment, c("X1", "X2", "unassigned")))
    logger::log_info(
        "[{donor}] XCI fit complete: {nrow(haplotypes)} genes retained, {nrow(assignments)} {unit_label} ",
        "(X1={counts[['X1']]}, X2={counts[['X2']]}, unassigned={counts[['unassigned']]})"
    )

    ref_mat_filtered <- ref_mat[xci_result$gene_keep, , drop = FALSE]
    alt_mat_filtered <- alt_mat[xci_result$gene_keep, , drop = FALSE]
    rownames(ref_mat_filtered) <- snp_info_filtered$snp_id
    rownames(alt_mat_filtered) <- snp_info_filtered$snp_id
    colnames(ref_mat_filtered) <- unit_ids
    colnames(alt_mat_filtered) <- unit_ids

    list(
        donor = donor,
        unit = "cell",
        assignments = assignments,
        haplotypes = haplotypes,
        ref_mat = ref_mat_filtered,
        alt_mat = alt_mat_filtered
    )
}

#' Write XCI diagnostics from a fit into a SNPData object
#'
#' Promotes the fit diagnostics into the object's indexable metadata slots so
#' they survive subsetting: barcode metadata gains \code{inactive_x} and
#' \code{xci_post_X1}; SNP metadata gains \code{xci_informative},
#' \code{xci_allele_on_x1}, and \code{xci_escape_fraction}. For a clonotype-level fit the
#' per-cell projection is used for barcode annotation.
#'
#' @keywords internal
.store_xci_fit <- function(x, fit) {
    donor_fits <- purrr::keep(fit, ~ !is.null(.x))

    barcode_diag <- purrr::map(donor_fits, function(f) {
        assignments <- if (!is.null(f$cell_assignments)) f$cell_assignments else f$assignments
        assignments %>%
            dplyr::transmute(
                cell_id,
                inactive_x = ifelse(assignment == "unassigned", NA_character_, assignment),
                xci_post_X1 = post_X1
            )
    }) %>%
        dplyr::bind_rows()

    snp_diag <- purrr::map(donor_fits, function(f) {
        f$haplotypes %>%
            dplyr::transmute(
                snp_id,
                xci_informative = TRUE,
                xci_allele_on_x1 = allele_on_x1,
                xci_escape_fraction = escape_fraction
            )
    }) %>%
        dplyr::bind_rows()

    if (nrow(barcode_diag) > 0) {
        x <- add_barcode_metadata(x, barcode_diag, join_by = "cell_id", overwrite = TRUE)
    }
    if (nrow(snp_diag) > 0) {
        x <- add_snp_metadata(x, snp_diag, join_by = "snp_id", overwrite = TRUE)
        # SNPs not in the final EM set are non-informative rather than NA
        snp_info <- get_snp_info(x)
        snp_info$xci_informative[is.na(snp_info$xci_informative)] <- FALSE
        x <- add_snp_metadata(
            x,
            dplyr::select(snp_info, snp_id, xci_informative),
            join_by = "snp_id",
            overwrite = TRUE
        )
    }
    x
}

#' @keywords internal
.assign_inactive_x_single_donor <- function(
    snp_data,
    n_inits = 10,
    confidence_threshold = 0.95,
    refit_after_filter = FALSE
) {
    fit <- .fit_xci_donor(snp_data, n_inits, confidence_threshold, refit_after_filter)
    if (is.null(fit)) {
        return(tibble::tibble(cell_id = character(), inactive_x = character()))
    }
    fit$assignments %>%
        dplyr::filter(assignment != "unassigned") %>%
        dplyr::select(cell_id, inactive_x = assignment)
}

#' @keywords internal
.assign_inactive_x_single_donor_by_clonotype <- function(
    snp_data,
    n_inits = 10,
    confidence_threshold = 0.95,
    refit_after_filter = FALSE
) {
    fit <- .fit_xci_donor_by_clonotype(snp_data, n_inits, confidence_threshold, refit_after_filter)
    if (is.null(fit)) {
        return(tibble::tibble(cell_id = character(), inactive_x = character()))
    }
    fit$cell_assignments %>%
        dplyr::filter(assignment != "unassigned") %>%
        dplyr::select(cell_id, inactive_x = assignment)
}

#' @keywords internal
.filter_to_informative_het_snps <- function(snp_data, donor = NULL) {
    # Suppress the generic per-filter logs from filter_snps; we emit a single
    # consolidated, donor-labelled summary instead.
    old_threshold <- logger::log_threshold()
    logger::log_threshold(logger::WARN)
    on.exit(logger::log_threshold(old_threshold), add = TRUE)

    n_start <- nrow(get_snp_info(snp_data))

    het_snp_ids <- snp_data %>%
        donor_het_status_df() %>%
        dplyr::filter(zygosity == "het") %>%
        dplyr::pull(snp_id)

    snp_data <- snp_data %>%
        filter_snps(snp_id %in% het_snp_ids)
    n_het <- nrow(get_snp_info(snp_data))

    top_snp_per_gene <- get_snp_info(snp_data) %>%
        dplyr::arrange(dplyr::desc(coverage)) %>%
        dplyr::slice_head(n = 1, by = "gene_name")

    snp_data <- snp_data %>%
        filter_snps(snp_id %in% top_snp_per_gene$snp_id)
    n_genes <- nrow(get_snp_info(snp_data))

    logger::log_info(
        "[{donor}] het-SNP filter: {n_start} X SNPs -> {n_het} het -> {n_genes} genes (top SNP per gene)"
    )

    snp_data
}

#' @keywords internal
.infer_xci <- function(
    ref_mat,
    alt_mat,
    n_inits = 10,
    confidence_threshold = 0.95,
    min_cells = 10,
    min_cov = 1,
    refit_after_filter = FALSE,
    donor = NULL
) {
    passes_outlier_filter <- .filter_outlier_genes(ref_mat, alt_mat, min_cells, min_cov)
    ref_mat <- ref_mat[passes_outlier_filter, , drop = FALSE]
    alt_mat <- alt_mat[passes_outlier_filter, , drop = FALSE]

    dat <- .pivot_counts_to_long(ref_mat, alt_mat, min_cov)
    n_genes <- nrow(ref_mat)

    logger::log_info(
        "[{donor}] Running EM: {n_inits} random restarts over {n_genes} genes, {nrow(dat)} observations"
    )
    fits <- lapply(seq_len(n_inits), function(s) {
        fit <- .run_em(dat, n_genes, init_seed = s)
        logger::log_debug("[{donor}] EM restart {s}/{n_inits} done (logLik = {round(fit$ll, 2)})")
        fit
    })
    best <- fits[[which.max(sapply(fits, `[[`, "ll"))]]
    logger::log_info("[{donor}] EM complete: best logLik = {round(best$ll, 2)}")

    # Post-convergence escapee filter: genes with LLR <= 0 are inconsistent with
    # current cell assignments; MAD filter on pi_g removes outlier escape fractions.
    passes_escapee_filter <- .filter_escapee_genes(dat, n_genes, best$h_g, best$pi_g, best$rho, best$post)
    escaped_gene_indices <- which(passes_escapee_filter)

    if (refit_after_filter && !all(passes_escapee_filter)) {
        # Re-run EM on the cleaned gene set for sharper posteriors.
        n_removed <- sum(!passes_escapee_filter)
        logger::log_info("[{donor}] Removing {n_removed} escape genes after initial EM; re-running")

        gene_indices <- which(passes_escapee_filter)
        dat <- dat %>%
            dplyr::filter(gene %in% gene_indices) %>%
            dplyr::mutate(gene = match(gene, gene_indices))
        n_genes <- length(gene_indices)

        fits <- lapply(seq_len(n_inits), function(s) .run_em(dat, n_genes, init_seed = s))
        best <- fits[[which.max(sapply(fits, `[[`, "ll"))]]

        # After refit, all genes in dat passed the escapee filter, so mark them all TRUE
        passes_escapee_filter <- rep(TRUE, n_genes)
    }

    best$post <- best$post %>%
        dplyr::mutate(
            assignment = dplyr::case_when(
                post_X1 >= confidence_threshold ~ "X1",
                post_X1 <= 1 - confidence_threshold ~ "X2",
                TRUE ~ "unassigned"
            )
        )

    # gene_keep: logical of length nrow(original ref_mat), TRUE = gene survived
    # both the outlier filter and the post-convergence escapee filter
    gene_keep <- passes_outlier_filter
    if (refit_after_filter) {
        # Map filtered gene indices back to original space
        gene_keep[escaped_gene_indices] <- passes_escapee_filter
    } else {
        # Direct assignment works when no refit (passes_escapee_filter has same length as filtered genes)
        gene_keep[passes_outlier_filter] <- passes_escapee_filter
    }

    best$h_g <- best$h_g[passes_escapee_filter]
    best$pi_g <- best$pi_g[passes_escapee_filter]

    c(best, list(gene_keep = gene_keep))
}

#' @keywords internal
.compute_gene_llr <- function(dat, post, h_g, pi_g, rho) {
    p_if_X1 <- ifelse(h_g[dat$gene] == 0, pi_g[dat$gene], 1 - pi_g[dat$gene])
    p_if_X2 <- 1 - p_if_X1

    ll_X1 <- .loglik_obs(dat$ref, dat$n, p_if_X1, rho)
    ll_X2 <- .loglik_obs(dat$ref, dat$n, p_if_X2, rho)

    dat %>%
        dplyr::left_join(post %>% dplyr::select(cell, post_X1), by = "cell") %>%
        # Expected LLR per observation: weighted by P(X1-inactive) for each cell.
        # Positive contribution means the observation is consistent with current assignments.
        dplyr::mutate(llr = post_X1 * (ll_X1 - ll_X2)) %>%
        dplyr::group_by(gene) %>%
        dplyr::summarise(llr = sum(llr), .groups = "drop")
}

#' @keywords internal
.filter_escapee_genes <- function(dat, n_genes, h_g, pi_g, rho, post, mad_threshold = 2) {
    gene_llr <- .compute_gene_llr(dat, post, h_g, pi_g, rho)

    # Genes absent from dat have no observations above min_cov — exclude them
    keep_llr <- rep(FALSE, n_genes)
    keep_llr[gene_llr$gene[gene_llr$llr > 0]] <- TRUE # LLR > 0: data supports current assignment

    # Secondary pass: remove genes with unusually high escape fraction relative
    # to the rest of the sample. Guard against zero MAD (all pi_g identical).
    pi_mad <- stats::mad(pi_g)
    keep_pi <- if (pi_mad > 0) {
        # Robust z-score relative to the sample's own escape distribution
        (pi_g - stats::median(pi_g)) / pi_mad <= mad_threshold
    } else {
        rep(TRUE, n_genes) # all pi_g identical — no outliers possible
    }

    keep_llr & keep_pi
}

#' @keywords internal
.run_em <- function(dat, n_genes, max_iter = 50, tol = 1e-4, init_seed = 1) {
    h_g <- withr::with_seed(init_seed, sample(0:1, n_genes, replace = TRUE)) # random phase initialisation
    pi_g <- rep(0.05, n_genes) # start conservatively: assume 5% escape fraction
    rho <- 0.05 # fixed beta-binomial overdispersion

    # Unique kernel arguments are fixed across iterations; index them once.
    dedup <- .build_ll_dedup(dat)

    ll_prev <- -Inf
    for (iter in seq_len(max_iter)) {
        # Both escape orientations, computed once and shared by E- and M-steps.
        ll <- .betabinom_ll_both(dedup, pi_g, rho)
        post <- .e_step(dat, h_g, ll)
        h_g <- .m_step_phase(dat, post, h_g, ll)
        pi_g <- .m_step_pi(dat, post, h_g)
        # log-likelihood: sum of log(sigmoid(lor)) for each cell
        # Numerically stable: log(sigmoid(lor)) = pmin(0, lor) - log(1 + exp(-|lor|))
        ll_current <- sum(pmin(0, post$lor) - log1p(exp(-abs(post$lor))))
        if (abs(ll_current - ll_prev) < tol) {
            break
        }
        ll_prev <- ll_current
    }
    list(post = post, h_g = h_g, pi_g = pi_g, rho = rho, ll = ll_current)
}

#' @keywords internal
.e_step <- function(dat, h_g, ll, prior = 0.5) {
    # ll$L0 = loglik at pi_g[gene] (REF silenced), ll$L1 = loglik at 1 - pi_g[gene].
    # p(REF | X1 inactive): if h=0 (X1 carries REF), REF is silenced → pi_g → L0
    #                        if h=1 (X1 carries ALT), REF is active  → 1 - pi_g → L1
    h_row <- h_g[dat$gene] == 0
    # Per-observation X1-vs-X2 log-likelihood contrast. X2 is the complement
    # orientation, so the contrast is +(L0 - L1) when h=0 and -(L0 - L1) when h=1.
    obs_lor <- ifelse(h_row, ll$L0 - ll$L1, ll$L1 - ll$L0)

    logit_prior <- log(prior / (1 - prior)) # log(1) = 0 for equal prior
    # Σ_g [loglik(X1) - loglik(X2)] per cell, via rowsum (faster than group_by).
    cell_lor <- rowsum(obs_lor, dat$cell)
    lor <- as.numeric(cell_lor) + logit_prior

    tibble::tibble(
        cell = as.integer(rownames(cell_lor)),
        post_X1 = 1 / (1 + exp(-lor)), # sigmoid converts LOR to posterior
        lor = lor
    )
}

#' @keywords internal
.m_step_phase <- function(dat, post, h_g, ll) {
    # Expected log-likelihood under each phase: E_q[log p(ref | h, pi_g)], built
    # from the shared orientation pair. ll$L0 = loglik at pi_g (REF silenced),
    # ll$L1 = loglik at 1 - pi_g (REF active).
    # h=0: X1-inactive cells (weight post_X1) see REF fraction pi_g → L0;
    #       X2-inactive cells (weight 1-post_X1) see 1 - pi_g → L1.
    # h=1: the two orientations swap.
    # Scatter cell posteriors into a lookup indexed by cell id. post$cell is a
    # sorted subset of cell ids (cells with no covered gene are absent), so a
    # positional index would misalign — index by id instead.
    post_lookup <- numeric(max(dat$cell))
    post_lookup[post$cell] <- post$post_X1
    post_X1 <- post_lookup[dat$cell]

    obs_h0 <- post_X1 * ll$L0 + (1 - post_X1) * ll$L1
    obs_h1 <- post_X1 * ll$L1 + (1 - post_X1) * ll$L0

    ll_h0 <- rowsum(obs_h0, dat$gene)
    ll_h1 <- rowsum(obs_h1, dat$gene)
    genes <- as.integer(rownames(ll_h0))

    h_g_new <- h_g
    h_g_new[genes] <- as.integer(ll_h1 > ll_h0)
    h_g_new
}

#' @keywords internal
.m_step_pi <- function(dat, post, h_g, pi_bounds = c(0.001, 0.499)) {
    counts_with_posterior <- dat %>%
        dplyr::left_join(post %>% dplyr::select(cell, post_X1), by = "cell") %>%
        dplyr::mutate(
            # Soft-assigned inactive-allele read count:
            # h=0 (X1 carries REF): X1-inactive cells contribute REF, X2-inactive contribute ALT
            # h=1 (X1 carries ALT): roles reversed
            xi_ref_count = ifelse(
                h_g[gene] == 0,
                post_X1 * ref + (1 - post_X1) * alt,
                post_X1 * alt + (1 - post_X1) * ref
            ),
            xi_total = n
        )

    # MLE: inactive-allele fraction = (soft inactive reads) / (total reads)
    pi_new <- counts_with_posterior %>%
        dplyr::group_by(gene) %>%
        dplyr::summarise(pi = sum(xi_ref_count) / sum(xi_total), .groups = "drop")

    pi_g_new <- rep(0.05, length(h_g))
    # Clamp to (0.001, 0.499) so pi_g stays interpretable as a minor fraction
    pi_g_new[pi_new$gene] <- pmax(pi_bounds[1], pmin(pi_bounds[2], pi_new$pi))
    pi_g_new
}

#' @keywords internal
.loglik_obs <- function(ref, n, p, rho) {
    # Beta-binomial: rho > 0 adds overdispersion relative to binomial; rho = 0 reduces to binomial
    VGAM::dbetabinom(ref, size = n, prob = p, rho = rho, log = TRUE)
}

#' Beta-binomial log-density without the binomial coefficient
#'
#' Inlined equivalent of \code{VGAM::dbetabinom(ref, n, p, rho, log = TRUE)}
#' with the \code{lchoose(n, ref)} term dropped. That term depends only on the
#' data (not on \code{p}), so it cancels in every comparison the EM makes: the
#' E-step log-odds ratio is a difference of log-likelihoods at the same
#' \code{(ref, n)}, and the M-phase compares the two phase orientations at the
#' same rows. Dropping it avoids a per-row \code{lchoose} on every iteration
#' and matches the VGAM result up to that additive constant.
#'
#' The beta-binomial shape parameters relate to the mean \code{p} and the
#' correlation \code{rho} by \eqn{a = p (1 - rho) / rho} and
#' \eqn{b = (1 - p)(1 - rho) / rho}.
#'
#' @keywords internal
.betabinom_ll_kernel <- function(ref, n, p, rho) {
    scale <- (1 - rho) / rho
    a <- p * scale
    b <- (1 - p) * scale
    lbeta(ref + a, n - ref + b) - lbeta(a, b)
}

#' Precompute the deduplication index for the log-likelihood kernel
#'
#' Read depths are small integers, so the kernel arguments \code{(ref, n, gene)}
#' take only a few thousand distinct values across the hundreds of thousands of
#' observations in \code{dat} (often ~40x fewer). Because \code{dat} is constant
#' across EM iterations, we identify the unique argument triples once here and
#' reuse the mapping every iteration, evaluating the kernel on the unique rows
#' and scattering the result back.
#'
#' @return A list with \code{idx} (row -> unique-row index) and the unique-row
#'   slices \code{ref}, \code{n}, \code{gene}.
#'
#' @keywords internal
.build_ll_dedup <- function(dat) {
    idx <- as.integer(factor(paste(dat$ref, dat$n, dat$gene)))
    keep <- !duplicated(idx)
    ord <- order(idx[keep])
    list(
        idx = idx,
        ref = dat$ref[keep][ord],
        n = dat$n[keep][ord],
        gene = dat$gene[keep][ord]
    )
}

#' Log-likelihood of each observation under both escape orientations
#'
#' Returns, for every row of \code{dat}, the beta-binomial log-likelihood of
#' the REF count under \code{pi_g[gene]} (\code{L0}) and under
#' \code{1 - pi_g[gene]} (\code{L1}). These are the only two probabilities the
#' EM ever evaluates, so the E-step and both M-phase orientations are built
#' from this single pair rather than repeated \code{dbetabinom} calls.
#'
#' The kernel is evaluated only on the unique \code{(ref, n, gene)} rows
#' identified by \code{dedup} and scattered back to full length, avoiding the
#' large redundancy in the raw observations.
#'
#' @keywords internal
.betabinom_ll_both <- function(dedup, pi_g, rho) {
    p0 <- pi_g[dedup$gene]
    L0u <- .betabinom_ll_kernel(dedup$ref, dedup$n, p0, rho)
    L1u <- .betabinom_ll_kernel(dedup$ref, dedup$n, 1 - p0, rho)
    list(
        L0 = L0u[dedup$idx],
        L1 = L1u[dedup$idx]
    )
}

#' @keywords internal
.filter_outlier_genes <- function(ref_mat, alt_mat, min_cells = 10, min_cov = 1, mad_threshold = 2) {
    # Coverage per cell-gene pair
    n_mat <- ref_mat + alt_mat
    covered <- n_mat >= min_cov

    # Per gene: how many cells have sufficient coverage, and how many of those are REF-majority
    n_expressing <- rowSums(covered)
    ref_majority <- rowSums(covered & (ref_mat > alt_mat))

    # Drop genes seen in too few cells — not enough information to estimate allelic skew
    passes_count_filter <- n_expressing >= min_cells

    # Allelic skew: fraction of covered cells where REF > ALT, folded to [0.5, 1]
    # so that genes skewed toward either allele score equally high
    skew <- ref_majority[passes_count_filter] / n_expressing[passes_count_filter]
    skew <- pmax(skew, 1 - skew)

    # Robust z-score: genes with unusually extreme skew relative to the rest are
    # likely systematic (e.g. mapping bias, escape from XCI) rather than informative.
    # Guard against zero MAD (all skew values identical).
    skew_mad <- stats::mad(skew)
    passes_skew_filter <- if (skew_mad > 0) {
        z <- (skew - stats::median(skew)) / skew_mad
        abs(z) <= mad_threshold
    } else {
        rep(TRUE, length(skew)) # all skew identical — no outliers possible
    }

    keep <- rep(FALSE, nrow(ref_mat))
    keep[passes_count_filter] <- passes_skew_filter
    keep
}

#' @keywords internal
.pivot_counts_to_long <- function(ref_mat, alt_mat, min_cov = 3) {
    n_mat <- ref_mat + alt_mat
    # Matrix::which handles sparse lgCMatrix; base which() does not dispatch S4
    idx <- Matrix::which(n_mat >= min_cov, arr.ind = TRUE)
    tibble::tibble(
        gene = idx[, 1],
        cell = idx[, 2],
        ref = ref_mat[idx],
        alt = alt_mat[idx],
        n = n_mat[idx]
    )
}

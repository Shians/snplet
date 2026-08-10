#' Plot assignment heatmap from a SNPData object with stored XCI diagnostics
#'
#' Visualizes the REF allele fraction at the genes this donor's model retained.
#' Columns (cells or clonotypes) are ordered by posterior and annotated at the
#' top. Rows are ordered by discriminative power (the difference in mean REF
#' fraction between X1-active and X2-active cells) so the genes that most
#' cleanly flip allele between the two states appear at the bottom. A right-hand
#' annotation shows how many cells cover each gene. In the "Active X" annotation, units the
#' model could not call are split into "unassigned" (covered but low posterior
#' confidence) and "no coverage" (no covered informative SNPs, hence never
#' scored), so the two reasons are visually distinct. When any no-coverage units
#' are present they are further isolated into their own column slice, separated
#' from the scored units by a gap.
#'
#' @inheritSection assign_xci Phase is inferred from expression, not genotyped
#'
#' @param x A SNPData object, required, that had XCI diagnostics stored by
#'   \code{\link{assign_xci}} or \code{\link{assign_xci_by_clonotype}}.
#' @param donor Character scalar, required. Donor to visualize.
#' @param min_coverage_cells Integer (default 1, dropping only genes with no
#'   coverage in this donor). Minimum number of cells that must cover a gene
#'   for it to be shown; genes covered in fewer cells carry little signal and
#'   are dropped.
#' @param max_genes Integer, optional (default \code{NULL}, showing all
#'   retained genes). Maximum number of genes (rows) to display; when set,
#'   only the top \code{max_genes} most discriminative genes are shown.
#' @param show_gene_names Logical (default \code{TRUE}). Whether to draw gene
#'   names as row labels; useful to turn off when many genes are shown.
#' @param show_posterior Logical (default \code{TRUE}). Whether to draw the
#'   "P(active = X1)" posterior probability annotation row.
#' @param mark_boundaries Logical (default \code{TRUE}). Whether to draw
#'   dotted vertical lines at the transitions between assignment groups (X1
#'   -> unassigned -> X2).
#' @param show_unassigned Logical (default \code{TRUE}). Whether to include
#'   low-confidence unassigned columns: cells or clonotypes the model scored
#'   but could not confidently assign to X1 or X2.
#' @param show_no_coverage Logical (default \code{FALSE}). Whether to include
#'   no-coverage columns: cells or clonotypes with no covered informative
#'   SNPs, which the model never scored. These carry no data (all-NA columns).
#' @param cluster_rows Logical (default \code{FALSE}). Whether to cluster
#'   genes hierarchically instead of ordering them by discriminative power.
#' @param ref_fraction_palette Length-3 character vector of colours (default
#'   a blue-red diverging ramp). REF fraction heatmap body colours, mapped to
#'   fractions 0, 0.5, and 1.
#' @param assignment_palette Named character vector of colours (default
#'   green / purple), with names \code{"X1"}, \code{"X2"}, and
#'   \code{"unassigned"}. Assignment annotation colours.
#' @param posterior_palette Length-3 character vector of colours (default a
#'   brown-orange diverging ramp). posterior_X1_active annotation ramp
#'   colours, mapped to posteriors 0, 0.5, and 1. The three scales use
#'   distinct hue families (blue-red, green/purple, brown-orange) so viewers
#'   do not relate them.
#' @param na_fill Character scalar, a colour (default \code{"grey90"}, chosen
#'   to read distinctly from a balanced REF fraction). Fill for uncovered
#'   (missing) cells.
#'
#' @return A drawn \code{HeatmapList} (the donor is the plot title and the
#'   column axis is labelled with the modelling unit, "Cells" or "Clonotypes").
#'   Prints normally but, being already drawn, is not composable with \code{+}.
#'
#' @family X-chromosome inactivation functions
#' @importFrom circlize colorRamp2
#' @export
#'
#' @examples
#' \dontrun{
#' # Diagnostics stored by assign_xci() are read back automatically
#' snp_data <- assign_xci(snp_data)
#' plot_xci_heatmap(snp_data, donor = "donor1")
#' }
setGeneric(
    "plot_xci_heatmap",
    function(
        x,
        donor,
        min_coverage_cells = 1,
        max_genes = NULL,
        show_gene_names = TRUE,
        show_posterior = TRUE,
        mark_boundaries = TRUE,
        show_unassigned = TRUE,
        show_no_coverage = FALSE,
        cluster_rows = FALSE,
        ref_fraction_palette = c("#2166ac", "#f7f7f7", "#b2182b"),
        assignment_palette = c(X1 = "#1b7837", X2 = "#762a83", unassigned = "grey70"),
        posterior_palette = c("#8c510a", "#f7f7f7", "#e08214"),
        na_fill = "grey90"
    ) {
        standardGeneric("plot_xci_heatmap")
    }
)

#' @rdname plot_xci_heatmap
#' @include SNPData-class.R
setMethod(
    "plot_xci_heatmap",
    signature(x = "SNPData"),
    function(
        x,
        donor,
        min_coverage_cells = 1,
        max_genes = NULL,
        show_gene_names = TRUE,
        show_posterior = TRUE,
        mark_boundaries = TRUE,
        show_unassigned = TRUE,
        show_no_coverage = FALSE,
        cluster_rows = FALSE,
        ref_fraction_palette = c("#2166ac", "#f7f7f7", "#b2182b"),
        assignment_palette = c(X1 = "#1b7837", X2 = "#762a83", unassigned = "grey70"),
        posterior_palette = c("#8c510a", "#f7f7f7", "#e08214"),
        na_fill = "grey90"
    ) {
        barcode_info <- barcode_info(x)
        if (!"active_x" %in% colnames(barcode_info) || !.has_xci_diagnostics(x)) {
            stop("No stored XCI diagnostics found. Run assign_xci(x) first.")
        }

        donor_data <- filter_samples(x, donor == !!donor)
        # Restrict to the genes this donor's model actually retained, not the union
        # of informative genes across all donors.
        donor_informative_snps <- donor_snp_info(x) %>%
            dplyr::filter(donor == !!donor, xci_informative) %>%
            dplyr::pull(snp_id)
        donor_data <- filter_snps(donor_data, snp_id %in% donor_informative_snps)

        filtered_snp_info <- snp_info(donor_data)
        donor_barcode_info <- barcode_info(donor_data)

        # Cells this donor's model actually assigned (drop NA active_x)
        assigned <- donor_barcode_info %>%
            dplyr::mutate(
                assignment = ifelse(is.na(active_x), "unassigned", active_x),
                unit_id = cell_id
            )

        ref_mat <- ref_count(donor_data)
        alt_mat <- alt_count(donor_data)
        rownames(ref_mat) <- filtered_snp_info$snp_id
        rownames(alt_mat) <- filtered_snp_info$snp_id

        # For a clonotype-level fit, aggregate the per-cell counts and assignments up
        # to the modelling unit (clonotype) so the heatmap columns match how the
        # model actually saw the data.
        fit_unit <- unique(donor_barcode_info$xci_fit_unit %||% "cell")
        by_clonotype <- "clonotype" %in% fit_unit && "clonotype" %in% colnames(assigned)
        if (by_clonotype) {
            agg <- .aggregate_xci_to_clonotype(assigned, ref_mat, alt_mat)
            ref_mat <- agg$ref_mat
            alt_mat <- agg$alt_mat
            assigned <- agg$assigned
        }

        # Optionally drop the two kinds of uncalled columns independently:
        # no-coverage units have an NA posterior (never scored); low-confidence
        # unassigned units were scored but not confidently assigned.
        no_coverage <- is.na(assigned$xci_post_X1_active)
        low_confidence <- assigned$assignment == "unassigned" & !no_coverage
        drop_cols <- (!show_unassigned & low_confidence) | (!show_no_coverage & no_coverage)
        if (any(drop_cols)) {
            assigned <- assigned[!drop_cols, , drop = FALSE]
            ref_mat <- ref_mat[, assigned$unit_id, drop = FALSE]
            alt_mat <- alt_mat[, assigned$unit_id, drop = FALSE]
        }

        gene_name_map <- stats::setNames(filtered_snp_info$gene_name, filtered_snp_info$snp_id)
        .plot_xci_heatmap_from_parts(
            ref_mat = ref_mat,
            alt_mat = alt_mat,
            assignment = assigned$assignment,
            post_X1_active = assigned$xci_post_X1_active,
            unit_ids = assigned$unit_id,
            gene_name_map = gene_name_map,
            min_coverage_cells = min_coverage_cells,
            max_genes = max_genes,
            show_gene_names = show_gene_names,
            show_posterior = show_posterior,
            mark_boundaries = mark_boundaries,
            cluster_rows = cluster_rows,
            ref_fraction_palette = ref_fraction_palette,
            assignment_palette = assignment_palette,
            posterior_palette = posterior_palette,
            na_fill = na_fill,
            title = donor,
            unit_label = if (by_clonotype) "clonotypes" else "cells"
        )
    }
)

#' Aggregate per-cell XCI counts and assignments up to the clonotype level
#'
#' Sums ALT/REF counts across the cells of each clonotype and collapses the
#' per-cell assignments (which are identical within a clonotype for a
#' clonotype-level fit) to one row per clonotype.
#'
#' @keywords internal
.aggregate_xci_to_clonotype <- function(assigned, ref_mat, alt_mat) {
    # Cells sharing a clonotype form the aggregation groups, in a stable order.
    clonotype <- assigned$clonotype
    keep <- !is.na(clonotype)
    assigned <- assigned[keep, , drop = FALSE]
    clonotype <- clonotype[keep]
    cell_ids <- assigned$cell_id

    ref_mat <- ref_mat[, cell_ids, drop = FALSE]
    alt_mat <- alt_mat[, cell_ids, drop = FALSE]

    # Group indicator matrix (cells x clonotypes); counts sum within group.
    clonotype <- factor(clonotype, levels = unique(clonotype))
    group_mat <- Matrix::sparse.model.matrix(~ 0 + clonotype)
    colnames(group_mat) <- levels(clonotype)
    ref_agg <- ref_mat %*% group_mat
    alt_agg <- alt_mat %*% group_mat

    # One assignment row per clonotype (constant within a clonotype-level fit).
    assigned_agg <- assigned %>%
        dplyr::distinct(clonotype, .keep_all = TRUE) %>%
        dplyr::mutate(unit_id = as.character(clonotype)) %>%
        dplyr::arrange(match(unit_id, levels(clonotype)))

    list(
        ref_mat = as.matrix(ref_agg),
        alt_mat = as.matrix(alt_agg),
        assigned = assigned_agg
    )
}

#' @keywords internal
.plot_xci_heatmap_from_parts <- function(
    ref_mat,
    alt_mat,
    assignment,
    post_X1_active,
    unit_ids,
    gene_name_map,
    title = NULL,
    min_coverage_cells = 1,
    max_genes = NULL,
    show_gene_names = TRUE,
    show_posterior = TRUE,
    mark_boundaries = TRUE,
    cluster_rows = FALSE,
    ref_fraction_palette = c("#2166ac", "#f7f7f7", "#b2182b"),
    assignment_palette = c(X1 = "#1b7837", X2 = "#762a83", unassigned = "grey70"),
    posterior_palette = c("#8c510a", "#f7f7f7", "#e08214"),
    na_fill = "grey90",
    unit_label = "cells"
) {
    # REF allele fraction; NA where uncovered
    cov_mat <- ref_mat + alt_mat
    frac_mat <- ref_mat / cov_mat
    frac_mat[cov_mat == 0] <- NA

    # Replace SNP ID row names with gene names
    rownames(frac_mat) <- gene_name_map[rownames(frac_mat)]

    # Per-gene coverage: number of cells with any reads at this SNP. Drop genes
    # covered in fewer than min_coverage_cells cells — they carry little signal
    # and dominate the plot with grey.
    covered_cells <- rowSums(cov_mat > 0)
    row_keep <- covered_cells >= min_coverage_cells
    frac_mat <- frac_mat[row_keep, , drop = FALSE]
    covered_cells <- covered_cells[row_keep]

    # A unit the model never scored (NA posterior) had no covered informative
    # SNPs; it is unassigned for lack of data rather than lack of confidence.
    # These sort to the far end of the plot (NA posteriors order last).
    no_coverage <- is.na(post_X1_active)

    # Order units: X1-active → X2-active → unassigned → no-coverage (by
    # posterior, NAs last). Dense matrix from here so rowMeans(na.rm) and
    # subsetting behave predictably.
    idx <- order(post_X1_active, decreasing = TRUE)
    frac_mat <- as.matrix(frac_mat[, unit_ids[idx], drop = FALSE])
    assignment <- assignment[idx]
    post_X1_active <- post_X1_active[idx]
    no_coverage <- no_coverage[idx]

    # Order genes by discriminative power rather than clustering a mostly-NA
    # matrix: the difference in mean REF fraction between X1-active and
    # X2-active cells. Genes that cleanly flip allele between the two states
    # rise to the top.
    #
    # Naively taking abs(mean_X1 - mean_X2) lets a gene covered in only a
    # handful of cells score a perfect flip and outrank a well-covered gene with
    # a little biological noise. Shrink each group mean toward 0.5 by its cell
    # count (a Beta(k0/2, k0/2) prior) so a clean flip only counts when it is
    # backed by enough observations.
    x1_cols <- assignment == "X1"
    x2_cols <- assignment == "X2"
    shrunk_mean <- function(m, k0 = 5) {
        if (ncol(m) == 0) {
            return(rep(0.5, nrow(m)))
        }
        n <- rowSums(!is.na(m))
        s <- rowSums(m, na.rm = TRUE)
        # posterior mean of a Beta(k0/2, k0/2) prior with n observations
        (s + k0 / 2) / (n + k0)
    }
    mean_x1 <- shrunk_mean(frac_mat[, x1_cols, drop = FALSE])
    mean_x2 <- shrunk_mean(frac_mat[, x2_cols, drop = FALSE])
    discrimination <- abs(mean_x1 - mean_x2)
    gene_order <- order(discrimination, decreasing = TRUE)
    frac_mat <- frac_mat[gene_order, , drop = FALSE]
    covered_cells <- covered_cells[gene_order]

    # Keep only the top max_genes most discriminative genes when requested.
    if (!is.null(max_genes) && max_genes < nrow(frac_mat)) {
        top <- seq_len(max_genes)
        frac_mat <- frac_mat[top, , drop = FALSE]
        covered_cells <- covered_cells[top]
    }

    # Distinguish the two reasons a unit is unassigned: low posterior confidence
    # (has coverage but no clear call) versus no covered informative SNPs at all.
    # The latter carries no data and is labelled separately.
    assignment[assignment == "unassigned" & no_coverage] <- "no coverage"

    # Ensure the palette carries a colour for the "no coverage" category. Use a
    # dark slate distinct from the light grey of low-confidence "unassigned" so
    # the two kinds of unassigned unit are easy to tell apart in the annotation.
    if (!"no coverage" %in% names(assignment_palette)) {
        assignment_palette <- c(assignment_palette, "no coverage" = "grey30")
    }

    # Assemble the column annotation. The posterior_X1_active row is optional.
    # Internal names stay code-friendly for the colour mapping; display labels
    # are set via annotation_label.
    ann_args <- list(assignment = assignment)
    ann_col <- list(assignment = assignment_palette)
    ann_label <- c(assignment = "Active X")
    if (show_posterior) {
        ann_args$posterior_X1_active <- post_X1_active
        ann_col$posterior_X1_active <- circlize::colorRamp2(c(0, 0.5, 1), posterior_palette)
        ann_label <- c(ann_label, posterior_X1_active = "P(active = X1)")
    }
    col_ann <- do.call(
        ComplexHeatmap::HeatmapAnnotation,
        c(
            ann_args,
            list(col = ann_col, annotation_label = ann_label, show_annotation_name = TRUE)
        )
    )

    # Row annotation: how many units (cells or clonotypes) cover each gene, so
    # the reader can weight sparse rows appropriately.
    row_ann_args <- list(
        ComplexHeatmap::anno_barplot(
            covered_cells,
            which = "row",
            gp = grid::gpar(fill = "grey40", col = NA),
            width = grid::unit(1.5, "cm")
        )
    )
    names(row_ann_args) <- paste0(unit_label, "_covered")
    row_ann <- do.call(
        ComplexHeatmap::rowAnnotation,
        c(row_ann_args, show_annotation_name = TRUE, annotation_name_rot = 0)
    )

    # The modelling unit (Cells / Clonotypes) labels the column axis along the
    # bottom.
    axis_label <- stringr::str_to_title(unit_label)

    # Diverging colour ramp for the REF fraction body, anchored at 0 / 0.5 / 1.
    ref_col <- circlize::colorRamp2(c(0, 0.5, 1), ref_fraction_palette)

    # Physically isolate the no-coverage units with a column split when any
    # exist. The larger scored slice keeps the modelling-unit axis label; the
    # small slice is titled "no coverage".
    has_no_coverage <- any(no_coverage)
    column_split <- if (has_no_coverage) {
        factor(
            ifelse(no_coverage, "no coverage", "scored"),
            levels = c("scored", "no coverage")
        )
    } else {
        NULL
    }
    column_title <- if (has_no_coverage) c(axis_label, "no coverage") else axis_label

    # Optionally mark the X1 -> unassigned and unassigned -> X2 transitions with
    # dotted vertical lines. Boundaries are the fractional x-positions of the
    # assignment changes within the scored columns; the scored columns are the
    # first slice (index 1) whether or not a no-coverage slice exists. A
    # layer_fun bakes the lines into the heatmap so they survive re-drawing.
    layer_fun <- NULL
    if (mark_boundaries) {
        scored_assignment <- assignment[!no_coverage]
        n_scored <- length(scored_assignment)
        changes <- which(scored_assignment[-1] != scored_assignment[-n_scored])
        boundaries <- changes / n_scored
        if (length(boundaries) > 0) {
            layer_fun <- function(j, i, x, y, w, h, fill, slice_r, slice_c) {
                # Only the scored slice (first column slice) carries boundaries.
                if (slice_c != 1) {
                    return(invisible(NULL))
                }
                for (b in boundaries) {
                    grid::grid.lines(
                        x = grid::unit(c(b, b), "npc"),
                        y = grid::unit(c(0, 1), "npc"),
                        gp = grid::gpar(col = "black", lty = "dotted", lwd = 1.5)
                    )
                }
            }
        }
    }

    ht <- ComplexHeatmap::Heatmap(
        as.matrix(frac_mat),
        col = ref_col,
        top_annotation = col_ann,
        right_annotation = row_ann,
        cluster_columns = FALSE,
        cluster_column_slices = FALSE,
        column_split = column_split,
        cluster_rows = cluster_rows,
        show_row_names = show_gene_names,
        show_column_names = FALSE,
        name = "REF fraction",
        layer_fun = layer_fun,
        # A darker grey for uncovered cells (na_fill) so absence of data reads
        # distinctly from a balanced (REF fraction ~ 0.5) measurement.
        na_col = na_fill,
        column_title = column_title,
        column_title_side = "bottom",
        column_gap = grid::unit(2, "mm")
    )

    # Draw so the donor can sit as the overall plot title at the top while the
    # bottom column_title serves as the x-axis label. Returns a drawn
    # HeatmapList (still prints normally; no longer composable with `+`/`%v%`).
    ComplexHeatmap::draw(
        ht,
        column_title = title,
        column_title_gp = grid::gpar(fontsize = 14, fontface = "bold")
    )
}

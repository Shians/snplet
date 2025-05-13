#' Plot gene annotation track
#'
#' This function creates a gene annotation track visualization showing gene positions and names
#' for a genomic region of interest.
#'
#' @param gene_anno A data frame containing gene annotation information with columns 'pos' (position)
#'        and 'gene_name'
#' @param x_range A numeric vector of length 2 specifying the genomic range limits to display
#' @return A ggplot object representing the gene annotation track
#' @export
#'
#' @examples
#' \dontrun{
#' gene_anno <- data.frame(pos = c(1000, 2000, 3000), gene_name = c("GENE1", "GENE2", "GENE3"))
#' plot_gene_anno_track(gene_anno, x_range = c(500, 3500))
#' }
plot_gene_anno_track <- function(gene_anno, x_range) {
    ggplot2::ggplot(gene_anno, ggplot2::aes(x = pos, y = 0, label = gene_name)) +
        ggplot2::geom_segment(ggplot2::aes(xend = pos, yend = 0.15)) +
        ggrepel::geom_text_repel(
            ggplot2::aes(y = 0.15),
            force_pull = 0,
            nudge_y = 0.5,
            direction = "x",
            angle = 90,
            hjust = 1
        ) +
        ggplot2::scale_x_continuous(limits = x_range, expand = c(0, 0)) +
        ggplot2::theme_void() +
        ggplot2::theme(
            axis.line.y = ggplot2::element_blank(),
            axis.ticks.y = ggplot2::element_blank(),
            axis.text.y = ggplot2::element_blank(),
            axis.title.y = ggplot2::element_blank()
        )
}

#' Plot minor allele frequency track
#'
#' This function creates a visualization of minor allele frequencies (MAF) across
#' a genomic region, with optional faceting by a grouping variable.
#'
#' @param allele_count_df A data frame containing SNP information with columns 'pos' (position),
#'        'maf' (minor allele frequency), and 'donor_id' for faceting
#' @param facet A variable to use for faceting the plot
#' @param x_range A numeric vector of length 2 specifying the genomic range limits to display
#' @return A ggplot object showing MAF values across the genomic region
#' @export
#'
#' @examples
#' \dontrun{
#' allele_df <- get_allele_counts(snp_data)
#' plot_maf_track(allele_df, facet = donor_id, x_range = c(1000000, 2000000))
#' }
plot_maf_track <- function(allele_count_df, facet, x_range) {
    allele_count_df %>%
        ggplot2::ggplot(aes(x = pos, y = maf)) +
        ggplot2::geom_hline(yintercept = 0.1, linetype = "dashed", color = "red", size = 0.5) +
        ggplot2::geom_point(size = 1, alpha = 0.10) +
        ggplot2::facet_grid(rows = ggplot2::vars(donor_id)) +
        ggplot2::scale_x_continuous(
            labels = scales::label_number(scale_cut = scales::cut_short_scale(), suffix = "b"),
            limits = x_range,
            expand = c(0, 0)
        ) +
        ggplot2::scale_y_continuous(
            breaks = c(0, 0.5),
            labels = c("0", "0.5"),
            expand = c(0, 0)
        ) +
        ggplot2::coord_cartesian(clip = 'off') +
        ggplot2::theme_classic() +
        ggplot2::theme(
            panel.grid = ggplot2::element_blank(),
            axis.line.x = ggplot2::element_blank(),
            text = ggplot2::element_text(size = 16),
            panel.spacing = grid::unit(1.5, "lines")
        ) +
        ggplot2::labs(
            title = "Clonotype level MAF",
            x = "Position",
            y = "Minor Allele Frequency"
        )
}

#' Plot p-value track for MAF significance
#'
#' This function creates a Manhattan-style plot showing the significance (-log10 p-value)
#' of minor allele frequency differences across a genomic region.
#'
#' @param allele_counts_df A data frame containing SNP information with columns 'pos' (position),
#'        'adj_p_val' (adjusted p-value), and a column for faceting
#' @param facet A variable to use for faceting the plot, supplied as a bare column name
#' @param x_range A numeric vector of length 2 specifying the genomic range limits to display
#' @return A ggplot object showing -log10(adjusted p-value) across the genomic region
#' @export
#'
#' @examples
#' \dontrun{
#' snp_stats <- calculate_snp_significance(snp_data)
#' plot_maf_pval_track(snp_stats, facet = donor_id, x_range = c(1000000, 2000000))
#' }
plot_maf_pval_track <- function(allele_counts_df, facet, x_range) {
    allele_counts_df %>%
        ggplot2::ggplot(ggplot2::aes(x = pos, y = -log10(adj_p_val))) +
        ggplot2::geom_hline(yintercept = 1.3, linetype = "dashed", color = "red") +
        ggplot2::geom_vline(
            data = signif_snps_clonotype,
            ggplot2::aes(xintercept = pos),
            linetype = "dotted",
            color = "gray") +
        ggplot2::geom_point(size = 1) +
        ggplot2::facet_grid(rows = vars({{facet}})) +
        ggplot2::scale_x_continuous(
            labels = scales::label_number(scale_cut = scales::cut_short_scale(), suffix = "b"),
            limits = x_range,
            expand = c(0, 0)
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(
            panel.grid = ggplot2::element_blank(),
            text = ggplot2::element_text(size = 16)
        ) +
        ggplot2::labs(
            x = "Position",
            y = "-log10(adj-p-value)"
        )
}

#' Plot distribution of library size across samples
#'
#' This function visualizes the distribution of library sizes (total read counts)
#' across all samples in the SNP data object using a density plot.
#'
#' @param snp_data A SNPData object containing sample information with library size data
#' @return A ggplot object showing the density distribution of library sizes
#' @export
#'
#' @examples
#' \dontrun{
#' snp_data <- read_cellsnp("/path/to/cellsnp")
#' plot_lib_size_distribution(snp_data)
#' }
plot_lib_size_distribution <- function(snp_data) {
    get_sample_info(snp_data) %>%
        ggplot(aes(x = library_size)) +
        geom_density() +
        scale_y_continuous(expand = expansion(0, 0)) +
        theme_density() +
        labs(
            x = "Library size",
            y = "Density",
            title = "Distribution of library size"
        )
}

#' Plot distribution of SNP coverage
#'
#' This function visualizes the distribution of coverage values across all SNPs
#' in the SNP data object using a density plot with a log10-transformed x-axis.
#'
#' @param snp_data A SNPData object containing SNP information with coverage data
#' @return A ggplot object showing the density distribution of SNP coverage on a log10 scale
#' @export
#'
#' @examples
#' \dontrun{
#' snp_data <- read_cellsnp("/path/to/cellsnp")
#' plot_snp_cov_distribution(snp_data)
#' }
plot_snp_cov_distribution <- function(snp_data) {
    get_snp_info(snp_data) %>%
        ggplot(aes(x = coverage)) +
        geom_density() +
        scale_y_continuous(expand = expansion(0, 0)) +
        scale_x_log10(labels = scales::label_comma()) +
        theme_density() +
        labs(
            x = "Coverage",
            y = "Density",
            title = "Distribution of SNP coverage"
        )
}

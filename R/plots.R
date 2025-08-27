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
#' snp_data <- get_example_snpdata()
#' plot_lib_size_distribution(snp_data)
#' }
plot_lib_size_distribution <- function(snp_data) {
    get_barcode_info(snp_data) %>%
        ggplot2::ggplot(ggplot2::aes(x = library_size)) +
        ggplot2::geom_density() +
        ggplot2::scale_y_continuous(expand = ggplot2::expansion(0, 0)) +
        theme_density() +
        ggplot2::labs(
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
        ggplot2::ggplot(ggplot2::aes(x = coverage)) +
        ggplot2::geom_density() +
        ggplot2::scale_y_continuous(expand = ggplot2::expansion(0, 0)) +
        ggplot2::scale_x_log10(labels = scales::label_comma()) +
        theme_density() +
        ggplot2::labs(
            x = "Coverage",
            y = "Density",
            title = "Distribution of SNP coverage"
        )
}

#' Plot distribution of minor allele frequency (MAF)
#'
#' This function visualizes the distribution of minor allele frequencies (MAF)
#' in a given data frame using a density plot.
#'
#' @param df A data frame containing a column named 'maf' with MAF values
#' @return A ggplot object showing the density distribution of MAF
#' @export
#'
#' @examples
#' df <- data.frame(maf = runif(1000, 0, 0.5))
#' plot_maf_distribution(df)
plot_maf_distribution <- function(df) {
    ggplot2::ggplot(df, ggplot2::aes(x = maf)) +
        ggplot2::geom_density() +
        ggplot2::scale_x_continuous(expand = ggplot2::expansion(0, 0)) +
        ggplot2::scale_y_continuous(expand = ggplot2::expansion(0, 0)) +
        theme_density() +
        ggplot2::labs(
            x = "Minor Allele Frequency (MAF)",
            y = "Density",
            title = "Distribution of Minor Allele Frequency"
        )
}

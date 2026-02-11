#' Assign inactive X chromosome to cells
#'
#' Identifies which X chromosome is inactive in female cells based on allelic
#' imbalance patterns at heterozygous SNPs. Uses hierarchical clustering to
#' assign cells into two groups representing the two possible X-inactivation
#' states.
#'
#' @details
#' X-chromosome inactivation (XCI) is a dosage compensation mechanism in female
#' mammals where one of the two X chromosomes is randomly silenced in each cell.
#' This function infers which X is inactive by analyzing allelic expression
#' patterns at heterozygous SNPs on the X chromosome.
#'
#' The algorithm works as follows:
#' \enumerate{
#'   \item Filters to heterozygous SNPs on the X chromosome for each donor
#'   \item Selects the SNP with highest coverage per gene to avoid redundancy
#'   \item Identifies informative SNPs with sufficient coverage across cells
#'   \item Computes an expression-like matrix capturing allelic imbalance
#'   \item Performs hierarchical clustering on the sign of allelic imbalance
#'   \item Assigns cells to two clusters (X1 or X2) representing inactive X states
#' }
#'
#' Cells within each cluster should show consistent allelic bias patterns,
#' reflecting which parental X chromosome is active vs. inactive.
#'
#' @param x SNPData object containing X chromosome SNP data with donor
#'   assignments and heterozygosity information
#'
#' @return SNPData object with an additional \code{inactive_x} column in
#'   barcode metadata, with values "X1" or "X2" indicating the inferred
#'   inactive X chromosome state
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
setGeneric("assign_inactive_x", function(x) standardGeneric("assign_inactive_x"))

#' @rdname assign_inactive_x
#' @include SNPData-class.R
setMethod("assign_inactive_x", signature(x = "SNPData"), function(x) {
    unique_donors <- unique(get_barcode_info(x)$donor)

    unique_donor_snp_data <- map(
        unique_donors,
        ~filter_samples(x, donor == .x)
    ) %>%
        magrittr::set_names(unique_donors)

    # assign inactive x for each donor and combine results
    barcode_inactive_x_df <- unique_donor_snp_data %>%
        map(.assign_inactive_x_single_donor) %>%
        bind_rows()

    # add inactive_x to barcode_info
    x %>%
        add_barcode_metadata(barcode_inactive_x_df)
})

.assign_inactive_x_single_donor <- function(
    snp_data
) {
    expr_matrix <- .prepare_expr_matrix(
        snp_data,
        min_coverage = 2,
        min_sample_prop = 0.01
    )

    .assign_cluster(expr_matrix, n_clusters = 2)
}

.assign_inactive_x_single_donor_by_clonotype <- function(
    snp_data
) {
    donor <- unique(get_barcode_info(snp_data)$donor)

    expr_matrix <- .prepare_expr_matrix_by_clonotype(
        snp_data,
        min_coverage = 2,
        min_sample_prop = 0.01
    )

    # Assign clusters to clonotypes
    clonotype_cluster_df <- .assign_cluster(expr_matrix, n_clusters = 2) %>%
        dplyr::rename(clonotype = cell_id)

    # Project clonotype assignments back to cells
    barcode_info <- get_barcode_info(snp_data)
    cell_cluster_df <- dplyr::select(barcode_info, cell_id, clonotype) %>%
        dplyr::left_join(
            clonotype_cluster_df,
            by = c("clonotype"),
            multiple = "first"
        ) %>%
        dplyr::select(cell_id, inactive_x)

    logger::log_success("Inactive X assignment complete for donor {donor}")
    cell_cluster_df
}

.get_informative_snps <- function(coverage_mat, min_coverage, min_sample_prop) {
    barcodes_with_sufficient_coverage <- rowSums(coverage_mat >= min_coverage)
    barcode_threshold <- ncol(coverage_mat) * min_sample_prop
    barcodes_with_sufficient_coverage >= barcode_threshold
}

.assign_cluster <- function(expr_mat, n_clusters) {
    dist_mat <- dist(t(as.matrix(expr_mat)), method = "euclidean")
    hc <- hclust(dist_mat, method = "ward.D2")
    cluster_assignments <- cutree(hc, k = n_clusters)
    tibble(
        cell_id = names(cluster_assignments),
        inactive_x = paste0("X", cluster_assignments)
    )
}

.prepare_expr_matrix <- function(
    snp_data,
    min_coverage = 2,
    min_sample_prop = 0.01
) {
    donor <- unique(get_barcode_info(snp_data)$donor)
    logger::log_info("Assigning inactive X for donor {donor}")

    # filter to heterozygous snps only
    het_snp_df <- snp_data %>%
        donor_het_status_df() %>%
        filter(zygosity == "het")

    het_snp_ids <- het_snp_df %>%
        pull(snp_id)

    logger::log_info("Filtering to heterozygous SNPs")
    snp_data <- snp_data %>%
        filter_snps(snp_id %in% het_snp_ids)

    # filter to top snp per gene
    top_snp_per_gene <- get_snp_info(snp_data) %>%
        arrange(desc(coverage)) %>%
        slice_head(n = 1, by = "gene_name")

    logger::log_info("Filtering to top SNP per gene")
    snp_data <- snp_data %>%
        filter_snps(snp_id %in% top_snp_per_gene$snp_id)

    # get expression matrix
    coverage_mat <- snplet::coverage(snp_data)

    min_coverage <- 2
    min_sample_prop <- 0.01
    logger::log_info("Selecting informative SNPs with at least {min_coverage} coverage in at least {min_sample_prop * 100}% of cells")
    informative_snps <- .get_informative_snps(
        coverage_mat,
        min_coverage = 2,
        min_sample_prop = 0.01
    )

    coverage_mat <- snplet::coverage(snp_data)
    informative_snps <- .get_informative_snps(
        coverage_mat,
        min_coverage = min_coverage,
        min_sample_prop = min_sample_prop
    )
    logger::log_info("Using {sum(informative_snps)} informative SNPs for inactive X assignment")
    expr_matrix <- to_expr_matrix(snp_data)[informative_snps, , drop = FALSE]
    nonzero_cells <- colSums(abs(expr_matrix)) > 0
    logger::log_info("Retaining {sum(nonzero_cells)} barcodes with non-zero expression across informative SNPs")
    expr_matrix[, nonzero_cells, drop = FALSE]
}

.prepare_expr_matrix_by_clonotype <- function(
    snp_data,
    min_coverage = 2,
    min_sample_prop = 0.01
) {
    donor <- unique(get_barcode_info(snp_data)$donor)
    logger::log_info("Assigning inactive X for donor {donor} at clonotype level")

    # filter to heterozygous snps only
    het_snp_df <- snp_data %>%
        donor_het_status_df() %>%
        filter(zygosity == "het")

    het_snp_ids <- het_snp_df %>%
        pull(snp_id)

    logger::log_info("Filtering to heterozygous SNPs")
    snp_data <- snp_data %>%
        filter_snps(snp_id %in% het_snp_ids)

    # filter to top snp per gene
    top_snp_per_gene <- get_snp_info(snp_data) %>%
        arrange(desc(coverage)) %>%
        slice_head(n = 1, by = "gene_name")

    logger::log_info("Filtering to top SNP per gene")
    snp_data <- snp_data %>%
        filter_snps(snp_id %in% top_snp_per_gene$snp_id)

    # Get clonotype grouping (remove NA clonotypes)
    barcode_info <- get_barcode_info(snp_data)
    clonotypes <- barcode_info$clonotype

    # Filter out cells with NA clonotypes
    non_na_clonotypes <- !is.na(clonotypes)
    if (!any(non_na_clonotypes)) {
        stop("No cells with non-NA clonotype values for donor {donor}")
    }

    snp_data_filtered <- snp_data[, non_na_clonotypes]
    clonotypes_filtered <- clonotypes[non_na_clonotypes]

    logger::log_info("Aggregating counts by clonotype ({length(unique(clonotypes_filtered))} unique clonotypes)")

    # Aggregate counts by clonotype
    ref_count_by_clonotype <- groupedRowSums(ref_count(snp_data_filtered), clonotypes_filtered)
    alt_count_by_clonotype <- groupedRowSums(alt_count(snp_data_filtered), clonotypes_filtered)

    # Calculate expression matrix (SNPs × clonotypes)
    # expr = (REF - ALT) / (REF + ALT)
    total_count <- ref_count_by_clonotype + alt_count_by_clonotype
    expr_matrix <- (ref_count_by_clonotype - alt_count_by_clonotype) / total_count

    # Replace NaN values (0/0) with 0
    expr_matrix[is.nan(expr_matrix)] <- 0

    # Get coverage for informative SNP selection
    coverage_mat <- total_count

    logger::log_info("Selecting informative SNPs with at least {min_coverage} coverage in at least {min_sample_prop * 100}% of clonotypes")
    informative_snps <- .get_informative_snps(
        coverage_mat,
        min_coverage = min_coverage,
        min_sample_prop = min_sample_prop
    )

    logger::log_info("Using {sum(informative_snps)} informative SNPs for inactive X assignment")
    expr_matrix_filtered <- expr_matrix[informative_snps, , drop = FALSE]

    # Filter to clonotypes with non-zero expression
    nonzero_clonotypes <- colSums(abs(expr_matrix_filtered)) > 0
    logger::log_info("Retaining {sum(nonzero_clonotypes)} clonotypes with non-zero expression across informative SNPs")
    expr_matrix_filtered[, nonzero_clonotypes, drop = FALSE]
}

#' Assign inactive X chromosome to cells by clonotype
#'
#' Identifies which X chromosome is inactive in female cells based on allelic
#' imbalance patterns at heterozygous SNPs, aggregating to the clonotype level
#' before assigning and then projecting back to individual cells.
#'
#' @details
#' This function is similar to \code{\link{assign_inactive_x}} but performs
#' clustering at the clonotype level rather than the cell level. This approach:
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
#'   \item Identifies informative SNPs with sufficient coverage across clonotypes
#'   \item Computes an expression-like matrix capturing allelic imbalance at clonotype level
#'   \item Performs hierarchical clustering on the sign of allelic imbalance
#'   \item Assigns clonotypes to two clusters (X1 or X2) representing inactive X states
#'   \item Projects clonotype assignments back to individual cells
#' }
#'
#' @param x SNPData object containing X chromosome SNP data with donor
#'   assignments, clonotype information, and heterozygosity information
#'
#' @return SNPData object with an additional \code{inactive_x} column in
#'   barcode metadata, with values "X1" or "X2" indicating the inferred
#'   inactive X chromosome state
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
setGeneric("assign_inactive_x_by_clonotype", function(x) standardGeneric("assign_inactive_x_by_clonotype"))

#' @rdname assign_inactive_x_by_clonotype
#' @include SNPData-class.R
setMethod("assign_inactive_x_by_clonotype", signature(x = "SNPData"), function(x) {
    # Check if clonotype column exists
    barcode_info <- get_barcode_info(x)
    if (!"clonotype" %in% colnames(barcode_info)) {
        stop("Clonotype information not available. Add clonotype data using add_barcode_metadata() or import_cellsnp() with vdj_file parameter.")
    }

    # Check if all clonotype values are NA
    if (all(is.na(barcode_info$clonotype))) {
        stop("All clonotype values are NA. Cannot perform clonotype-level X-inactivation assignment. Add clonotype data using add_barcode_metadata() or import_cellsnp() with vdj_file parameter.")
    }

    unique_donors <- unique(get_barcode_info(x)$donor)

    unique_donor_snp_data <- map(
        unique_donors,
        ~filter_samples(x, donor == .x)
    ) %>%
        magrittr::set_names(unique_donors)

    # assign inactive x for each donor and combine results
    barcode_inactive_x_df <- unique_donor_snp_data %>%
        map(.assign_inactive_x_single_donor_by_clonotype) %>%
        bind_rows()

    # add inactive_x to barcode_info
    x %>%
        add_barcode_metadata(barcode_inactive_x_df)
})

#' Plot inactive X chromosome heatmap for a donor
#'
#' Visualizes allelic imbalance patterns at informative heterozygous SNPs to
#' reveal X-inactivation structure. Cells are clustered based on their allelic
#' bias signatures, which should separate into two groups corresponding to the
#' two inactive X states.
#'
#' @details
#' This function creates a heatmap showing the sign of allelic expression at
#' heterozygous SNPs across cells from a single donor. The heatmap uses:
#' \itemize{
#'   \item Rows: Informative heterozygous SNPs on the X chromosome (one per gene)
#'   \item Columns: Individual cells/barcodes
#'   \item Values: Sign of allelic imbalance (+1 for REF bias, -1 for ALT bias, 0 for balanced)
#'   \item Clustering: Hierarchical clustering with complete linkage and Euclidean distance
#' }
#'
#' The clustering should reveal two main groups of cells, each showing opposite
#' allelic bias patterns across SNPs. This bimodal pattern reflects which of the
#' two parental X chromosomes is inactive in each cell.
#'
#' Only SNPs with at least 2 reads coverage in at least 1% of cells are included
#' to focus on informative loci.
#'
#' @param x SNPData object containing X chromosome SNP data
#' @param donor Character string specifying which donor to visualize
#'
#' @return A ComplexHeatmap object showing the inactive X clustering pattern
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Plot inactive X pattern for a specific donor
#' plot_inactive_x_heatmap(snp_data, donor = "donor1")
#'
#' # Save the plot
#' pdf("inactive_x_heatmap.pdf", width = 10, height = 8)
#' plot_inactive_x_heatmap(snp_data, donor = "donor1")
#' dev.off()
#' }
setGeneric("plot_inactive_x_heatmap", function(x, donor) standardGeneric("plot_inactive_x_heatmap"))

#' @rdname plot_inactive_x_heatmap
#' @include SNPData-class.R
setMethod("plot_inactive_x_heatmap", signature(x = "SNPData"), function(x, donor) {
    expr_matrix <- .prepare_expr_matrix(
        filter_barcodes(x, donor == donor),
        min_coverage = 2,
        min_sample_prop = 0.01
    )

    ComplexHeatmap::Heatmap(
        as.matrix(expr_matrix),
        clustering_method_columns = "complete",
        clustering_distance_columns = "ward.D2",
        show_row_names = FALSE,
        show_column_names = FALSE,
        name = paste("Donor", donor, "Inactive X Heatmap")
    )
})

#' Plot inactive X chromosome heatmap for a donor at clonotype level
#'
#' Visualizes allelic imbalance patterns at informative heterozygous SNPs to
#' reveal X-inactivation structure at the clonotype level. Clonotypes are
#' clustered based on their allelic bias signatures, which should separate into
#' two groups corresponding to the two inactive X states.
#'
#' @details
#' This function creates a heatmap showing the sign of allelic expression at
#' heterozygous SNPs across clonotypes from a single donor. The heatmap uses:
#' \itemize{
#'   \item Rows: Informative heterozygous SNPs on the X chromosome (one per gene)
#'   \item Columns: Clonotypes (aggregated from individual cells)
#'   \item Values: Sign of allelic imbalance (+1 for REF bias, -1 for ALT bias, 0 for balanced)
#'   \item Clustering: Hierarchical clustering with complete linkage and Euclidean distance
#' }
#'
#' The clustering should reveal two main groups of clonotypes, each showing opposite
#' allelic bias patterns across SNPs. This bimodal pattern reflects which of the
#' two parental X chromosomes is inactive in cells from each clonotype.
#'
#' Only SNPs with at least 2 reads coverage in at least 1% of clonotypes are included
#' to focus on informative loci. Counts are aggregated at the clonotype level before
#' computing allelic imbalance.
#'
#' @param x SNPData object containing X chromosome SNP data with clonotype information
#' @param donor Character string specifying which donor to visualize
#'
#' @return A ComplexHeatmap object showing the inactive X clustering pattern at clonotype level
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Plot inactive X pattern for a specific donor at clonotype level
#' plot_inactive_x_heatmap_by_clonotype(snp_data, donor = "donor1")
#'
#' # Save the plot
#' pdf("inactive_x_heatmap_clonotype.pdf", width = 10, height = 8)
#' plot_inactive_x_heatmap_by_clonotype(snp_data, donor = "donor1")
#' dev.off()
#' }
setGeneric("plot_inactive_x_heatmap_by_clonotype", function(x, donor) standardGeneric("plot_inactive_x_heatmap_by_clonotype"))

#' @rdname plot_inactive_x_heatmap_by_clonotype
#' @include SNPData-class.R
setMethod("plot_inactive_x_heatmap_by_clonotype", signature(x = "SNPData"), function(x, donor) {
    expr_matrix <- .prepare_expr_matrix_by_clonotype(
        filter_barcodes(x, donor == donor),
        min_coverage = 2,
        min_sample_prop = 0.01
    )

    ComplexHeatmap::Heatmap(
        as.matrix(expr_matrix),
        clustering_method_columns = "complete",
        clustering_distance_columns = "ward.D2",
        show_row_names = FALSE,
        show_column_names = FALSE,
        name = paste("Donor", donor, "Inactive X Heatmap (Clonotype)")
    )
})

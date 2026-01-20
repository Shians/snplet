#' Assign inactive x
assign_inactive_x <- function(
    snp_data
) {
    unique_donors <- unique(get_barcode_info(snp_data)$donor)

    unique_donor_snp_data <- map(
        unique_donors,
        ~filter_samples(snp_data, donor == .x)
    ) %>%
        magrittr::set_names(unique_donors)

    # assign inactive x for each donor and combine results
    barcode_inactive_x_df <- unique_donor_snp_data %>%
        map(.assign_inactive_x_single_donor) %>%
        bind_rows()

    # add inactive_x to barcode_info
    snp_data %>%
        add_barcode_metadata(barcode_inactive_x_df)
}

#' Assign inactive x for a single donor
#' @returns data.frame with barcode and inactive_x assignment
.assign_inactive_x_single_donor <- function(
    snp_data
) {
    expr_matrix <- .prepare_expr_matrix(
        snp_data,
        min_coverage = 2,
        min_sample_prop = 0.01
    )

    logger::log_success("Inactive X assignment complete for donor {donor}")
    .assign_cluster(expr_matrix, n_clusters = 2)
}

.get_informative_snps <- function(coverage_mat, min_coverage, min_sample_prop) {
    barcodes_with_sufficient_coverage <- rowSums(coverage_mat >= min_coverage)
    barcode_threshold <- ncol(coverage_mat) * min_sample_prop
    barcodes_with_sufficient_coverage >= barcode_threshold
}

.assign_cluster <- function(expr_mat, n_clusters) {
    dist_mat <- dist(t(as.matrix(sign(expr_mat))), method = "euclidean")
    hc <- hclust(dist_mat, method = "complete")
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

plot_inactive_x_heatmap <- function(snp_data, donor) {
    expr_matrix <- .prepare_expr_matrix(
        filter_barcodes(snp_data, donor == donor),
        min_coverage = 2,
        min_sample_prop = 0.01
    )

    ComplexHeatmap::Heatmap(
        sign(as.matrix(expr_matrix)),
        clustering_method_columns = "complete",
        clustering_distance_columns = "euclidean",
        show_row_names = FALSE,
        show_column_names = FALSE,
        name = paste("Donor", donor, "Inactive X Heatmap")
    )
}

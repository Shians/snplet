#' Assign Inactive X Chromosome
#'
#' Infers which X chromosome is inactive in female cells based on allele-specific
#' expression patterns. Uses hierarchical clustering on high-coverage chrX SNPs
#' to identify two groups of cells with correlated expression patterns, then
#' fits a Gaussian Mixture Model to classify all cells.
#'
#' @param x A SNPData object
#' @param snp_quantile Quantile threshold for SNP coverage filtering (default 0.9).
#'   Only chrX SNPs with coverage >= this quantile are used for clustering.
#' @param cell_quantile Quantile threshold for cell library size filtering (default 0.8).
#'   Only cells with library size >= this quantile are used for training the model.
#' @param max_snps Maximum number of SNPs to use for clustering (default 200).
#'   If more SNPs pass the quantile filter, the top max_snps by coverage are selected.
#' @param max_cells Maximum number of cells to use for training (default 1000).
#'   If more cells pass the quantile filter, max_cells are randomly sampled.
#' @param min_total_count Minimum total read depth per donor for heterozygosity testing (default 10).
#'   Passed to donor_het_status_df for identifying heterozygous SNPs.
#' @param het_p_value P-value threshold for heterozygosity testing (default 0.05).
#'   Passed to donor_het_status_df for identifying heterozygous SNPs.
#'
#' @return A SNPData object with added barcode metadata columns:
#'   \itemize{
#'     \item inactive_x: Predicted inactive X chromosome ("X1" or "X2")
#'     \item inactive_x_prob: Probability of the assignment (max of the two class probabilities)
#'   }
#'
#' @details
#' The function follows these steps:
#' \enumerate{
#'   \item Filter to chrX SNPs with high coverage (>= snp_quantile)
#'   \item Filter to heterozygous SNPs using donor_het_status_df
#'   \item Subsample to at most max_snps SNPs (highest coverage)
#'   \item Filter to cells with high library size (>= cell_quantile)
#'   \item Subsample to at most max_cells cells (random sampling)
#'   \item Convert to expression matrix using sign(REF-ALT) * log1p(|REF-ALT|)
#'   \item Cluster cells using hierarchical clustering (Ward's method) on correlation distance
#'   \item Cut dendrogram into 2 clusters
#'   \item Calculate cluster assignment scores for all cells
#'   \item Fit Gaussian Mixture Model on training scores
#'   \item Predict inactive X for all cells using the GMM
#' }
#'
#' Note: This function requires donor information in barcode metadata, as it filters
#' to heterozygous SNPs based on donor-level genotypes. Homozygous SNPs cannot be
#' used for X inactivation analysis since they lack allele-specific signal.
#'
#' @examples
#' \dontrun{
#' # Assign inactive X chromosome
#' snpdata <- assign_inactive_x(snpdata)
#'
#' # View results
#' barcode_info <- get_barcode_info(snpdata)
#' table(barcode_info$inactive_x)
#'
#' # Use more stringent filtering
#' snpdata <- assign_inactive_x(
#'     snpdata,
#'     snp_quantile = 0.95,
#'     cell_quantile = 0.9
#' )
#'
#' # Use fewer SNPs and cells for faster computation
#' snpdata <- assign_inactive_x(
#'     snpdata,
#'     max_snps = 100,
#'     max_cells = 500
#' )
#'
#' # Adjust heterozygosity testing parameters
#' snpdata <- assign_inactive_x(
#'     snpdata,
#'     min_total_count = 20,
#'     het_p_value = 0.01
#' )
#' }
#'
#' @export
assign_inactive_x <- function(
    x,
    snp_quantile = 0.9,
    cell_quantile = 0.8,
    max_snps = 200,
    max_cells = 1000,
    min_total_count = 10,
    het_p_value = 0.05
) {
    # Validate input object type
    if (!methods::is(x, "SNPData")) {
        stop("assign_inactive_x expects a SNPData object")
    }

    # Validate that chromosome style is known
    .validate_chr_style(x, "assign_inactive_x")

    # Validate snp_quantile parameter range
    if (snp_quantile < 0 || snp_quantile > 1) {
        stop("snp_quantile must be between 0 and 1")
    }

    # Validate cell_quantile parameter range
    if (cell_quantile < 0 || cell_quantile > 1) {
        stop("cell_quantile must be between 0 and 1")
    }

    # Validate max_snps parameter
    if (!is.numeric(max_snps) || max_snps < 2 || max_snps != as.integer(max_snps)) {
        stop("max_snps must be an integer >= 2")
    }

    # Validate max_cells parameter
    if (!is.numeric(max_cells) || max_cells < 3 || max_cells != as.integer(max_cells)) {
        stop("max_cells must be an integer >= 3")
    }

    # Validate min_total_count parameter
    if (!is.numeric(min_total_count) || min_total_count < 1) {
        stop("min_total_count must be a numeric value >= 1")
    }

    # Validate het_p_value parameter
    if (!is.numeric(het_p_value) || het_p_value <= 0 || het_p_value > 1) {
        stop("het_p_value must be between 0 and 1")
    }

    # Check for required chromosome annotation
    snp_info <- get_snp_info(x)
    if (!"chrom_canonical" %in% colnames(snp_info)) {
        stop("snp_info must contain a 'chrom_canonical' column. This should be automatically added during SNPData creation.")
    }

    # Verify presence of chrX SNPs in the dataset
    chrx_snps <- sum(snp_info$chrom_canonical == "chrX", na.rm = TRUE)
    if (chrx_snps == 0) {
        stop("No chrX SNPs found in dataset")
    }

    # Check for donor information (required for heterozygosity testing)
    barcode_info <- get_barcode_info(x)
    if (!"donor" %in% colnames(barcode_info)) {
        stop("assign_inactive_x requires donor information in barcode metadata. Please add donor assignments using Vireo or similar tools.")
    }

    # Filter to high-coverage chrX SNPs for clustering
    # These SNPs provide the most reliable signal for X inactivation patterns
    logger::with_log_threshold({
        snp_subset <- x %>%
            filter_snps(chrom_canonical == "chrX") %>%
            filter_snps(coverage >= stats::quantile(coverage, snp_quantile))
        },
        threshold = logger::INFO
    )

    # Filter to heterozygous SNPs (required for X inactivation analysis)
    het_status <- donor_het_status_df(
        snp_subset,
        min_total_count = min_total_count,
        p_value_threshold = het_p_value,
        minor_allele_prop = 0.1
    )

    # Identify SNPs that are heterozygous in at least one donor
    het_snp_ids <- het_status %>%
        dplyr::filter(zygosity == "het") %>%
        dplyr::pull(snp_id) %>%
        unique()

    if (length(het_snp_ids) == 0) {
        stop("No heterozygous chrX SNPs found. Try adjusting het_p_value or min_total_count.")
    }

    logger::with_log_threshold({
        snp_subset <- snp_subset %>%
            filter_snps(snp_id %in% het_snp_ids)
        },
        threshold = logger::INFO
    )

    # Subsample to at most max_snps SNPs (select highest coverage SNPs)
    snp_info_subset <- get_snp_info(snp_subset)
    if (nrow(snp_info_subset) > max_snps) {
        top_snp_ids <- snp_info_subset %>%
            dplyr::arrange(dplyr::desc(coverage)) %>%
            dplyr::slice_head(n = max_snps) %>%
            dplyr::pull(snp_id)

        logger::with_log_threshold({
            snp_subset <- snp_subset %>%
                filter_snps(snp_id %in% top_snp_ids)
            },
            threshold = logger::INFO
        )
    }

    # Filter to high library size cells for training
    # These cells have sufficient reads to reliably estimate expression patterns
    logger::with_log_threshold({
        high_cells <- snp_subset %>%
            filter_barcodes(library_size >= stats::quantile(library_size, cell_quantile))
        },
        threshold = logger::INFO
    )

    # Subsample to at most max_cells cells (random sampling)
    barcode_info_subset <- get_barcode_info(high_cells)
    if (nrow(barcode_info_subset) > max_cells) {
        set.seed(42)
        sampled_cell_ids <- barcode_info_subset %>%
            dplyr::slice_sample(n = max_cells) %>%
            dplyr::pull(cell_id)

        logger::with_log_threshold({
            high_cells <- high_cells %>%
                filter_barcodes(cell_id %in% sampled_cell_ids)
            },
            threshold = logger::INFO
        )
    }

    # Convert ALT/REF counts to signed expression values
    # sign(REF-ALT) * log1p(|REF-ALT|) captures both direction and magnitude
    expr_high <- snplet::to_expr_matrix(high_cells) %>%
        as("matrix")

    # Verify sufficient SNPs remain after filtering
    if (nrow(expr_high) < 2) {
        stop("At least 2 SNPs required after filtering. Try lower snp_quantile.")
    }

    # Verify sufficient cells remain after filtering
    if (ncol(expr_high) < 3) {
        stop("At least 3 cells required after filtering. Try lower cell_quantile.")
    }

    # Compute pairwise correlation between cells across chrX SNPs
    # High correlation indicates similar X inactivation patterns
    cor_mat <- tryCatch(
        stats::cor(expr_high, method = "pearson"),
        error = function(e) {
            stop(paste0(
                "Failed to compute correlation matrix. ",
                "This may be due to constant or near-constant SNP values. ",
                "Original error: ",
                e$message
            ))
        }
    )

    # Validate correlation matrix quality
    if (any(is.na(cor_mat)) || any(is.nan(cor_mat))) {
        stop("Correlation matrix contains NA or NaN values")
    }

    # Convert correlation to distance and perform hierarchical clustering
    # Ward's method minimizes within-cluster variance
    dist_mat <- stats::as.dist(1 - cor_mat)
    hc <- stats::hclust(dist_mat, method = "ward.D2")
    cluster_id <- stats::cutree(hc, k = 2)

    # Verify that clustering produced two distinct groups
    if (length(unique(cluster_id)) != 2) {
        stop("Clustering failed to produce 2 distinct clusters")
    }

    # Compute cluster centroids by averaging expression within each cluster
    # Sign normalization creates two opposite reference vectors
    assignment_score_matrix <- DelayedArray::colsum(expr_high, group = cluster_id) %>%
        as.matrix() %>%
        sign()

    # Project training cells onto cluster centroids to obtain scores
    training_scores <- t(assignment_score_matrix) %*% expr_high

    # Fit Gaussian Mixture Model on training scores to learn cluster distributions
    # This allows probabilistic assignment of all cells including low-coverage ones
    gmm <- tryCatch(
        mclust::Mclust(training_scores[1, ], G = 2, verbose = FALSE),
        error = function(e) {
            stop(paste0(
                "GMM fitting failed. ",
                "This may occur with insufficient data separation. ",
                "Original error: ",
                e$message
            ))
        }
    )

    # Verify GMM convergence
    if (is.null(gmm)) {
        stop("GMM fitting failed to converge")
    }

    # Filter full dataset to SNPs used for scoring
    logger::with_log_threshold({
            scored_data <- x %>%
                filter_snps(snp_id %in% rownames(assignment_score_matrix))
        },
        threshold = logger::INFO
    )

    # Convert all cells to expression values using the same transformation
    expr_all <- snplet::to_expr_matrix(scored_data) %>%
        as("matrix")

    # Project all cells onto cluster centroids
    all_scores <- t(assignment_score_matrix) %*% expr_all

    # Use trained GMM to predict inactive X for all cells
    preds <- mclust::predict.Mclust(gmm, all_scores[1, ])

    # Package predictions as metadata
    # inactive_x: predicted inactive chromosome (X1 or X2)
    # inactive_x_prob: confidence of assignment (max posterior probability)
    metadata <- tibble::tibble(
        cell_id = colnames(expr_all),
        inactive_x = paste0("X", preds$classification),
        inactive_x_prob = pmax(preds$z[, 1], preds$z[, 2])
    )

    # Add predictions to SNPData object and return
    snplet::add_barcode_metadata(
        x,
        metadata
    )
}

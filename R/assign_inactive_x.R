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
#'   \item Filter to cells with high library size (>= cell_quantile)
#'   \item Convert to expression matrix using sign(REF-ALT) * log1p(|REF-ALT|)
#'   \item Cluster cells using hierarchical clustering (Ward's method) on correlation distance
#'   \item Cut dendrogram into 2 clusters
#'   \item Calculate cluster assignment scores for all cells
#'   \item Fit Gaussian Mixture Model on training scores
#'   \item Predict inactive X for all cells using the GMM
#' }
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
#' }
#'
#' @export
assign_inactive_x <- function(
    x,
    snp_quantile = 0.9,
    cell_quantile = 0.8
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

    # Filter to high-coverage chrX SNPs for clustering
    # These SNPs provide the most reliable signal for X inactivation patterns
    logger::with_log_threshold({
        snp_subset <- x %>%
            filter_snps(chrom_canonical == "chrX") %>%
            filter_snps(coverage >= stats::quantile(coverage, snp_quantile))
        },
        threshold = logger::INFO
    )

    # Filter to high library size cells for training
    # These cells have sufficient reads to reliably estimate expression patterns
    logger::with_log_threshold({
        high_cells <- snp_subset %>%
            filter_barcodes(library_size >= stats::quantile(library_size, cell_quantile))
        },
        threshold = logger::INFO
    )

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

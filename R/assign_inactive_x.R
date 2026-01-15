#' Assign Inactive X Chromosome
#'
#' Infers which X chromosome is inactive in female cells based on allele-specific
#' expression patterns. X inactivation is determined independently for each donor,
#' as each female individual has her own pattern of X inactivation. Uses hierarchical
#' clustering on high-coverage heterozygous chrX SNPs to identify two groups of cells
#' with correlated expression patterns, then fits a Gaussian Mixture Model to classify
#' all cells from that donor.
#'
#' @param x A SNPData object with donor information in barcode metadata
#' @param min_cell_prop Minimum proportion of cells with non-zero counts for a SNP
#'   to be included in clustering (default 0.2). After filtering cells by cell_quantile,
#'   only SNPs with non-zero counts in at least this proportion of remaining cells are used.
#'   Applied per donor.
#' @param cell_quantile Quantile threshold for cell library size filtering (default 0.8).
#'   Only cells with library size >= this quantile are used for training the model.
#'   Applied per donor.
#' @param max_snps Maximum number of SNPs to use for clustering (default 200).
#'   If more SNPs pass the quantile filter, the top max_snps by coverage are selected.
#'   Applied per donor.
#' @param max_cells Maximum number of cells (or clonotypes if level = "clonotype")
#'   to use for training (default 1000). If more pass the quantile filter, the top
#'   max_cells by library size are selected. Applied per donor.
#' @param min_total_count Minimum total read depth per donor for heterozygosity testing (default 10).
#'   Passed to donor_het_status_df for identifying heterozygous SNPs.
#' @param het_p_value P-value threshold for heterozygosity testing (default 0.05).
#'   Passed to donor_het_status_df for identifying heterozygous SNPs.
#' @param level Either "barcode" (default) or "clonotype". When "clonotype",
#'   counts are aggregated per clonotype before assignment, and clonotype-level
#'   assignments are then assigned back to individual cells.
#'   Requires a "clonotype" column in barcode metadata.
#'
#' @return A SNPData object with added barcode metadata columns:
#'   \itemize{
#'     \item inactive_x: Predicted inactive X chromosome ("X1" or "X2")
#'     \item inactive_x_prob: Probability of the assignment (max of the two class probabilities)
#'   }
#'
#' @details
#' The function first removes doublets and unassigned cells, then processes each valid
#' donor independently with the following steps:
#' \enumerate{
#'   \item Remove cells with donor = "doublet" or "unassigned"
#'   \item For each valid donor:
#'   \itemize{
#'     \item When level = "barcode":
#'       \itemize{
#'         \item Filter to cells with high library size (>= cell_quantile)
#'         \item Select top max_cells cells by library size
#'       }
#'     \item When level = "clonotype":
#'       \itemize{
#'         \item Filter to clonotypes with high total library size (>= cell_quantile)
#'         \item Select top max_cells clonotypes by total library size
#'       }
#'     \item Identify heterozygous chrX SNPs for this donor using donor_het_status_df
#'     \item Filter to SNPs with non-zero counts in at least min_cell_prop of the filtered cells
#'     \item Subsample to at most max_snps SNPs (highest coverage)
#'     \item Convert to expression matrix using sign(REF-ALT) * log1p(|REF-ALT|)
#'     \item Cluster using hierarchical clustering (Ward's method) on correlation distance
#'     \item Cut dendrogram into 2 clusters
#'     \item Calculate cluster assignment scores
#'     \item Fit Gaussian Mixture Model on training scores
#'     \item Predict inactive X using the GMM
#'     \item When level = "clonotype": Map clonotype predictions back to individual cells
#'   }
#' }
#'
#' Note: This function requires donor information in barcode metadata, as X inactivation
#' is donor-specific and heterozygous SNPs are identified per donor. Doublets and
#' unassigned cells are automatically excluded from analysis. Homozygous SNPs cannot
#' be used for X inactivation analysis since they lack allele-specific signal.
#'
#' @examples
#' \dontrun{
#' # Assign inactive X chromosome (analyzed per donor)
#' snpdata <- assign_inactive_x(snpdata)
#'
#' # View results by donor
#' barcode_info <- get_barcode_info(snpdata)
#' table(barcode_info$donor, barcode_info$inactive_x)
#'
#' # Use more stringent filtering
#' snpdata <- assign_inactive_x(
#'     snpdata,
#'     min_cell_prop = 0.3,
#'     cell_quantile = 0.9
#' )
#'
#' # Assign by clonotype instead of individual cells
#' snpdata <- assign_inactive_x(
#'     snpdata,
#'     level = "clonotype"
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
    min_cell_prop = 0.2,
    cell_quantile = 0.8,
    max_snps = 200,
    max_cells = 1000,
    min_total_count = 10,
    het_p_value = 0.05,
    level = c("barcode", "clonotype")
) {
    # Validate input object type
    if (!methods::is(x, "SNPData")) {
        stop("assign_inactive_x expects a SNPData object")
    }

    # Validate level parameter
    level <- match.arg(level)

    # Validate that chromosome style is known
    .validate_chr_style(x, "assign_inactive_x")

    # Validate min_cell_prop parameter range
    if (min_cell_prop < 0 || min_cell_prop > 1) {
        stop("min_cell_prop must be between 0 and 1")
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
        stop(
            "snp_info must contain a 'chrom_canonical' column. This should be automatically added during SNPData creation."
        )
    }

    # Verify presence of chrX SNPs in the dataset
    chrx_snps <- sum(snp_info$chrom_canonical == "chrX", na.rm = TRUE)
    if (chrx_snps == 0) {
        stop("No chrX SNPs found in dataset")
    }

    # Check for donor information (required for heterozygosity testing and per-donor analysis)
    barcode_info <- get_barcode_info(x)
    if (!"donor" %in% colnames(barcode_info)) {
        stop(
            "assign_inactive_x requires donor information in barcode metadata. Please add donor assignments using Vireo or similar tools."
        )
    }

    # Check for clonotype information if aggregating by clonotype
    if (level == "clonotype" && !"clonotype" %in% colnames(barcode_info)) {
        stop("level = 'clonotype' requires a 'clonotype' column in barcode metadata")
    }

    # Remove doublets and unassigned cells before processing
    n_cells_before <- nrow(barcode_info)
    logger::with_log_threshold(
        {
            x <- x %>%
                filter_barcodes(!donor %in% c("doublet", "unassigned"))
        },
        threshold = logger::INFO
    )

    barcode_info <- get_barcode_info(x)
    n_cells_after <- nrow(barcode_info)
    n_removed <- n_cells_before - n_cells_after

    if (n_removed > 0) {
        logger::log_info(
            "Removed {n_removed} doublet/unassigned cells ({round(100 * n_removed / n_cells_before, 1)}% of total)"
        )
    }

    # Get list of valid donors (excluding doublets and unassigned)
    donors <- unique(barcode_info$donor)
    donors <- donors[!is.na(donors)]
    donors <- donors[!donors %in% c("doublet", "unassigned")]

    if (length(donors) == 0) {
        stop("No valid donors found in barcode metadata after removing doublets and unassigned cells")
    }

    logger::log_info("Processing {length(donors)} donor(s) for X inactivation assignment")

    # Process each donor independently
    all_metadata <- list()
    for (donor_id in donors) {
        logger::log_info("Processing donor: {donor_id}")

        # Process this donor
        donor_metadata <- tryCatch(
            .assign_inactive_x_single_donor(
                x = x,
                donor_id = donor_id,
                min_cell_prop = min_cell_prop,
                cell_quantile = cell_quantile,
                max_snps = max_snps,
                max_cells = max_cells,
                min_total_count = min_total_count,
                het_p_value = het_p_value,
                level = level
            ),
            error = function(e) {
                logger::log_warn("Failed to process donor {donor_id}: {e$message}")
                NULL
            }
        )

        if (!is.null(donor_metadata)) {
            all_metadata[[donor_id]] <- donor_metadata
        }
    }

    if (length(all_metadata) == 0) {
        stop("Failed to process any donors. Check logs for details.")
    }

    # Combine metadata from all donors
    combined_metadata <- dplyr::bind_rows(all_metadata)

    # Add predictions to SNPData object and return
    snplet::add_barcode_metadata(x, combined_metadata)
}

#' Remove constant rows and columns from a matrix
#'
#' Removes rows and columns with zero variance from a numeric matrix.
#' This is useful for preparing data for correlation analysis or other
#' statistical methods that fail with constant values.
#'
#' @param mat A numeric matrix
#'
#' @return A list with components:
#'   \itemize{
#'     \item matrix: The filtered matrix with constant rows/columns removed
#'     \item n_constant_rows: Number of rows removed
#'     \item n_constant_cols: Number of columns removed
#'   }
#'
#' @keywords internal
.remove_constant_rows_and_columns <- function(mat) {
    # Calculate variance for rows and columns
    row_variance <- apply(mat, 1, stats::var)
    col_variance <- apply(mat, 2, stats::var)

    # Identify non-constant rows and columns
    non_constant_rows <- row_variance > 0 & !is.na(row_variance)
    non_constant_cols <- col_variance > 0 & !is.na(col_variance)

    # Count how many were removed
    n_constant_rows <- sum(!non_constant_rows)
    n_constant_cols <- sum(!non_constant_cols)

    # Filter matrix
    filtered_mat <- mat[non_constant_rows, non_constant_cols, drop = FALSE]

    list(
        matrix = filtered_mat,
        n_constant_rows = n_constant_rows,
        n_constant_cols = n_constant_cols
    )
}

#' Filter SNPs to heterozygous chrX SNPs for a donor
#'
#' @param x A SNPData object
#' @param donor_id The donor ID to process
#' @param min_total_count Minimum total count for heterozygosity testing
#' @param het_p_value P-value threshold for heterozygosity testing
#' @param min_cell_prop Minimum proportion of cells with non-zero counts
#' @param max_snps Maximum number of SNPs to use
#' @param filtered_cells SNPData object already filtered to the desired cells
#'
#' @return A SNPData object filtered to heterozygous chrX SNPs
#' @keywords internal
.filter_to_het_chrx_snps <- function(
    x,
    donor_id,
    min_total_count,
    het_p_value,
    min_cell_prop,
    max_snps,
    filtered_cells
) {
    logger::with_log_threshold(
        {
            donor_data <- x %>%
                filter_barcodes(donor == donor_id)
        },
        threshold = logger::INFO
    )

    logger::with_log_threshold(
        {
            chrx_data <- donor_data %>%
                filter_snps(chrom_canonical == "chrX")
        },
        threshold = logger::INFO
    )

    het_status <- donor_het_status_df(
        chrx_data,
        min_total_count = min_total_count,
        p_value_threshold = het_p_value,
        minor_allele_prop = 0.1
    )

    het_snp_ids <- het_status %>%
        dplyr::filter(donor == donor_id, zygosity == "het") %>%
        dplyr::pull(snp_id) %>%
        unique()

    if (length(het_snp_ids) == 0) {
        stop(glue::glue("No heterozygous chrX SNPs found for donor {donor_id}"))
    }

    logger::log_info("Found {length(het_snp_ids)} heterozygous SNPs for donor {donor_id}")

    logger::with_log_threshold(
        {
            het_data <- chrx_data %>%
                filter_snps(snp_id %in% het_snp_ids)
        },
        threshold = logger::INFO
    )

    # Filter to the same cells that will be used for clustering
    logger::with_log_threshold(
        {
            het_data_filtered_cells <- het_data %>%
                filter_barcodes(cell_id %in% get_barcode_info(filtered_cells)$cell_id)
        },
        threshold = logger::INFO
    )

    # Calculate total coverage for each SNP in the filtered cells
    total_coverage <- alt_count(het_data_filtered_cells) + ref_count(het_data_filtered_cells)

    # Count cells with non-zero coverage for each SNP
    cells_with_coverage <- Matrix::rowSums(total_coverage > 0)
    n_cells <- ncol(total_coverage)

    # Filter to SNPs with non-zero counts in at least min_cell_prop of cells
    prop_cells_with_coverage <- cells_with_coverage / n_cells
    passing_snp_indices <- which(prop_cells_with_coverage >= min_cell_prop)

    if (length(passing_snp_indices) == 0) {
        stop(glue::glue(
            "No SNPs found with non-zero counts in at least {min_cell_prop * 100}% of cells for donor {donor_id}"
        ))
    }

    passing_snp_ids <- rownames(total_coverage)[passing_snp_indices]

    logger::log_info("Found {length(passing_snp_ids)} SNPs with non-zero counts in >= {min_cell_prop * 100}% of cells")

    logger::with_log_threshold(
        {
            snp_subset <- het_data %>%
                filter_snps(snp_id %in% passing_snp_ids)
        },
        threshold = logger::INFO
    )

    # If we have more SNPs than max_snps, select top by coverage
    snp_info_subset <- get_snp_info(snp_subset)
    if (nrow(snp_info_subset) > max_snps) {
        top_snp_ids <- snp_info_subset %>%
            dplyr::arrange(dplyr::desc(coverage)) %>%
            dplyr::slice_head(n = max_snps) %>%
            dplyr::pull(snp_id)

        logger::with_log_threshold(
            {
                snp_subset <- snp_subset %>%
                    filter_snps(snp_id %in% top_snp_ids)
            },
            threshold = logger::INFO
        )
    }

    snp_subset
}

#' Cluster cells/clonotypes and fit GMM
#'
#' @param expr_mat Expression matrix (SNPs x cells/clonotypes)
#' @param donor_id The donor ID (for error messages)
#'
#' @return List with cluster assignments and fitted GMM model
#' @keywords internal
.cluster_and_fit_gmm <- function(expr_mat, donor_id) {
    filtered_result <- .remove_constant_rows_and_columns(expr_mat)

    if (filtered_result$n_constant_rows > 0) {
        logger::log_info("Donor {donor_id}: Removing {filtered_result$n_constant_rows} constant SNP(s)")
    }
    if (filtered_result$n_constant_cols > 0) {
        logger::log_info("Donor {donor_id}: Removing {filtered_result$n_constant_cols} constant cell(s)")
    }

    expr_mat <- filtered_result$matrix

    if (nrow(expr_mat) < 2) {
        stop(glue::glue("Donor {donor_id}: At least 2 non-constant SNPs required. Try lower min_cell_prop."))
    }
    if (ncol(expr_mat) < 3) {
        stop(glue::glue("Donor {donor_id}: At least 3 non-constant cells required. Try lower cell_quantile."))
    }

    cor_mat <- tryCatch(
        stats::cor(expr_mat, method = "pearson"),
        error = function(e) {
            stop(glue::glue(
                "Donor {donor_id}: Failed to compute correlation matrix. ",
                "Original error: {e$message}"
            ))
        }
    )

    if (any(is.na(cor_mat)) || any(is.nan(cor_mat))) {
        stop(glue::glue("Donor {donor_id}: Correlation matrix contains NA or NaN values"))
    }

    dist_mat <- stats::as.dist(1 - cor_mat)
    hc <- stats::hclust(dist_mat, method = "ward.D2")
    cluster_id <- stats::cutree(hc, k = 2)

    if (length(unique(cluster_id)) != 2) {
        stop(glue::glue("Donor {donor_id}: Clustering failed to produce 2 distinct clusters"))
    }

    assignment_score_matrix <- DelayedArray::colsum(expr_mat, group = cluster_id) %>%
        as.matrix() %>%
        sign()

    training_scores <- t(assignment_score_matrix) %*% expr_mat

    gmm <- tryCatch(
        mclust::Mclust(training_scores[1, ], G = 2, verbose = FALSE),
        error = function(e) {
            stop(glue::glue(
                "Donor {donor_id}: GMM fitting failed. ",
                "This may occur with insufficient data separation. ",
                "Original error: {e$message}"
            ))
        }
    )

    if (is.null(gmm)) {
        stop(glue::glue("Donor {donor_id}: GMM fitting failed to converge"))
    }

    list(
        assignment_score_matrix = assignment_score_matrix,
        gmm = gmm
    )
}

#' Predict inactive X for all cells (clonotype aggregation)
#'
#' @param x A SNPData object
#' @param donor_id The donor ID to process
#' @param assignment_score_matrix Matrix for scoring cells
#' @param gmm Fitted GMM model
#'
#' @return Tibble with cell_id, inactive_x, and inactive_x_prob
#' @keywords internal
.predict_clonotype_assignments <- function(
    x,
    donor_id,
    assignment_score_matrix,
    gmm
) {
    logger::with_log_threshold(
        {
            scored_data <- x %>%
                filter_barcodes(donor == donor_id) %>%
                filter_snps(snp_id %in% rownames(assignment_score_matrix))
        },
        threshold = logger::INFO
    )

    # Get clonotype-level expression matrix
    expr_all <- snplet::to_expr_matrix(scored_data, level = "clonotype") %>%
        as("matrix")

    all_scores <- t(assignment_score_matrix) %*% expr_all
    preds <- mclust::predict.Mclust(gmm, all_scores[1, ])

    clonotype_metadata <- tibble::tibble(
        clonotype = colnames(expr_all),
        inactive_x = paste0("X", preds$classification),
        inactive_x_prob = pmax(preds$z[, 1], preds$z[, 2])
    )

    get_barcode_info(scored_data) %>%
        dplyr::select(cell_id, clonotype) %>%
        dplyr::left_join(clonotype_metadata, by = "clonotype") %>%
        dplyr::select(cell_id, inactive_x, inactive_x_prob)
}

#' Predict inactive X for all cells (cell-level)
#'
#' @param x A SNPData object
#' @param donor_id The donor ID to process
#' @param assignment_score_matrix Matrix for scoring cells
#' @param gmm Fitted GMM model
#'
#' @return Tibble with cell_id, inactive_x, and inactive_x_prob
#' @keywords internal
.predict_cell_assignments <- function(
    x,
    donor_id,
    assignment_score_matrix,
    gmm
) {
    logger::with_log_threshold(
        {
            scored_data <- x %>%
                filter_barcodes(donor == donor_id) %>%
                filter_snps(snp_id %in% rownames(assignment_score_matrix))
        },
        threshold = logger::INFO
    )

    expr_all <- snplet::to_expr_matrix(scored_data) %>%
        as("matrix")

    all_scores <- t(assignment_score_matrix) %*% expr_all
    preds <- mclust::predict.Mclust(gmm, all_scores[1, ])

    tibble::tibble(
        cell_id = colnames(expr_all),
        inactive_x = paste0("X", preds$classification),
        inactive_x_prob = pmax(preds$z[, 1], preds$z[, 2])
    )
}

#' Filter cells for clustering based on library size
#'
#' @param donor_data SNPData filtered to a single donor
#' @param cell_quantile Quantile threshold for library size filtering
#' @param max_cells Maximum number of cells/clonotypes to retain
#' @param level Either "barcode" or "clonotype"
#'
#' @return SNPData object with filtered cells
#' @keywords internal
.filter_cells_for_clustering <- function(
    donor_data,
    cell_quantile,
    max_cells,
    level
) {
    if (level == "clonotype") {
        barcode_info <- get_barcode_info(donor_data)

        clonotype_lib_size <- barcode_info %>%
            dplyr::group_by(clonotype) %>%
            dplyr::summarise(total_library_size = sum(library_size, na.rm = TRUE), .groups = "drop")

        lib_size_threshold <- stats::quantile(clonotype_lib_size$total_library_size, cell_quantile)

        high_clonotypes <- clonotype_lib_size %>%
            dplyr::filter(total_library_size >= lib_size_threshold) %>%
            dplyr::pull(clonotype)

        if (length(high_clonotypes) > max_cells) {
            high_clonotypes <- clonotype_lib_size %>%
                dplyr::filter(clonotype %in% high_clonotypes) %>%
                dplyr::arrange(dplyr::desc(total_library_size)) %>%
                dplyr::slice_head(n = max_cells) %>%
                dplyr::pull(clonotype)
        }

        logger::with_log_threshold(
            {
                filtered_cells <- donor_data %>%
                    filter_barcodes(clonotype %in% high_clonotypes)
            },
            threshold = logger::INFO
        )
    } else {
        barcode_lib_size <- get_barcode_info(donor_data)$library_size
        lib_size_threshold <- stats::quantile(barcode_lib_size, cell_quantile)

        logger::with_log_threshold(
            {
                filtered_cells <- donor_data %>%
                    filter_barcodes(library_size >= lib_size_threshold)
            },
            threshold = logger::INFO
        )

        barcode_info <- get_barcode_info(filtered_cells)
        if (nrow(barcode_info) > max_cells) {
            top_cell_ids <- barcode_info %>%
                dplyr::arrange(dplyr::desc(library_size)) %>%
                dplyr::slice_head(n = max_cells) %>%
                dplyr::pull(cell_id)

            logger::with_log_threshold(
                {
                    filtered_cells <- filtered_cells %>%
                        filter_barcodes(cell_id %in% top_cell_ids)
                },
                threshold = logger::INFO
            )
        }
    }

    filtered_cells
}

#' Prepare expression matrix for clustering
#'
#' @param snp_subset SNPData with selected SNPs
#' @param filtered_cells SNPData with filtered cells
#' @param level Either "barcode" or "clonotype"
#'
#' @return Expression matrix (SNPs x cells/clonotypes)
#' @keywords internal
.prepare_expr_matrix_for_clustering <- function(
    snp_subset,
    filtered_cells,
    level
) {
    logger::with_log_threshold(
        {
            snp_subset_filtered <- snp_subset %>%
                filter_barcodes(cell_id %in% get_barcode_info(filtered_cells)$cell_id)
        },
        threshold = logger::INFO
    )

    # Use to_expr_matrix with the specified level
    expr_matrix <- snplet::to_expr_matrix(snp_subset_filtered, level = level) %>%
        as("matrix")

    expr_matrix
}

#' Assign inactive X chromosome for a single donor
#'
#' @param x A SNPData object
#' @param donor_id The donor ID to process
#' @param min_cell_prop Minimum proportion of cells with non-zero counts
#' @param cell_quantile Quantile threshold for cell library size filtering
#' @param max_snps Maximum number of SNPs to use
#' @param max_cells Maximum number of cells to use
#' @param min_total_count Minimum total count for heterozygosity testing
#' @param het_p_value P-value threshold for heterozygosity testing
#' @param level Either "barcode" or "clonotype" for aggregation level
#'
#' @return A tibble with cell_id, inactive_x, and inactive_x_prob columns
#' @keywords internal
.assign_inactive_x_single_donor <- function(
    x,
    donor_id,
    min_cell_prop,
    cell_quantile,
    max_snps,
    max_cells,
    min_total_count,
    het_p_value,
    level
) {
    # Filter to donor
    logger::with_log_threshold(
        {
            donor_data <- x %>%
                filter_barcodes(donor == donor_id)
        },
        threshold = logger::INFO
    )

    # Filter cells based on cell_quantile and max_cells
    filtered_cells <- .filter_cells_for_clustering(
        donor_data = donor_data,
        cell_quantile = cell_quantile,
        max_cells = max_cells,
        level = level
    )

    # Filter SNPs based on the filtered cells
    snp_subset <- .filter_to_het_chrx_snps(
        x = x,
        donor_id = donor_id,
        min_total_count = min_total_count,
        het_p_value = het_p_value,
        min_cell_prop = min_cell_prop,
        max_snps = max_snps,
        filtered_cells = filtered_cells
    )

    # Prepare expression matrix
    expr_high <- .prepare_expr_matrix_for_clustering(
        snp_subset = snp_subset,
        filtered_cells = filtered_cells,
        level = level
    )

    if (nrow(expr_high) < 2) {
        stop(glue::glue("Donor {donor_id}: At least 2 SNPs required after filtering. Try lower min_cell_prop."))
    }

    if (ncol(expr_high) < 3) {
        stop(glue::glue("Donor {donor_id}: At least 3 cells required after filtering. Try lower cell_quantile."))
    }

    if (level == "clonotype") {
        logger::log_info(
            "Donor {donor_id}: Using {nrow(expr_high)} SNPs and {ncol(expr_high)} clonotypes for clustering"
        )
    } else {
        logger::log_info("Donor {donor_id}: Using {nrow(expr_high)} SNPs and {ncol(expr_high)} cells for clustering")
    }

    model <- .cluster_and_fit_gmm(expr_high, donor_id)

    if (level == "clonotype") {
        metadata <- .predict_clonotype_assignments(
            x = x,
            donor_id = donor_id,
            assignment_score_matrix = model$assignment_score_matrix,
            gmm = model$gmm
        )
    } else {
        metadata <- .predict_cell_assignments(
            x = x,
            donor_id = donor_id,
            assignment_score_matrix = model$assignment_score_matrix,
            gmm = model$gmm
        )
    }

    logger::log_success("Donor {donor_id}: Successfully assigned inactive X for {nrow(metadata)} cells")

    metadata
}

# Merge helpers for SNPData

.expand_subset_matrix <- function(mat, retained_rows, retained_cols, col_mapping = NULL) {
    # Use hash-based lookup for faster matching with large vectors
    retained_rows_set <- setNames(seq_along(retained_rows), retained_rows)
    mat_rows <- rownames(mat)

    # Find rows that exist in both mat and retained_rows
    rows_in_both <- mat_rows %in% names(retained_rows_set)
    rows_to_keep <- mat_rows[rows_in_both]
    row_idx_in_mat <- which(rows_in_both)
    row_idx_in_expanded <- retained_rows_set[rows_to_keep]

    if (is.null(col_mapping)) {
        # Original behavior: match by column names
        retained_cols_set <- setNames(seq_along(retained_cols), retained_cols)
        mat_cols <- colnames(mat)

        # Find columns that exist in both
        cols_in_both <- mat_cols %in% names(retained_cols_set)
        cols_to_keep <- mat_cols[cols_in_both]
        col_idx_in_mat <- which(cols_in_both)
        col_idx_in_expanded <- retained_cols_set[cols_to_keep]

        # Subset efficiently
        mat_subset <- mat[row_idx_in_mat, col_idx_in_mat, drop = FALSE]

        # Build sparse matrix using triplet format (much faster)
        if (nrow(mat_subset) > 0 && ncol(mat_subset) > 0) {
            # Convert to triplet format (works for both sparse and dense matrices)
            mat_T <- as(mat_subset, "TsparseMatrix")

            # Remap indices (TsparseMatrix uses 0-based indexing, but slot access gives 1-based)
            new_i <- row_idx_in_expanded[mat_T@i + 1]
            new_j <- col_idx_in_expanded[mat_T@j + 1]

            # Create sparse matrix from triplets
            expanded <- Matrix::sparseMatrix(
                i = new_i,
                j = new_j,
                x = mat_T@x,
                dims = c(length(retained_rows), length(retained_cols))
            )
        } else {
            # Empty matrix case
            expanded <- Matrix::Matrix(
                0,
                nrow = length(retained_rows),
                ncol = length(retained_cols),
                sparse = TRUE
            )
        }

        rownames(expanded) <- retained_rows
        colnames(expanded) <- retained_cols
    } else {
        # Barcode-based behavior: use column mapping
        # Only process columns that map to retained positions
        valid_cols <- !is.na(col_mapping)

        if (any(valid_cols) && length(rows_to_keep) > 0) {
            # Subset rows and valid columns
            col_idx_in_mat <- which(valid_cols)
            col_idx_in_expanded <- col_mapping[valid_cols]

            mat_subset <- mat[row_idx_in_mat, col_idx_in_mat, drop = FALSE]

            # Convert to triplet format for efficient sparse matrix construction
            mat_T <- as(mat_subset, "TsparseMatrix")

            # Remap indices (TsparseMatrix uses 0-based indexing, but slot access gives 1-based)
            new_i <- row_idx_in_expanded[mat_T@i + 1]
            new_j <- col_idx_in_expanded[mat_T@j + 1]

            # Create sparse matrix from triplets
            expanded <- Matrix::sparseMatrix(
                i = new_i,
                j = new_j,
                x = mat_T@x,
                dims = c(length(retained_rows), length(retained_cols))
            )
        } else {
            # Empty matrix case
            expanded <- Matrix::Matrix(
                0,
                nrow = length(retained_rows),
                ncol = length(retained_cols),
                sparse = TRUE
            )
        }

        rownames(expanded) <- retained_rows
        colnames(expanded) <- retained_cols
    }

    return(expanded)
}

.merge_snp_info <- function(x, y, snp_ids_retained, snp_join) {
    join_fun <- switch(
        snp_join,
        "union" = dplyr::full_join,
        "intersect" = dplyr::inner_join,
        "left" = dplyr::left_join,
        "right" = dplyr::right_join
    )

    auto_cols <- c("coverage", "non_zero_samples")
    snp_info_x <- x@snp_info %>% dplyr::select(-dplyr::any_of(auto_cols))
    snp_info_y <- y@snp_info %>% dplyr::select(-dplyr::any_of(auto_cols))

    merged <- join_fun(
        snp_info_x,
        snp_info_y,
        by = "snp_id",
        suffix = c(".x", ".y")
    )

    x_conflicts <- grep("\\.x$", colnames(merged), value = TRUE)
    for (x_col in x_conflicts) {
        base_col <- sub("\\.x$", "", x_col)
        y_col <- paste0(base_col, ".y")
        merged[[base_col]] <- dplyr::coalesce(merged[[x_col]], merged[[y_col]])
        merged[[x_col]] <- NULL
        merged[[y_col]] <- NULL
    }

    merged <- merged[merged$snp_id %in% snp_ids_retained, , drop = FALSE]
    order_idx <- match(snp_ids_retained, merged$snp_id)
    if (any(is.na(order_idx))) {
        stop("Some snp_ids could not be matched during merge_snp_info()")
    }
    merged <- merged[order_idx, , drop = FALSE]

    tibble::as_tibble(merged)
}

.merge_barcode_info <- function(x, y, barcodes_retained, cell_join) {
    join_fun <- switch(
        cell_join,
        "union" = dplyr::full_join,
        "intersect" = dplyr::inner_join,
        "left" = dplyr::left_join,
        "right" = dplyr::right_join
    )

    auto_cols <- c("library_size", "non_zero_snps")
    barcode_info_x <- x@barcode_info %>% dplyr::select(-dplyr::any_of(auto_cols))
    barcode_info_y <- y@barcode_info %>% dplyr::select(-dplyr::any_of(auto_cols))

    # Check if barcode column exists in both objects
    has_barcode_x <- "barcode" %in% colnames(barcode_info_x)
    has_barcode_y <- "barcode" %in% colnames(barcode_info_y)

    # Determine join key: use barcode if available, otherwise fall back to cell_id
    join_by <- if (has_barcode_x && has_barcode_y) "barcode" else "cell_id"

    merged <- join_fun(
        barcode_info_x,
        barcode_info_y,
        by = join_by,
        suffix = c(".x", ".y")
    )

    x_conflicts <- grep("\\.x$", colnames(merged), value = TRUE)
    for (x_col in x_conflicts) {
        base_col <- sub("\\.x$", "", x_col)
        y_col <- paste0(base_col, ".y")
        merged[[base_col]] <- dplyr::coalesce(merged[[x_col]], merged[[y_col]])
        merged[[x_col]] <- NULL
        merged[[y_col]] <- NULL
    }

    # Filter and order by barcodes_retained
    # Use the appropriate column for filtering
    filter_col <- if (join_by == "barcode") "barcode" else "cell_id"
    merged <- merged[merged[[filter_col]] %in% barcodes_retained, , drop = FALSE]
    order_idx <- match(barcodes_retained, merged[[filter_col]])
    if (any(is.na(order_idx))) {
        stop(paste0("Some ", filter_col, "s could not be matched during merge_barcode_info()"))
    }
    merged <- merged[order_idx, , drop = FALSE]

    # Regenerate cell_id as sequential identifiers if we joined by barcode
    if (join_by == "barcode") {
        merged$cell_id <- paste0("cell_", seq_len(nrow(merged)))
    }

    tibble::as_tibble(merged)
}

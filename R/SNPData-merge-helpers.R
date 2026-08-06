# Merge helpers for SNPData

# Retained position (1-based) of every entry of `names_vec` in `retained`,
# or NA where it is not retained. A hash-based lookup shared by both the
# row and (name-based) column mapping in .merge_count_matrix().
.retained_position <- function(names_vec, retained) {
    lookup <- setNames(seq_along(retained), retained)
    unname(lookup[names_vec])
}

# Merges one allele channel (ref/alt/oth) of x and y directly from their
# triplet forms into the retained dimensions, rather than expanding x and y
# to the merged dimensions separately and adding the results -- sparseMatrix()
# sums duplicate (i, j) entries automatically, so concatenating both objects'
# remapped triplets and building the merged matrix in one call does the same
# work as expand-then-add without materialising two full-size intermediates.
# `col_target_x`/`col_target_y` are either name-based lookups (built the same
# way as `row_target_x`/`row_target_y`) or, in the barcode-matched case, the
# caller's precomputed `match(barcodes, barcodes_retained)` mapping -- both
# give the retained column position (or NA) for each column of the matrix.
.merge_count_matrix <- function(
    mat_x,
    mat_y,
    row_target_x,
    col_target_x,
    row_target_y,
    col_target_y,
    retained_rows,
    retained_cols
) {
    mat_x_T <- methods::as(mat_x, "TsparseMatrix")
    mat_y_T <- methods::as(mat_y, "TsparseMatrix")

    # Remap indices (TsparseMatrix uses 0-based indexing, but slot access gives 1-based)
    new_i_x <- row_target_x[mat_x_T@i + 1L]
    new_j_x <- col_target_x[mat_x_T@j + 1L]
    keep_x <- !is.na(new_i_x) & !is.na(new_j_x)

    new_i_y <- row_target_y[mat_y_T@i + 1L]
    new_j_y <- col_target_y[mat_y_T@j + 1L]
    keep_y <- !is.na(new_i_y) & !is.na(new_j_y)

    merged <- Matrix::sparseMatrix(
        i = c(new_i_x[keep_x], new_i_y[keep_y]),
        j = c(new_j_x[keep_x], new_j_y[keep_y]),
        x = c(mat_x_T@x[keep_x], mat_y_T@x[keep_y]),
        dims = c(length(retained_rows), length(retained_cols))
    )

    rownames(merged) <- retained_rows
    colnames(merged) <- retained_cols

    merged
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

.merge_donor_info <- function(x, y, donors_retained) {
    auto_cols <- c("n_cells")
    donor_info_x <- donor_info(x) %>% dplyr::select(-dplyr::any_of(auto_cols))
    donor_info_y <- donor_info(y) %>% dplyr::select(-dplyr::any_of(auto_cols))

    merged <- dplyr::full_join(
        donor_info_x,
        donor_info_y,
        by = "donor",
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

    merged <- merged[merged$donor %in% donors_retained, , drop = FALSE]
    order_idx <- match(donors_retained, merged$donor)
    if (any(is.na(order_idx))) {
        stop("Some donors could not be matched during merge_donor_info()")
    }
    merged <- merged[order_idx, , drop = FALSE]

    tibble::as_tibble(merged)
}

# Merging donor_snp_info is not a plain coalesce: the same (snp_id, donor,
# zygosity_source) call disagreeing between x and y on a "hazard" column
# (zygosity, the XCI phase, or informativeness) is a genuine reproducibility
# problem and must stop the merge rather than silently pick a side. Numeric
# estimate columns (p-values, GT posteriors, escape fractions) are exempt:
# two independently fit/tested values can differ by floating point alone,
# which is not a hazard, so those keep the x-priority coalesce used
# elsewhere in this file. zygosity_source is part of the join key (a pair
# may carry one row per source), not a value to compare -- the object-level
# *active* source is a separate concern, resolved by .merge_zygosity_source().
.merge_donor_snp_info <- function(x, y, snp_ids_retained, donors_retained) {
    hazard_cols <- c("zygosity", "allele_on_x1", "xci_informative")
    key_cols <- c("snp_id", "donor", "zygosity_source")

    dsi_x <- donor_snp_info(x, source = "all")
    dsi_y <- donor_snp_info(y, source = "all")
    value_cols <- setdiff(colnames(dsi_x), key_cols)

    merged <- dplyr::full_join(
        dsi_x,
        dsi_y,
        by = key_cols,
        suffix = c(".x", ".y")
    )

    for (col in value_cols) {
        x_col <- paste0(col, ".x")
        y_col <- paste0(col, ".y")
        if (!all(c(x_col, y_col) %in% colnames(merged))) {
            next
        }

        if (col %in% hazard_cols) {
            conflict <- !is.na(merged[[x_col]]) & !is.na(merged[[y_col]]) & merged[[x_col]] != merged[[y_col]]
            if (any(conflict)) {
                conflict_labels <- paste0(merged$snp_id[conflict], "/", merged$donor[conflict])
                conflict_msg <- paste(head(conflict_labels, 5), collapse = ", ")
                if (length(conflict_labels) > 5) {
                    conflict_msg <- paste0(conflict_msg, ", ... (", length(conflict_labels), " conflicts)")
                }
                stop(sprintf(
                    "merge_snpdata(): conflicting '%s' values in donor_snp_info for snp_id/donor %s. %s",
                    col,
                    conflict_msg,
                    "Resolve the disagreement before merging."
                ))
            }
        }

        merged[[col]] <- dplyr::coalesce(merged[[x_col]], merged[[y_col]])
        merged[[x_col]] <- NULL
        merged[[y_col]] <- NULL
    }

    merged <- merged[merged$snp_id %in% snp_ids_retained & merged$donor %in% donors_retained, , drop = FALSE]
    tibble::as_tibble(merged)
}

# Resolves the object-level *active* zygosity_source slot for a merge, using
# the same error-on-conflict/coalesce-on-agreement philosophy as the hazard
# columns above: if both objects have an active source and they differ,
# merging would silently change which source drives downstream analysis, so
# it errors instead. If only one is set, that one carries over.
.merge_zygosity_source <- function(x, y) {
    source_x <- zygosity_source(x)
    source_y <- zygosity_source(y)

    if (!is.na(source_x) && !is.na(source_y) && source_x != source_y) {
        stop(sprintf(
            "merge_snpdata(): conflicting active zygosity_source ('%s' vs '%s'). Resolve with zygosity_source<-() before merging.",
            source_x,
            source_y
        ))
    }

    dplyr::coalesce(source_x, source_y)
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

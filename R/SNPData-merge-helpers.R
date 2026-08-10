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

# Decides, for every cell present in both x and y, which object's copy survives.
#
# Summing a collided cell's counts is right only when the two columns really are
# the same physical cell sequenced twice. Across independently barcoded libraries
# they are not: a shared barcode string is a coincidence between two different
# cells, often from different donors, and summing them fabricates a cell that
# expresses both donors' alleles (see merge_snpdata()'s `cell_collision`).
#
# The loser is suppressed by setting its column target to NA, which
# .merge_count_matrix() already drops -- so no counts from the losing column
# reach the merged matrix.
#
# Depth is the full-object library size (ref + alt over every SNP, the same
# quantity as `barcode_info$library_size`), not just the retained SNPs: which of
# two cells is better sequenced is a property of the cell, and making it depend
# on `snp_join` would let the SNP-side join silently change which cells survive.
#
# @return A list with the two adjusted column-target vectors, plus `collided`
#   and `winner_is_y` -- both logical, indexed by retained column position, for
#   .merge_barcode_info() to resolve metadata to the same winner.
.resolve_cell_collisions <- function(x, y, col_target_x, col_target_y, n_retained) {
    collided <- rep(FALSE, n_retained)
    winner_is_y <- rep(FALSE, n_retained)

    in_x <- !is.na(col_target_x)
    in_y <- !is.na(col_target_y)
    collided[intersect(col_target_x[in_x], col_target_y[in_y])] <- TRUE

    if (!any(collided)) {
        return(list(
            col_target_x = col_target_x,
            col_target_y = col_target_y,
            collided = collided,
            winner_is_y = winner_is_y
        ))
    }

    # Depth per retained position, NA where that object does not contribute.
    depth_x <- rep(NA_real_, n_retained)
    depth_y <- rep(NA_real_, n_retained)
    depth_x[col_target_x[in_x]] <- (Matrix::colSums(x@ref_count) + Matrix::colSums(x@alt_count))[in_x]
    depth_y[col_target_y[in_y]] <- (Matrix::colSums(y@ref_count) + Matrix::colSums(y@alt_count))[in_y]

    # Strict `>` so an exact tie goes to x, the left-hand object.
    winner_is_y <- collided & !is.na(depth_x) & !is.na(depth_y) & depth_y > depth_x

    target_x <- col_target_x[in_x]
    col_target_x[in_x] <- ifelse(winner_is_y[target_x], NA_integer_, target_x)
    target_y <- col_target_y[in_y]
    col_target_y[in_y] <- ifelse(collided[target_y] & !winner_is_y[target_y], NA_integer_, target_y)

    logger::log_warn(
        "{sum(collided)} cell(s) share an identifier between the two objects; ",
        "keeping the higher-coverage copy ({sum(winner_is_y)} from y, ",
        "{sum(collided) - sum(winner_is_y)} from x) and discarding the other. ",
        "Pass cell_collision = 'sum' if these are technical replicates of the same cells."
    )

    list(
        col_target_x = col_target_x,
        col_target_y = col_target_y,
        collided = collided,
        winner_is_y = winner_is_y
    )
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

# `winner_is_y` is the per-retained-cell decision made by
# .resolve_cell_collisions(), or NULL under cell_collision = "sum". When it is
# supplied, a collided cell's metadata must follow the same winner as its
# counts: coalescing x-first would otherwise label the surviving column with the
# discarded cell's donor, which is the very mislabelling the resolution exists
# to prevent.
.merge_barcode_info <- function(x, y, barcodes_retained, cell_join, winner_is_y = NULL) {
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

    # Filter and order by barcodes_retained
    # Use the appropriate column for filtering
    # Ordering happens before the conflicting columns are resolved so that
    # `winner_is_y`, which is indexed by retained position, lines up row-wise.
    filter_col <- if (join_by == "barcode") "barcode" else "cell_id"
    merged <- merged[merged[[filter_col]] %in% barcodes_retained, , drop = FALSE]
    order_idx <- match(barcodes_retained, merged[[filter_col]])
    if (any(is.na(order_idx))) {
        stop(paste0("Some ", filter_col, "s could not be matched during merge_barcode_info()"))
    }
    merged <- merged[order_idx, , drop = FALSE]

    x_conflicts <- grep("\\.x$", colnames(merged), value = TRUE)
    for (x_col in x_conflicts) {
        base_col <- sub("\\.x$", "", x_col)
        y_col <- paste0(base_col, ".y")
        merged[[base_col]] <- if (is.null(winner_is_y)) {
            dplyr::coalesce(merged[[x_col]], merged[[y_col]])
        } else {
            # A non-collided cell contributes from one object only, so the
            # coalesce is still what resolves it; a collided one takes its
            # winner's value outright, with no fallback to the discarded copy.
            dplyr::if_else(winner_is_y, merged[[y_col]], dplyr::coalesce(merged[[x_col]], merged[[y_col]]))
        }
        merged[[x_col]] <- NULL
        merged[[y_col]] <- NULL
    }

    # Regenerate cell_id as sequential identifiers if we joined by barcode
    if (join_by == "barcode") {
        merged$cell_id <- paste0("cell_", seq_len(nrow(merged)))
    }

    tibble::as_tibble(merged)
}

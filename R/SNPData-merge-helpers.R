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

# The composite cell identity used throughout a merge: a barcode is unique only
# within its own library, so (library_id, barcode) is what actually names a cell.
#
# Pasted into a single string rather than kept as a pair because the retained
# cells are carried around as a character vector and looked up with match(). The
# separator is the ASCII unit separator, which cannot occur in a 10x barcode and
# is vanishingly unlikely in a library label, so two distinct pairs cannot
# collapse onto one key.
.CELL_KEY_SEP <- "\u001f"

.cell_keys <- function(barcode_info, id_col) {
    paste(barcode_info$library_id, barcode_info[[id_col]], sep = .CELL_KEY_SEP)
}

# The identifier half of each retained key, recovered from the inputs rather
# than by splitting the key, so a separator appearing inside a label cannot
# corrupt the result.
.key_identifiers <- function(keys_retained, all_keys, all_ids) {
    all_ids[match(keys_retained, all_keys)]
}

# Every cell must carry a library label before it can be merged.
#
# Without one there is no way to tell a barcode shared *within* a library (the
# same physical cell, whose counts should be summed) from one shared *across*
# libraries (two different cells, often from different donors, whose counts must
# not be). 10x barcodes come from a whitelist of roughly 737,000, so two runs of
# n cells are expected to share about n^2 / 737000 barcodes by chance -- around
# 34 for two 5,000-cell runs. Summing such a pair fabricates a cell carrying both
# donors' alleles, which at a SNP where their phase differs reads as biallelic
# expression: for X-linked genes that is indistinguishable from escape from
# X-chromosome inactivation, so the artefact lands precisely on the signal
# assign_xci() and test_escape() exist to measure. Guessing is not safe, so an
# unlabelled object is refused outright.
.check_library_ids <- function(x, y) {
    n_missing_x <- sum(is.na(x@barcode_info$library_id))
    n_missing_y <- sum(is.na(y@barcode_info$library_id))

    if (n_missing_x == 0 && n_missing_y == 0) {
        return(invisible(NULL))
    }

    stop(sprintf(
        paste0(
            "merge_snpdata() requires a library_id for every cell, but %d cell(s) in x and %d in y have none. ",
            "Cells are matched on (library_id, barcode) because a barcode only identifies a cell within its own ",
            "library. Set it at import with import_cellsnp(..., library_id = ) or afterwards with ",
            "barcode_info(obj)$library_id <- ."
        ),
        n_missing_x,
        n_missing_y
    ))
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

# The merged object's library_info is already derived from its cells, so only
# the BAM paths need carrying over. They are unioned rather than given to one
# side: two cellSNP runs over the same library legitimately share a BAM, and a
# library present in both objects has the same reads behind it either way. A
# path that is genuinely wrong is caught where it is used, by
# add_molecule_phase()'s existence and duplicate checks, rather than guessed at
# here.
.merge_library_bams <- function(merged, x, y) {
    if (nrow(merged@library_info) == 0) {
        return(merged)
    }
    from_x <- .lookup_bam_files(merged@library_info$library_id, library_info(x))
    from_y <- .lookup_bam_files(merged@library_info$library_id, library_info(y))
    merged@library_info$bam_files <- purrr::map2(from_x, from_y, union)
    merged
}

# The map is a property of the genome annotation, not of either object's cells,
# so the two sides' rows are unioned and cut to the SNPs the merge retained. An
# object whose annotation lacked strand contributes nothing, leaving the other's
# rows intact rather than blanking them.
.merge_snp_gene_map <- function(merged, x, y) {
    combined <- dplyr::distinct(dplyr::bind_rows(snp_gene_map(x), snp_gene_map(y)))
    merged@snp_gene_map <- combined[combined$snp_id %in% merged@snp_info$snp_id, , drop = FALSE]
    merged
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

# Metadata follows the same (library_id, `id_col`) identity the counts do, so a
# row's annotation always describes the cell whose counts sit in that column. A
# cell that both objects carry contributes one row and resolves conflicting
# columns x-first; two cells that merely share a barcode across libraries never
# meet, since the library halves of their keys differ.
#
# `cell_ids_retained` is decided by the caller (which knows whether the retained
# identifiers stayed unique) rather than regenerated here, so barcode_info$cell_id
# and the count matrices' column names cannot drift apart.
.merge_barcode_info <- function(x, y, keys_retained, cell_join, id_col, cell_ids_retained) {
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

    merged <- join_fun(
        barcode_info_x,
        barcode_info_y,
        by = c("library_id", id_col),
        suffix = c(".x", ".y")
    )

    merged_keys <- .cell_keys(merged, id_col)
    merged <- merged[merged_keys %in% keys_retained, , drop = FALSE]
    order_idx <- match(keys_retained, .cell_keys(merged, id_col))
    if (any(is.na(order_idx))) {
        stop("Some cells could not be matched during merge_barcode_info()")
    }
    merged <- merged[order_idx, , drop = FALSE]

    x_conflicts <- grep("\\.x$", colnames(merged), value = TRUE)
    for (x_col in x_conflicts) {
        base_col <- sub("\\.x$", "", x_col)
        y_col <- paste0(base_col, ".y")
        merged[[base_col]] <- dplyr::coalesce(merged[[x_col]], merged[[y_col]])
        merged[[x_col]] <- NULL
        merged[[y_col]] <- NULL
    }

    merged$cell_id <- cell_ids_retained

    tibble::as_tibble(merged)
}

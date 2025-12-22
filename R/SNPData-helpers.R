# Internal helpers for SNPData construction and validation

.validate_count_dims <- function(ref_count, alt_count, oth_count) {
    stopifnot(nrow(alt_count) == nrow(ref_count))
    stopifnot(ncol(alt_count) == ncol(ref_count))

    if (is.null(oth_count)) {
        oth_count <- Matrix::Matrix(
            0,
            nrow = nrow(ref_count),
            ncol = ncol(ref_count),
            sparse = TRUE
        )
    } else {
        stopifnot(nrow(oth_count) == nrow(ref_count))
        stopifnot(ncol(oth_count) == ncol(ref_count))
    }

    oth_count
}

.validate_info_dims <- function(ref_count, alt_count, snp_info, barcode_info) {
    stopifnot(ncol(alt_count) == nrow(barcode_info))
    stopifnot(nrow(ref_count) == nrow(snp_info))
}

.assign_snp_ids <- function(snp_info) {
    if (!"snp_id" %in% colnames(snp_info)) {
        if (all(c("chrom", "pos", "ref", "alt") %in% colnames(snp_info))) {
            snp_info$snp_id <- make_snp_id(
                snp_info$chrom,
                snp_info$pos,
                snp_info$ref,
                snp_info$alt
            )
        } else {
            snp_info$snp_id <- paste0("snp_", seq_len(nrow(snp_info)))
        }
    }

    snp_info
}

.assign_cell_ids <- function(barcode_info) {
    if (!"cell_id" %in% colnames(barcode_info)) {
        barcode_info$cell_id <- paste0("cell_", seq_len(nrow(barcode_info)))
    }

    barcode_info
}

.dedupe_snps <- function(ref_count, alt_count, oth_count, snp_info) {
    if (!any(duplicated(snp_info$snp_id))) {
        return(list(
            ref_count = ref_count,
            alt_count = alt_count,
            oth_count = oth_count,
            snp_info = snp_info
        ))
    }

    dup_positions <- which(duplicated(snp_info$snp_id))
    dup_labels <- paste0(
        snp_info$snp_id[dup_positions],
        " (row ",
        dup_positions,
        ")"
    )
    dup_snps_msg <- paste(head(dup_labels, 5), collapse = ", ")
    if (length(dup_positions) > 5) {
        dup_snps_msg <- paste0(
            dup_snps_msg,
            ", ... (",
            length(dup_positions),
            " duplicates)"
        )
    }

    warning(
        sprintf(
            "Duplicate SNP IDs detected (%s). Keeping first occurrence and dropping duplicates.",
            dup_snps_msg
        ),
        call. = FALSE
    )

    keep_snps <- !duplicated(snp_info$snp_id)
    list(
        ref_count = ref_count[keep_snps, , drop = FALSE],
        alt_count = alt_count[keep_snps, , drop = FALSE],
        oth_count = oth_count[keep_snps, , drop = FALSE],
        snp_info = snp_info[keep_snps, , drop = FALSE]
    )
}

.set_dimnames <- function(ref_count, alt_count, oth_count, snp_info, barcode_info) {
    colnames(ref_count) <- barcode_info$cell_id
    colnames(alt_count) <- barcode_info$cell_id
    colnames(oth_count) <- barcode_info$cell_id
    rownames(ref_count) <- snp_info$snp_id
    rownames(alt_count) <- snp_info$snp_id
    rownames(oth_count) <- snp_info$snp_id

    list(
        ref_count = ref_count,
        alt_count = alt_count,
        oth_count = oth_count
    )
}

.recompute_snp_stats <- function(snp_info, ref_count, alt_count) {
    snp_info$coverage <- Matrix::rowSums(alt_count + ref_count)
    snp_info$non_zero_samples <- Matrix::rowSums(alt_count + ref_count > 0)
    snp_info
}

.recompute_barcode_stats <- function(barcode_info, ref_count, alt_count) {
    barcode_info$library_size <- Matrix::colSums(alt_count + ref_count)
    barcode_info$non_zero_snps <- Matrix::colSums(alt_count + ref_count > 0)
    barcode_info
}

.recompute_metrics <- function(snp_info, barcode_info, ref_count, alt_count) {
    list(
        snp_info = .recompute_snp_stats(snp_info, ref_count, alt_count),
        barcode_info = .recompute_barcode_stats(barcode_info, ref_count, alt_count)
    )
}

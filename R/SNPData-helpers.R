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

.validate_donor_dims <- function(donor_info, donor_snp_info, snp_info) {
    problems <- character(0)

    if (!"donor" %in% colnames(donor_info)) {
        problems <- c(problems, "donor_info must contain a 'donor' column")
    } else if (any(duplicated(donor_info$donor))) {
        problems <- c(problems, "donor_info contains duplicate donor values")
    }

    missing_cols <- setdiff(c("snp_id", "donor"), colnames(donor_snp_info))
    if (length(missing_cols) > 0) {
        problems <- c(
            problems,
            paste0(
                "donor_snp_info is missing required column(s): ",
                paste(missing_cols, collapse = ", ")
            )
        )
        return(problems)
    }

    # A (snp_id, donor) pair may carry one row per zygosity_source (e.g. a
    # "vireo_gt" row and a "binomial" row coexisting) -- so the key only
    # widens to include zygosity_source when that column is actually present.
    dup_key <- if ("zygosity_source" %in% colnames(donor_snp_info)) {
        c("snp_id", "donor", "zygosity_source")
    } else {
        c("snp_id", "donor")
    }
    if (any(duplicated(donor_snp_info[dup_key]))) {
        problems <- c(problems, sprintf("donor_snp_info contains duplicate (%s) rows", paste(dup_key, collapse = ", ")))
    }

    unknown_snps <- setdiff(donor_snp_info$snp_id, snp_info$snp_id)
    if (length(unknown_snps) > 0) {
        problems <- c(
            problems,
            sprintf(
                "donor_snp_info references %d snp_id(s) not present in snp_info",
                length(unknown_snps)
            )
        )
    }

    unknown_donors <- setdiff(donor_snp_info$donor, donor_info$donor)
    if (length(unknown_donors) > 0) {
        problems <- c(
            problems,
            sprintf(
                "donor_snp_info references %d donor(s) not present in donor_info",
                length(unknown_donors)
            )
        )
    }

    if (all(c("zygosity", "zygosity_source") %in% colnames(donor_snp_info))) {
        missing_source <- !is.na(donor_snp_info$zygosity) & is.na(donor_snp_info$zygosity_source)
        if (any(missing_source)) {
            problems <- c(
                problems,
                sprintf(
                    "donor_snp_info has %d row(s) with a zygosity call but no zygosity_source",
                    sum(missing_source)
                )
            )
        }
    }

    problems
}

.empty_donor_snp_info <- function() {
    tibble::tibble(
        snp_id = character(0),
        donor = character(0),
        zygosity = character(0),
        zygosity_source = character(0),
        zygosity_p_val = double(0),
        zygosity_adj_p_val = double(0),
        zygosity_gt_prob = double(0),
        xci_informative = logical(0),
        allele_on_x1 = character(0),
        xci_escape_fraction = double(0)
    )
}

.apply_donor_map <- function(donor, donor_map) {
    if (is.null(donor_map)) {
        return(donor)
    }
    if (
        !is.character(donor_map) ||
            is.null(names(donor_map)) ||
            any(names(donor_map) == "") ||
            any(duplicated(names(donor_map)))
    ) {
        stop(
            "donor_map must be a named character vector (new donor label = old donor label) with unique, non-empty names"
        )
    }
    if (any(duplicated(donor_map))) {
        stop("donor_map must not map more than one new label to the same old donor")
    }

    idx <- match(donor, donor_map)
    matched <- !is.na(idx)
    # Not `ifelse()`: on zero-length input it returns `logical(0)` regardless
    # of `donor`'s type (a base R quirk), which breaks a downstream join
    # expecting `donor` to stay character.
    result <- donor
    result[matched] <- names(donor_map)[idx[matched]]
    result
}

.derive_zygosity_source <- function(donor_snp_info) {
    if (!"zygosity_source" %in% colnames(donor_snp_info)) {
        return(NA_character_)
    }
    sources <- unique(stats::na.omit(donor_snp_info$zygosity_source))
    if (length(sources) == 1) sources else NA_character_
}

.default_donor_info <- function(barcode_info) {
    if ("donor" %in% colnames(barcode_info)) {
        tibble::tibble(donor = sort(unique(stats::na.omit(barcode_info$donor))))
    } else {
        tibble::tibble(donor = character(0))
    }
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

.dedupe_snps <- function(ref_count, alt_count, oth_count, snp_info, donor_snp_info) {
    if (!any(duplicated(snp_info$snp_id))) {
        return(list(
            ref_count = ref_count,
            alt_count = alt_count,
            oth_count = oth_count,
            snp_info = snp_info,
            donor_snp_info = donor_snp_info
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
    kept_snp_info <- snp_info[keep_snps, , drop = FALSE]
    list(
        ref_count = ref_count[keep_snps, , drop = FALSE],
        alt_count = alt_count[keep_snps, , drop = FALSE],
        oth_count = oth_count[keep_snps, , drop = FALSE],
        snp_info = kept_snp_info,
        donor_snp_info = donor_snp_info[donor_snp_info$snp_id %in% kept_snp_info$snp_id, , drop = FALSE]
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

.recompute_donor_stats <- function(donor_info, barcode_info) {
    if (nrow(donor_info) == 0) {
        donor_info$n_cells <- integer(0)
        return(donor_info)
    }

    if ("donor" %in% colnames(barcode_info)) {
        cell_counts <- table(barcode_info$donor)
        donor_info$n_cells <- as.integer(cell_counts[donor_info$donor])
        donor_info$n_cells[is.na(donor_info$n_cells)] <- 0L
    } else {
        donor_info$n_cells <- 0L
    }

    donor_info
}

.recompute_metrics <- function(snp_info, barcode_info, donor_info, ref_count, alt_count) {
    list(
        snp_info = .recompute_snp_stats(snp_info, ref_count, alt_count),
        barcode_info = .recompute_barcode_stats(barcode_info, ref_count, alt_count),
        donor_info = .recompute_donor_stats(donor_info, barcode_info)
    )
}

#' Get barcode/cell-level SNP count summary
#'
#' Returns a long-format data frame of reference and alternate allele counts per SNP and cell/barcode.
#'
#' @param x A SNPData object, required.
#' @param test_maf Logical (default \code{TRUE}). Whether to include a test_maf column.
#' @return A tibble with columns: snp_id, gene_name, chrom, pos, strand (if available in snp_info), cell_id, ref_count, alt_count, total_count, ref_ratio, maf, (optionally test_maf)
#' @export
#'
#' @examples
#' \dontrun{
#' snp_data <- get_example_snpdata()
#' barcode_count_df(snp_data)
#' }
#' @family count aggregation functions
#' @rdname barcode_count_df
setGeneric("barcode_count_df", function(x, test_maf = TRUE) standardGeneric("barcode_count_df"))

.build_long_count_df <- function(ref_count_mat, alt_count_mat, group_name, snp_info = NULL) {
    ref_count_df <- ref_count_mat %>%
        as.matrix() %>%
        tibble::as_tibble() %>%
        dplyr::mutate(snp_id = rownames(ref_count_mat), .before = 1) %>%
        tidyr::pivot_longer(
            cols = -snp_id,
            names_to = group_name,
            values_to = "ref_count"
        )

    alt_count_df <- alt_count_mat %>%
        as.matrix() %>%
        tibble::as_tibble() %>%
        dplyr::mutate(snp_id = rownames(alt_count_mat), .before = 1) %>%
        tidyr::pivot_longer(
            cols = -snp_id,
            names_to = group_name,
            values_to = "alt_count"
        )

    result <- dplyr::inner_join(ref_count_df, alt_count_df, by = c("snp_id", group_name)) %>%
        dplyr::mutate(
            total_count = ref_count + alt_count,
            ref_ratio = ref_count / total_count,
            maf = pmin(ref_count, alt_count) / total_count
        ) %>%
        dplyr::filter(total_count > 0)

    # Join SNP metadata if provided
    if (!is.null(snp_info)) {
        # Select relevant columns from snp_info
        snp_cols <- c("snp_id")
        optional_cols <- c("gene_name", "chrom", "pos", "strand")
        available_cols <- optional_cols[optional_cols %in% colnames(snp_info)]
        snp_metadata <- snp_info %>%
            dplyr::select(dplyr::all_of(c(snp_cols, available_cols)))

        result <- result %>%
            dplyr::left_join(snp_metadata, by = "snp_id") %>%
            dplyr::relocate(dplyr::all_of(available_cols), .after = snp_id)
    }

    result
}

.prepare_grouped_metadata <- function(barcode_info, group_by) {
    missing_cols <- group_by[!group_by %in% colnames(barcode_info)]
    if (length(missing_cols) > 0) {
        if (identical(group_by, "donor")) {
            missing_message <- "Donor information not available. Add donor data using add_barcode_metadata() or import_cellsnp() with vireo_folder parameter."
        } else if (identical(group_by, "clonotype")) {
            missing_message <- "Clonotype information not available. Add clonotype data using add_barcode_metadata() or import_cellsnp() with vdj_file parameter."
        } else if (length(group_by) == 1) {
            missing_message <- glue::glue("Column '{group_by}' not found in barcode_info.")
        } else {
            missing_message <- glue::glue(
                "Column(s) '{paste(missing_cols, collapse = \"', '\")}' not found in barcode_info."
            )
        }
        stop(missing_message)
    }

    na_mask <- Reduce(`|`, lapply(group_by, function(col) is.na(barcode_info[[col]])))
    if (any(na_mask)) {
        na_count <- sum(na_mask)
        logger::log_warn(
            "Found {na_count} rows with NA values in '{paste(group_by, collapse = \"', '\")}' column(s). These will be excluded from aggregation."
        )
        barcode_info <- barcode_info[!na_mask, , drop = FALSE]
    }

    if ("donor" %in% group_by) {
        keep_donor <- !barcode_info$donor %in% c("unassigned", "doublet")
        barcode_info <- barcode_info[keep_donor, , drop = FALSE]
    }

    if (nrow(barcode_info) == 0) {
        if (identical(group_by, "donor")) {
            all_na_message <- "All donor values are NA. Cannot perform donor-level aggregation. Add donor data using add_barcode_metadata() or import_cellsnp() with vireo_folder parameter."
        } else if (identical(group_by, "clonotype")) {
            all_na_message <- "All clonotype values are NA. Cannot perform clonotype-level aggregation. Add clonotype data using add_barcode_metadata() or import_cellsnp() with vdj_file parameter."
        } else {
            all_na_message <- glue::glue("All values are NA for '{paste(group_by, collapse = \"', '\")}'.")
        }
        stop(all_na_message)
    }

    group_key <- as.character(as.integer(interaction(barcode_info[group_by], drop = TRUE)))
    group_key_lookup <- dplyr::bind_cols(barcode_info[group_by], group_key = group_key) %>%
        dplyr::distinct()

    list(barcode_info = barcode_info, groups = group_key, group_key_lookup = group_key_lookup)
}

.aggregate_grouped_counts <- function(x, group_by, test_maf = TRUE, add_most_likely_donor = FALSE) {
    if (add_most_likely_donor && !identical(group_by, "clonotype")) {
        stop("add_most_likely_donor is only supported when group_by is 'clonotype'.")
    }

    barcode_info <- barcode_info(x)
    metadata <- .prepare_grouped_metadata(barcode_info = barcode_info, group_by = group_by)

    logger::log_info("Extracting reference counts")
    ref_count_grouped <- groupedRowSums(
        ref_count(x[, metadata$barcode_info$cell_id, drop = FALSE]),
        metadata$groups
    )

    logger::log_info("Extracting alternate counts")
    alt_count_grouped <- groupedRowSums(
        alt_count(x[, metadata$barcode_info$cell_id, drop = FALSE]),
        metadata$groups
    )

    logger::log_info("Processing reference and alternate counts")
    out <- .build_long_count_df(ref_count_grouped, alt_count_grouped, "group_key", snp_info(x))
    out <- out %>%
        dplyr::left_join(metadata$group_key_lookup, by = "group_key") %>%
        dplyr::select(-group_key) %>%
        dplyr::relocate(dplyr::all_of(group_by), .after = snp_id)

    if (identical(group_by, "clonotype") && isTRUE(add_most_likely_donor)) {
        most_likely_donor <- barcode_info %>%
            dplyr::filter(!is.na(clonotype) & !is.na(donor)) %>%
            dplyr::select(clonotype, donor) %>%
            dplyr::count(donor, clonotype) %>%
            dplyr::arrange(dplyr::desc(n), .by = clonotype) %>%
            dplyr::slice(1, .by = clonotype) %>%
            dplyr::select(clonotype, donor)

        out <- out %>%
            dplyr::left_join(most_likely_donor, by = "clonotype")
    }

    if (isTRUE(test_maf)) {
        out <- test_maf(out)
    }

    group_label <- paste(stringr::str_to_title(group_by), collapse = " x ")
    logger::log_success("{group_label} level counts calculated")
    out
}

barcode_count_df_impl <- function(x, test_maf = TRUE) {
    logger::log_info("Calculating barcode/cell level counts")

    out <- .build_long_count_df(ref_count(x), alt_count(x), "cell_id", snp_info(x))
    if (isTRUE(test_maf)) {
        out <- test_maf(out)
    }

    logger::log_success("Barcode/cell level counts calculated")
    out
}

#' @rdname barcode_count_df
#' @include SNPData-class.R
setMethod(
    "barcode_count_df",
    signature(x = "SNPData"),
    barcode_count_df_impl
)

#' Get donor-level SNP count summary
#'
#' Returns a long-format data frame of reference and alternate allele counts per SNP and donor.
#'
#' @param x A SNPData object, required.
#' @param test_maf Logical (default \code{TRUE}). Whether to include a test_maf column.
#' @return A tibble with columns: snp_id, gene_name, chrom, pos, strand (if available in snp_info), donor, ref_count, alt_count, total_count, ref_ratio, maf, (optionally test_maf)
#' @export
#'
#' @examples
#' \dontrun{
#' snp_data <- get_example_snpdata()
#' donor_count_df(snp_data)
#' }
#' @family count aggregation functions
#' @rdname donor_count_df
setGeneric("donor_count_df", function(x, test_maf = TRUE) standardGeneric("donor_count_df"))

donor_count_df_impl <- function(x, test_maf = TRUE) {
    aggregate_count_df(x, "donor", test_maf = test_maf)
}

#' @rdname donor_count_df
setMethod(
    "donor_count_df",
    signature(x = "SNPData"),
    donor_count_df_impl
)


#' Get clonotype-level SNP count summary
#'
#' Returns a long-format data frame of reference and alternate allele counts per SNP and clonotype, with most likely donor annotation.
#'
#' @param x A SNPData object, required.
#' @param test_maf Logical (default \code{TRUE}). Whether to include a test_maf column.
#' @return A tibble with columns: snp_id, gene_name, chrom, pos, strand (if available in snp_info), clonotype, ref_count, alt_count, total_count, ref_ratio, maf, donor, (optionally test_maf)
#' @export
#'
#' @examples
#' \dontrun{
#' snp_data <- get_example_snpdata()
#' clonotype_count_df(snp_data)
#' }
#' @family count aggregation functions
#' @rdname clonotype_count_df
setGeneric("clonotype_count_df", function(x, test_maf = TRUE) standardGeneric("clonotype_count_df"))

clonotype_count_df_impl <- function(x, test_maf = TRUE) {
    aggregate_count_df(x, "clonotype", test_maf = test_maf)
}

#' @rdname clonotype_count_df
setMethod(
    "clonotype_count_df",
    signature(x = "SNPData"),
    clonotype_count_df_impl
)


#' Get aggregated SNP count summary by any barcode_info column
#'
#' Returns a long-format data frame of reference and alternate allele counts per SNP and aggregated by the specified grouping column from barcode_info.
#'
#' @param x A SNPData object, required.
#' @param group_by Character vector, required. Column name(s) in \code{barcode_info} to group by.
#'   When more than one column is given, counts are aggregated over each unique combination of
#'   values across those columns, and each grouping column is retained in the output.
#' @param test_maf Logical (default \code{TRUE}). Whether to include a test_maf column.
#' @return A tibble with columns: snp_id, gene_name, chrom, pos, strand (if available in snp_info), [group_by], ref_count, alt_count, total_count, ref_ratio, maf, (optionally test_maf)
#' @export
#'
#' @examples
#' \dontrun{
#' snp_data <- get_example_snpdata()
#' # Aggregate by donor
#' aggregate_count_df(snp_data, "donor")
#'
#' # Aggregate by clonotype
#' aggregate_count_df(snp_data, "clonotype")
#'
#' # Aggregate by donor and cell_type combined
#' aggregate_count_df(snp_data, c("donor", "cell_type"))
#' }
#' @family count aggregation functions
#' @rdname aggregate_count_df
setGeneric("aggregate_count_df", function(x, group_by, test_maf = TRUE) standardGeneric("aggregate_count_df"))

aggregate_count_df_impl <- function(x, group_by, test_maf = TRUE) {
    logger::log_info("Calculating {paste(group_by, collapse = ' x ')} level counts")
    if (identical(group_by, "clonotype")) {
        .aggregate_grouped_counts(
            x = x,
            group_by = group_by,
            test_maf = test_maf,
            add_most_likely_donor = TRUE
        )
    } else {
        .aggregate_grouped_counts(
            x = x,
            group_by = group_by,
            test_maf = test_maf,
            add_most_likely_donor = FALSE
        )
    }
}

#' @rdname aggregate_count_df
setMethod(
    "aggregate_count_df",
    signature(x = "SNPData"),
    aggregate_count_df_impl
)

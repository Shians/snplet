#' Get barcode/cell-level SNP count summary
#'
#' Returns a long-format data frame of reference and alternate allele counts per SNP and cell/barcode.
#'
#' @param x A SNPData object
#' @param test_maf Logical, whether to include a test_maf column (default TRUE)
#' @return A tibble with columns: snp_id, cell_id, ref_count, alt_count, total_count, ref_ratio, maf, (optionally test_maf)
#' @export
#'
#' @examples
#' \dontrun{
#' snp_data <- get_example_snpdata()
#' barcode_count_df(snp_data)
#' }
#' @rdname barcode_count_df
setGeneric("barcode_count_df", function(x, test_maf = TRUE) standardGeneric("barcode_count_df"))

.build_long_count_df <- function(ref_count_mat, alt_count_mat, group_name) {
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

    dplyr::inner_join(ref_count_df, alt_count_df, by = c("snp_id", group_name)) %>%
        dplyr::mutate(
            total_count = ref_count + alt_count,
            ref_ratio = ref_count / total_count,
            maf = pmin(ref_count, alt_count) / total_count
        ) %>%
        dplyr::filter(total_count > 0)
}

barcode_count_df_impl <- function(x, test_maf = TRUE) {
    logger::log_info("Calculating barcode/cell level counts")

    out <- .build_long_count_df(ref_count(x), alt_count(x), "cell_id")
    if (isTRUE(test_maf)) {
        out <- test_maf(out)
    }

    logger::log_success("Barcode/cell level counts calculated")
    out
}

#' @rdname barcode_count_df
setMethod(
    "barcode_count_df",
    signature(x = "SNPData"),
    barcode_count_df_impl
)

#' Get donor-level SNP count summary
#'
#' Returns a long-format data frame of reference and alternate allele counts per SNP and donor.
#'
#' @param x A SNPData object
#' @param test_maf Logical, whether to include a test_maf column (default TRUE)
#' @return A tibble with columns: snp_id, donor, ref_count, alt_count, total_count, ref_ratio, maf, (optionally test_maf)
#' @export
#'
#' @examples
#' \dontrun{
#' snp_data <- get_example_snpdata()
#' donor_count_df(snp_data)
#' }
#' @rdname donor_count_df
setGeneric("donor_count_df", function(x, test_maf = TRUE) standardGeneric("donor_count_df"))

donor_count_df_impl <- function(x, test_maf = TRUE) {
    logger::log_info("Calculating donor level counts")

    logger::log_info("Extracting reference counts")
    ref_count_grouped <- groupedRowSums(ref_count(x), get_barcode_info(x)$donor)
    ref_count_grouped <- ref_count_grouped[
        ,
        !colnames(ref_count_grouped) %in% c("unassigned", "doublet"),
        drop = FALSE
    ]

    logger::log_info("Extracting alternate counts")
    alt_count_grouped <- groupedRowSums(alt_count(x), get_barcode_info(x)$donor)
    alt_count_grouped <- alt_count_grouped[
        ,
        !colnames(alt_count_grouped) %in% c("unassigned", "doublet"),
        drop = FALSE
    ]

    logger::log_info("Processing reference and alternate counts")
    out <- .build_long_count_df(ref_count_grouped, alt_count_grouped, "donor")
    if (isTRUE(test_maf)) {
        out <- test_maf(out)
    }

    logger::log_success("Donor level counts calculated")
    out
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
#' @param x A SNPData object
#' @param test_maf Logical, whether to include a test_maf column (default TRUE)
#' @return A tibble with columns: snp_id, clonotype, ref_count, alt_count, total_count, ref_ratio, maf, donor, (optionally test_maf)
#' @export
#'
#' @examples
#' \dontrun{
#' snp_data <- get_example_snpdata()
#' clonotype_count_df(snp_data)
#' }
#' @rdname clonotype_count_df
setGeneric("clonotype_count_df", function(x, test_maf = TRUE) standardGeneric("clonotype_count_df"))

clonotype_count_df_impl <- function(x, test_maf = TRUE) {
    logger::log_info("Calculating clonotype level counts")

    # Check if clonotype column exists
    barcode_info <- get_barcode_info(x)
    if (!"clonotype" %in% colnames(barcode_info)) {
        stop("Clonotype information not available. Add clonotype data using add_barcode_metadata() or import_cellsnp() with vdj_file parameter.")
    }

    # Check if all clonotype values are NA
    if (all(is.na(barcode_info$clonotype))) {
        stop("All clonotype values are NA. Cannot perform clonotype-level aggregation. Add clonotype data using add_barcode_metadata() or import_cellsnp() with vdj_file parameter.")
    }

    logger::log_info("Extracting reference counts")
    ref_count_grouped <- groupedRowSums(ref_count(x), barcode_info$clonotype)

    logger::log_info("Extracting alternate counts")
    alt_count_grouped <- groupedRowSums(alt_count(x), barcode_info$clonotype)

    logger::log_info("Processing reference and alternate counts")
    most_likely_donor <- barcode_info %>%
        dplyr::filter(!is.na(clonotype) & !is.na(donor)) %>%
        dplyr::select(clonotype, donor) %>%
        dplyr::count(donor, clonotype) %>%
        dplyr::arrange(dplyr::desc(n), .by = clonotype) %>%
        dplyr::slice(1, .by = clonotype) %>%
        dplyr::select(clonotype, donor)

    out <- .build_long_count_df(ref_count_grouped, alt_count_grouped, "clonotype") %>%
        dplyr::left_join(most_likely_donor, by = c("clonotype" = "clonotype"))
    if (isTRUE(test_maf)) {
        out <- test_maf(out)
    }

    logger::log_success("Clonotype level counts calculated")
    out
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
#' @param x A SNPData object
#' @param group_by Character string specifying the column name in barcode_info to group by
#' @param test_maf Logical, whether to include a test_maf column (default TRUE)
#' @return A tibble with columns: snp_id, [group_by], ref_count, alt_count, total_count, ref_ratio, maf, (optionally test_maf)
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
#' }
#' @rdname aggregate_count_df
setGeneric("aggregate_count_df", function(x, group_by, test_maf = TRUE) standardGeneric("aggregate_count_df"))

aggregate_count_df_impl <- function(x, group_by, test_maf = TRUE) {
    logger::log_info("Calculating {group_by} level counts")

    # Check that group_by column exists in barcode_info
    if (!group_by %in% colnames(get_barcode_info(x))) {
        available_cols <- paste0(colnames(get_barcode_info(x)), collapse = ", ")
        stop(glue::glue("Column '{group_by}' not found in barcode_info. Available columns: {available_cols}"))
    }

    # Get grouping variable
    groups <- get_barcode_info(x)[[group_by]]

    # Check for missing values in grouping variable
    if (any(is.na(groups))) {
        logger::log_warn(
            "Found {sum(is.na(groups))} NA values in '{group_by}' column. These will be excluded from aggregation."
        )
        # Filter out samples with NA values
        keep_samples <- !is.na(groups)
        x_filtered <- x[, keep_samples]
        groups <- groups[keep_samples]
    } else {
        x_filtered <- x
    }

    logger::log_info("Extracting reference counts")
    ref_count_grouped <- groupedRowSums(ref_count(x_filtered), groups)

    logger::log_info("Extracting alternate counts")
    alt_count_grouped <- groupedRowSums(alt_count(x_filtered), groups)

    logger::log_info("Processing reference and alternate counts")
    out <- .build_long_count_df(ref_count_grouped, alt_count_grouped, group_by)
    if (isTRUE(test_maf)) {
        out <- test_maf(out)
    }

    logger::log_success("{stringr::str_to_title(group_by)} level counts calculated")
    out
}

#' @rdname aggregate_count_df
setMethod(
    "aggregate_count_df",
    signature(x = "SNPData"),
    aggregate_count_df_impl
)

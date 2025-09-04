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
setGeneric("barcode_count_df", function(x, test_maf = TRUE) standardGeneric("barcode_count_df"))
#' @rdname barcode_count_df
setMethod(
  "barcode_count_df", signature(x = "SNPData"),
  function(x, test_maf = TRUE) {
    logger::log_info("Calculating barcode/cell level counts")

    ref_count_mat <- ref_count(x)
    ref_count_df <- ref_count_mat %>%
      as.matrix() %>%
      tibble::as_tibble() %>%
      dplyr::mutate(snp_id = rownames(ref_count_mat), .before = 1) %>%
      tidyr::pivot_longer(
        cols = -snp_id,
        names_to = "cell_id",
        values_to = "ref_count"
      )

    alt_count_mat <- alt_count(x)
    alt_count_df <- alt_count_mat %>%
      as.matrix() %>%
      tibble::as_tibble() %>%
      dplyr::mutate(snp_id = rownames(alt_count_mat), .before = 1) %>%
      tidyr::pivot_longer(
        cols = -snp_id,
        names_to = "cell_id",
        values_to = "alt_count"
      )

    out <- dplyr::inner_join(ref_count_df, alt_count_df, by = c("snp_id", "cell_id")) %>%
      dplyr::mutate(
        total_count = ref_count + alt_count,
        ref_ratio = ref_count / total_count,
        maf = pmin(ref_count, alt_count) / total_count
      ) %>%
      dplyr::filter(total_count > 0)
    if (isTRUE(test_maf)) {
      out <- out %>% dplyr::mutate(test_maf = !is.na(maf) & maf > 0)
    }

    logger::log_success("Barcode/cell level counts calculated")
    out
  }
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
setGeneric("donor_count_df", function(x, test_maf = TRUE) standardGeneric("donor_count_df"))
#' @rdname donor_count_df
setMethod(
  "donor_count_df", signature(x = "SNPData"),
  function(x, test_maf = TRUE) {
    logger::log_info("Calculating donor level counts")

    logger::log_info("Extracting reference counts")
    ref_count_grouped <- groupedRowSums(ref_count(x), get_barcode_info(x)$donor)
    ref_count_df <- ref_count_grouped %>%
      tibble::as_tibble() %>%
      dplyr::select(-any_of(c("unassigned", "doublet"))) %>%
      dplyr::mutate(snp_id = rownames(ref_count_grouped), .before = 1) %>%
      tidyr::pivot_longer(contains("donor"), names_to = "donor", values_to = "ref_count")

    logger::log_info("Extracting alternate counts")
    alt_count_grouped <- groupedRowSums(alt_count(x), get_barcode_info(x)$donor)
    alt_count_df <- alt_count_grouped %>%
      tibble::as_tibble() %>%
      dplyr::select(-any_of(c("unassigned", "doublet"))) %>%
      dplyr::mutate(snp_id = rownames(alt_count_grouped), .before = 1) %>%
      tidyr::pivot_longer(contains("donor"), names_to = "donor", values_to = "alt_count")

    logger::log_info("Processing reference and alternate counts")
    out <- dplyr::inner_join(ref_count_df, alt_count_df, by = c("snp_id", "donor")) %>%
      dplyr::mutate(
        total_count = ref_count + alt_count,
        ref_ratio = ref_count / total_count,
        maf = pmin(ref_count, alt_count) / total_count
      ) %>%
      dplyr::filter(total_count > 0)
    if (isTRUE(test_maf)) {
      out <- out %>% dplyr::mutate(test_maf = !is.na(maf) & maf > 0)
    }

    logger::log_success("Donor level counts calculated")
    out
  }
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
setGeneric("clonotype_count_df", function(x, test_maf = TRUE) standardGeneric("clonotype_count_df"))
#' @rdname clonotype_count_df
setMethod(
  "clonotype_count_df", signature(x = "SNPData"),
  function(x, test_maf = TRUE) {
    logger::log_info("Calculating clonotype level counts")

    logger::log_info("Extracting reference counts")
    ref_count_grouped <- groupedRowSums(ref_count(x), get_barcode_info(x)$clonotype)
    ref_count_df <- ref_count_grouped %>%
      tibble::as_tibble() %>%
      dplyr::mutate(snp_id = rownames(ref_count_grouped), .before = 1) %>%
      tidyr::pivot_longer(contains("clonotype"), names_to = "clonotype", values_to = "ref_count")

    logger::log_info("Extracting alternate counts")
    alt_count_grouped <- groupedRowSums(alt_count(x), get_barcode_info(x)$clonotype)
    alt_count_df <- alt_count_grouped %>%
      tibble::as_tibble() %>%
      dplyr::mutate(snp_id = rownames(alt_count_grouped), .before = 1) %>%
      tidyr::pivot_longer(contains("clonotype"), names_to = "clonotype", values_to = "alt_count")

    logger::log_info("Processing reference and alternate counts")
    most_likely_donor <- get_barcode_info(x) %>%
      dplyr::filter(!is.na(clonotype) & !is.na(donor)) %>%
      dplyr::select(clonotype, donor) %>%
      dplyr::count(donor, clonotype) %>%
      dplyr::arrange(dplyr::desc(n), .by = clonotype) %>%
      dplyr::slice(1, .by = clonotype) %>%
      dplyr::select(clonotype, donor)

    out <- dplyr::inner_join(ref_count_df, alt_count_df, by = c("snp_id", "clonotype")) %>%
      dplyr::mutate(
        total_count = ref_count + alt_count,
        ref_ratio = ref_count / total_count,
        maf = pmin(ref_count, alt_count) / total_count
      ) %>%
      dplyr::left_join(
        most_likely_donor,
        by = c("clonotype" = "clonotype")
      ) %>%
      dplyr::filter(total_count > 0)
    if (isTRUE(test_maf)) {
      out <- out %>% dplyr::mutate(test_maf = !is.na(maf) & maf > 0)
    }

    logger::log_success("Clonotype level counts calculated")
    out
  }
)

#' Check that columns in filter expressions exist in the data.frame
#'
#' Helper function to validate that all columns referenced in filter expressions
#' exist in the target data.frame or in the parent environment. Used internally
#' by filtering functions to provide informative error messages.
#'
#' @param df A data.frame to check column existence against
#' @param dots A list of quosures containing filter expressions
#' @param df_name Character string naming the data.frame for error messages
#' @return Invisibly returns NULL if all columns exist, otherwise throws an error
#' @keywords internal
check_filter_expr <- function(df, dots, df_name = "data.frame") {
  vars <- unique(unlist(lapply(dots, function(q) all.vars(rlang::get_expr(q)))))
  missing_vars <- setdiff(vars, colnames(df))
  # Only error if missing_vars are not found in parent.frame()
  still_missing <- missing_vars[!vapply(missing_vars, exists, logical(1), envir = parent.frame())]
  if (length(still_missing) > 0) {
    stop(paste0(
      "The following columns are not present in ", df_name, " or parent environment: ",
      paste(still_missing, collapse = ", ")
    ))
  }
}

#' Filter SNPData object by SNP information
#'
#' @param .data A SNPData object
#' @param ... Logical expressions to filter by, based on snp_info columns
#' @return A new filtered SNPData object
#' @export
#'
#' @examples
#' \dontrun{
#' snp_data <- get_example_snpdata()
#' # Filter SNPs with coverage > 10
#' filtered_snps <- filter_snps(snp_data, coverage > 10)
#' }
setGeneric("filter_snps", function(.data, ...) standardGeneric("filter_snps"))
#' @rdname filter_snps
setMethod(
  "filter_snps", signature(.data = "SNPData"),
  function(.data, ...) {
    # Capture NSE expressions
    dots <- rlang::enquos(...)

    # Check filter expressions
    check_filter_expr(.data@snp_info, dots, "snp_info")

    # Apply filter to snp_info
    selected_snps <- .data@snp_info %>%
      dplyr::filter(!!!dots)

    if (nrow(selected_snps) == 0) {
      stop("No SNPs remain after filtering. Please adjust your filter criteria.")
    }

    # Get indices of selected SNPs
    selected_indices <- match(selected_snps$snp_id, rownames(.data@alt_count))

    snp_total <- nrow(.data)
    snp_selected <- length(selected_indices)
    snp_removed <- snp_total - snp_selected
    removed_perc <- scales::percent(snp_removed / snp_total, accuracy = 0.01)
    logger::log_info("{snp_removed} ({removed_perc}) SNPs removed. {snp_selected} SNPs remaining.")
    # Subset the SNPData object
    .data[selected_indices, ]
  }
)

#' Filter SNPData object by sample/cell/barcode information
#'
#' @param .data A SNPData object
#' @param ... Logical expressions to filter by, based on barcode_info columns
#' @return A new filtered SNPData object
#' @export
#'
#' @examples
#' \dontrun{
#' # Filter cells with library size > 1000
#' filtered_cells <- filter_barcodes(snp_data, library_size > 1000)
#'
#' # Filter cells from a specific donor
#' filtered_cells <- filter_barcodes(snp_data, donor_id == "donor_1")
#'
#' # Filter cells with multiple conditions
#' filtered_cells <- filter_barcodes(snp_data, library_size > 1000, non_zero_snps > 50)
#' }
setGeneric("filter_barcodes", function(.data, ...) standardGeneric("filter_barcodes"))
#' @rdname filter_barcodes
setMethod(
  "filter_barcodes", signature(.data = "SNPData"),
  function(.data, ...) {
    # Capture NSE expressions
    dots <- rlang::enquos(...)

    # Check filter expressions
    check_filter_expr(.data@barcode_info, dots, "barcode_info")

    # Apply filter to barcode_info
    selected_samples <- .data@barcode_info %>%
      dplyr::filter(!!!dots)

    if (nrow(selected_samples) == 0) {
      stop("No barcodes remain after filtering. Please adjust your filter criteria.")
    }

    # Get indices of selected samples
    selected_indices <- match(selected_samples$cell_id, colnames(.data@alt_count))

    sample_total <- ncol(.data)
    sample_selected <- length(selected_indices)
    sample_removed <- sample_total - sample_selected
    removed_perc <- scales::percent(sample_removed / sample_total, accuracy = 0.01)
    logger::log_info("{sample_removed} ({removed_perc}) barcodes removed. {sample_selected} barcodes remaining.")

    # Subset the SNPData object
    .data[, selected_indices]
  }
)

#' @rdname filter_barcodes
#' @export
setGeneric("filter_samples", function(.data, ...) standardGeneric("filter_samples"))
#' @rdname filter_barcodes
setMethod(
  "filter_samples", signature(.data = "SNPData"),
  function(.data, ...) filter_barcodes(.data, ...)
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
#'
setGeneric("aggregate_count_df", function(x, group_by, test_maf = TRUE) standardGeneric("aggregate_count_df"))
#' @rdname aggregate_count_df
setMethod(
  "aggregate_count_df", signature(x = "SNPData"),
  function(x, group_by, test_maf = TRUE) {
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
      logger::log_warn("Found {sum(is.na(groups))} NA values in '{group_by}' column. These will be excluded from aggregation.")
      # Filter out samples with NA values
      keep_samples <- !is.na(groups)
      x_filtered <- x[, keep_samples]
      groups <- groups[keep_samples]
    } else {
      x_filtered <- x
    }

    logger::log_info("Extracting reference counts")
    ref_count_grouped <- groupedRowSums(ref_count(x_filtered), groups)
    ref_count_df <- ref_count_grouped %>%
      tibble::as_tibble() %>%
      dplyr::mutate(snp_id = rownames(ref_count_grouped), .before = 1) %>%
      tidyr::pivot_longer(-snp_id, names_to = group_by, values_to = "ref_count")

    logger::log_info("Extracting alternate counts")
    alt_count_grouped <- groupedRowSums(alt_count(x_filtered), groups)
    alt_count_df <- alt_count_grouped %>%
      tibble::as_tibble() %>%
      dplyr::mutate(snp_id = rownames(alt_count_grouped), .before = 1) %>%
      tidyr::pivot_longer(-snp_id, names_to = group_by, values_to = "alt_count")

    logger::log_info("Processing reference and alternate counts")
    out <- dplyr::inner_join(ref_count_df, alt_count_df, by = c("snp_id", group_by)) %>%
      dplyr::mutate(
        total_count = ref_count + alt_count,
        ref_ratio = ref_count / total_count,
        maf = pmin(ref_count, alt_count) / total_count
      ) %>%
      dplyr::filter(total_count > 0)
    if (isTRUE(test_maf)) {
      out <- out %>% dplyr::mutate(test_maf = !is.na(maf) & maf > 0)
    }

    logger::log_success("{stringr::str_to_title(group_by)} level counts calculated")
    out
  }
)

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

barcode_count_df_impl <- function(x, test_maf = TRUE) {
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

#' Get donor-level SNP heterozygosity status
#'
#' Returns a long-format data frame with heterozygosity status per SNP and donor.
#'
#' @param x A SNPData object
#' @param min_total_count Minimum total read depth (ref + alt) required per donor to test for heterozygosity (default: 10)
#' @param p_value_threshold P-value threshold for binomial test (default: 0.05). P-values are multiple-testing corrected and SNPs with p < threshold reject monoallelic expression.
#' @param minor_allele_prop Minor allele proportion used as the null threshold for monoallelic expression testing (default: 0.1).
#' @return A tibble with columns: snp_id, donor, ref_count, alt_count, total_count, minor_allele_count, p_val, adj_p_val, tested, status
#' @export
#'
#' @examples
#' \dontrun{
#' snp_data <- get_example_snpdata()
#' donor_het_status_df(snp_data)
#' }
#' @rdname donor_het_status_df
setGeneric("donor_het_status_df", function(
    x,
    min_total_count = 10,
    p_value_threshold = 0.05,
    minor_allele_prop = 0.1
) standardGeneric("donor_het_status_df"))

donor_het_status_df_impl <- function(
    x,
    min_total_count = 10,
    p_value_threshold = 0.05,
    minor_allele_prop = 0.1
) {
    stopifnot(min_total_count >= 1)
    stopifnot(p_value_threshold > 0 && p_value_threshold <= 1)
    stopifnot(minor_allele_prop > 0 && minor_allele_prop < 0.5)

    donor_counts <- donor_count_df(x, test_maf = FALSE) %>%
        dplyr::mutate(tested = total_count >= min_total_count)

    tested_counts <- donor_counts %>%
        dplyr::filter(tested)

    if (nrow(tested_counts) > 0) {
        tested_counts <- tested_counts %>%
            test_maf(p = minor_allele_prop) %>%
            dplyr::select(
                snp_id,
                donor,
                minor_allele_count,
                p_val,
                adj_p_val
            )
    } else {
        tested_counts <- donor_counts %>%
            dplyr::select(snp_id, donor) %>%
            dplyr::mutate(
                minor_allele_count = NA_real_,
                p_val = NA_real_,
                adj_p_val = NA_real_
            )
    }

    donor_counts %>%
        dplyr::left_join(tested_counts, by = c("snp_id", "donor")) %>%
        dplyr::mutate(
            status = dplyr::if_else(
                tested & !is.na(adj_p_val) & adj_p_val < p_value_threshold,
                "het",
                "hom"
            )
        )
}

#' @rdname donor_het_status_df
setMethod(
    "donor_het_status_df",
    signature(x = "SNPData"),
    donor_het_status_df_impl
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
    ref_count_df <- ref_count_grouped %>%
        tibble::as_tibble() %>%
        dplyr::mutate(snp_id = rownames(ref_count_grouped), .before = 1) %>%
        tidyr::pivot_longer(-snp_id, names_to = "clonotype", values_to = "ref_count")

    logger::log_info("Extracting alternate counts")
    alt_count_grouped <- groupedRowSums(alt_count(x), barcode_info$clonotype)
    alt_count_df <- alt_count_grouped %>%
        tibble::as_tibble() %>%
        dplyr::mutate(snp_id = rownames(alt_count_grouped), .before = 1) %>%
        tidyr::pivot_longer(-snp_id, names_to = "clonotype", values_to = "alt_count")

    logger::log_info("Processing reference and alternate counts")
    most_likely_donor <- barcode_info %>%
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

#' Check that columns in filter expressions exist in the data.frame
#'
#' Helper function to validate that all columns referenced in filter expressions
#' exist in the target data.frame or in the parent environment. Used internally
#' by filtering functions to provide informative error messages.
#'
#' @param df A data.frame to check column existence against
#' @param dots A list of quosures containing filter expressions
#' @param df_name Character string naming the data.frame for error messages (default "data.frame")
#' @return Invisibly returns NULL if all columns exist, otherwise throws an error
#' @keywords internal
check_filter_expr <- function(df, dots, df_name = "data.frame") {
    vars <- unique(unlist(lapply(dots, function(q) all.vars(rlang::get_expr(q)))))
    missing_vars <- setdiff(vars, colnames(df))

    # Check if variables exist in any accessible parent frame environments
    still_missing <- missing_vars
    if (length(still_missing) > 0) {
        # Start from parent.frame(2) to skip the calling filtering function
        # and check the environment where the user called the filtering function
        for (i in 2:sys.nframe()) {
            tryCatch(
                {
                    parent_frame <- parent.frame(i)

                    # Check which missing variables exist in this frame
                    found_vars <- still_missing[vapply(
                        still_missing,
                        exists,
                        logical(1),
                        envir = parent_frame,
                        inherits = TRUE # Allow inheritance to check enclosing environments
                    )]
                    still_missing <- setdiff(still_missing, found_vars)

                    # Stop if all variables found
                    if (length(still_missing) == 0) {
                        break
                    }
                },
                error = function(e) {
                    # Continue to next frame if this one is inaccessible
                }
            )
        }
    }

    if (length(still_missing) > 0) {
        stop(paste0(
            "The following columns are not present in ",
            df_name,
            " or accessible environment: ",
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
#' @rdname filter_snps
setGeneric("filter_snps", function(.data, ...) standardGeneric("filter_snps"))

filter_snps_impl <- function(.data, ...) {
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

#' @rdname filter_snps
setMethod(
    "filter_snps",
    signature(.data = "SNPData"),
    filter_snps_impl
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
#' snp_data <- get_example_snpdata()
#' # Filter cells with library size > 1000
#' filtered_cells <- filter_barcodes(snp_data, library_size > 1000)
#'
#' # Filter cells from a specific donor
#' filtered_cells <- filter_barcodes(snp_data, donor_id == "donor_1")
#'
#' # Filter cells with multiple conditions
#' filtered_cells <- filter_barcodes(snp_data, library_size > 1000, non_zero_snps > 50)
#' }
#' @rdname filter_barcodes
setGeneric("filter_barcodes", function(.data, ...) standardGeneric("filter_barcodes"))

filter_barcodes_impl <- function(.data, ...) {
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

#' @rdname filter_barcodes
setMethod(
    "filter_barcodes",
    signature(.data = "SNPData"),
    filter_barcodes_impl
)

#' @rdname filter_barcodes
#' @export
setGeneric("filter_samples", function(.data, ...) standardGeneric("filter_samples"))

filter_samples_impl <- function(.data, ...) {
    filter_barcodes(.data, ...)
}

#' @rdname filter_barcodes
setMethod(
    "filter_samples",
    signature(.data = "SNPData"),
    filter_samples_impl
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

#' Add metadata to SNPData barcode_info or snp_info
#'
#' These functions provide a standardized interface for adding new columns
#' to the barcode_info or snp_info data frames of a SNPData object. The functions
#' ensure data integrity by validating dimensions and preserving automatically
#' computed summary statistics.
#'
#' @param x A SNPData object
#' @param metadata A data.frame containing new columns to add to barcode_info or snp_info
#' @param join_by Character string specifying the column name to join by.
#'   For barcode_info: "cell_id" or "barcode"
#'   For snp_info: "snp_id"
#' @param overwrite Logical, whether to overwrite existing columns (default FALSE).
#'   Set to TRUE to update existing columns in barcode_info or snp_info.
#' @param validate Logical, whether to validate that all rows have matches (default TRUE).
#'   Set to FALSE when adding partial data (subset of barcodes or SNPs).
#'
#' @return A SNPData object with updated barcode_info or snp_info
#'
#' @examples
#' \dontrun{
#' # Add new columns to barcode_info
#' new_barcode_info <- data.frame(
#'   cell_id = c("cell_1", "cell_2"),
#'   cell_type = c("T_cell", "B_cell"),
#'   treatment = c("control", "treated")
#' )
#' updated_snpdata <- add_barcode_metadata(snpdata, new_barcode_info)
#'
#' # Update existing columns in barcode_info by setting overwrite=TRUE
#' clonotype_info <- data.frame(
#'   cell_id = c("cell_1", "cell_2"),
#'   clonotype = c("clonotype_1", "clonotype_2")
#' )
#' updated_snpdata <- add_barcode_metadata(
#'   snpdata,
#'   clonotype_info,
#'   overwrite = TRUE
#' )
#'
#' # Add columns to barcode_info for subset of cells (validate=FALSE)
#' partial_barcode_info <- data.frame(
#'   cell_id = c("cell_1", "cell_2"),
#'   annotation = c("A", "B")
#' )
#' updated_snpdata <- add_barcode_metadata(
#'   snpdata,
#'   partial_barcode_info,
#'   validate = FALSE
#' )
#'
#' # Add new columns to snp_info
#' new_snp_info <- data.frame(
#'   snp_id = c("snp_1", "snp_2"),
#'   gene_name = c("GENE1", "GENE2"),
#'   gene_biotype = c("protein_coding", "lncRNA")
#' )
#' updated_snpdata <- add_snp_metadata(snpdata, new_snp_info)
#' }
#'
#' @name add_metadata
NULL

#' @rdname add_metadata
#' @export
add_barcode_metadata <- function(x, metadata, join_by = "cell_id", overwrite = FALSE, validate = TRUE) {
    # Validate input
    if (!methods::is(x, "SNPData")) {
        stop("Input must be a SNPData object")
    }

    if (!is.data.frame(metadata)) {
        stop("metadata must be a data.frame")
    }

    if (!join_by %in% colnames(metadata)) {
        stop(paste0("Column '", join_by, "' not found in metadata data.frame"))
    }

    current_barcode_info <- get_barcode_info(x)

    # Check if join column exists in current barcode_info
    if (!join_by %in% colnames(current_barcode_info)) {
        stop(paste0("Column '", join_by, "' not found in current barcode_info"))
    }

    # Check for duplicate join keys in metadata
    if (any(duplicated(metadata[[join_by]]))) {
        stop(paste0("Duplicate values found in join column '", join_by, "' of metadata data.frame"))
    }

    # Check for column name conflicts in barcode_info
    if (!overwrite) {
        conflicting_cols <- intersect(
            setdiff(colnames(metadata), join_by),
            colnames(current_barcode_info)
        )
        if (length(conflicting_cols) > 0) {
            stop(paste0(
                "Column(s) already exist in barcode_info: ",
                paste(conflicting_cols, collapse = ", "),
                ". Set overwrite=TRUE to replace existing columns."
            ))
        }
    }

    # Validate that all barcodes in barcode_info have matches if requested
    if (validate) {
        missing_barcodes <- setdiff(current_barcode_info[[join_by]], metadata[[join_by]])
        if (length(missing_barcodes) > 0) {
            stop(paste0(
                "Some barcodes are missing from metadata data.frame: ",
                paste(head(missing_barcodes, 5), collapse = ", "),
                if (length(missing_barcodes) > 5) paste0(" (and ", length(missing_barcodes) - 5, " more)")
            ))
        }
    }

    # Perform the join
    updated_barcode_info <- current_barcode_info %>%
        dplyr::left_join(metadata, by = join_by, suffix = c("", ".new"))

    # Handle overwrite logic
    if (overwrite && any(grepl("\\.new$", colnames(updated_barcode_info)))) {
        new_cols <- grep("\\.new$", colnames(updated_barcode_info), value = TRUE)
        base_cols <- sub("\\.new$", "", new_cols)

        for (i in seq_along(new_cols)) {
            # Replace original column with new one
            updated_barcode_info[[base_cols[i]]] <- updated_barcode_info[[new_cols[i]]]
            # Remove the .new column
            updated_barcode_info[[new_cols[i]]] <- NULL
        }
    }

    # Preserve automatically computed columns in barcode_info
    # These are computed from the count matrices and should not be overwritten
    preserved_cols <- c("library_size", "non_zero_snps")
    auto_computed <- current_barcode_info[preserved_cols[preserved_cols %in% colnames(current_barcode_info)]]

    # Remove auto-computed columns from updated barcode_info if present
    for (col in preserved_cols) {
        if (col %in% colnames(updated_barcode_info)) {
            updated_barcode_info[[col]] <- NULL
        }
    }

    # Restore the original auto-computed columns
    for (col in names(auto_computed)) {
        updated_barcode_info[[col]] <- auto_computed[[col]]
    }

    # Create new SNPData object
    new(
        "SNPData",
        ref_count = x@ref_count,
        alt_count = x@alt_count,
        oth_count = x@oth_count,
        snp_info = x@snp_info,
        barcode_info = updated_barcode_info
    )
}

#' @rdname add_metadata
#' @export
add_snp_metadata <- function(x, metadata, join_by = "snp_id", overwrite = FALSE, validate = TRUE) {
    # Validate input
    if (!methods::is(x, "SNPData")) {
        stop("Input must be a SNPData object")
    }

    if (!is.data.frame(metadata)) {
        stop("metadata must be a data.frame")
    }

    if (!join_by %in% colnames(metadata)) {
        stop(paste0("Column '", join_by, "' not found in metadata data.frame"))
    }

    current_snp_info <- get_snp_info(x)

    # Check if join column exists in current snp_info
    if (!join_by %in% colnames(current_snp_info)) {
        stop(paste0("Column '", join_by, "' not found in current snp_info"))
    }

    # Check for duplicate join keys in metadata
    if (any(duplicated(metadata[[join_by]]))) {
        stop(paste0("Duplicate values found in join column '", join_by, "' of metadata data.frame"))
    }

    # Check for column name conflicts in snp_info
    if (!overwrite) {
        conflicting_cols <- intersect(
            setdiff(colnames(metadata), join_by),
            colnames(current_snp_info)
        )
        if (length(conflicting_cols) > 0) {
            stop(paste0(
                "Column(s) already exist in snp_info: ",
                paste(conflicting_cols, collapse = ", "),
                ". Set overwrite=TRUE to replace existing columns."
            ))
        }
    }

    # Validate that all SNPs in snp_info have matches if requested
    if (validate) {
        missing_snps <- setdiff(current_snp_info[[join_by]], metadata[[join_by]])
        if (length(missing_snps) > 0) {
            stop(paste0(
                "Some SNPs are missing from metadata data.frame: ",
                paste(head(missing_snps, 5), collapse = ", "),
                if (length(missing_snps) > 5) paste0(" (and ", length(missing_snps) - 5, " more)")
            ))
        }
    }

    # Perform the join
    updated_snp_info <- current_snp_info %>%
        dplyr::left_join(metadata, by = join_by, suffix = c("", ".new"))

    # Handle overwrite logic
    if (overwrite && any(grepl("\\.new$", colnames(updated_snp_info)))) {
        new_cols <- grep("\\.new$", colnames(updated_snp_info), value = TRUE)
        base_cols <- sub("\\.new$", "", new_cols)

        for (i in seq_along(new_cols)) {
            # Replace original column with new one
            updated_snp_info[[base_cols[i]]] <- updated_snp_info[[new_cols[i]]]
            # Remove the .new column
            updated_snp_info[[new_cols[i]]] <- NULL
        }
    }

    # Preserve automatically computed columns in snp_info
    # These are computed from the count matrices and should not be overwritten
    preserved_cols <- c("coverage", "non_zero_samples")
    auto_computed <- current_snp_info[preserved_cols[preserved_cols %in% colnames(current_snp_info)]]

    # Remove auto-computed columns from updated snp_info if present
    for (col in preserved_cols) {
        if (col %in% colnames(updated_snp_info)) {
            updated_snp_info[[col]] <- NULL
        }
    }

    # Restore the original auto-computed columns
    for (col in names(auto_computed)) {
        updated_snp_info[[col]] <- auto_computed[[col]]
    }

    # Create new SNPData object
    new(
        "SNPData",
        ref_count = x@ref_count,
        alt_count = x@alt_count,
        oth_count = x@oth_count,
        snp_info = updated_snp_info,
        barcode_info = x@barcode_info
    )
}

#' Expand and subset a sparse matrix to new dimensions
#'
#' Internal helper function to expand a matrix to include specified rows and columns,
#' filling missing entries with zeros. Used by merge_snpdata.
#'
#' @param mat A sparse Matrix
#' @param retained_rows Character vector of row names to retain
#' @param retained_cols Character vector of column names to retain
#' @param col_mapping Optional integer vector mapping original column indices to new positions
#' @return A sparse Matrix with dimensions matching retained_rows x retained_cols
#' @keywords internal
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

#' Merge SNP metadata according to join strategy
#'
#' @keywords internal
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

#' Merge barcode metadata according to join strategy
#'
#' @keywords internal
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

#' Merge two SNPData objects
#'
#' Combines two SNPData objects with flexible join strategies for SNPs and cells.
#' Overlapping entries are summed; unique entries are retained based on join type.
#'
#' @param x A SNPData object
#' @param y A SNPData object to merge with x
#' @param snp_join Join type for SNPs (rows). One of:
#'   \describe{
#'     \item{\code{"union"}}{Keep all SNPs from both objects (default)}
#'     \item{\code{"intersect"}}{Keep only SNPs present in both objects}
#'     \item{\code{"left"}}{Keep all SNPs from x, add y data where available}
#'     \item{\code{"right"}}{Keep all SNPs from y, add x data where available}
#'   }
#' @param cell_join Join type for cells (columns). One of:
#'   \describe{
#'     \item{\code{"union"}}{Keep all cells from both objects (default)}
#'     \item{\code{"intersect"}}{Keep only cells present in both objects}
#'     \item{\code{"left"}}{Keep all cells from x, add y data where available}
#'     \item{\code{"right"}}{Keep all cells from y, add x data where available}
#'   }
#'
#' @return A new SNPData object containing the merged data
#'
#' @details
#' This function provides independent control over which SNPs (rows) and cells
#' (columns) to retain when merging two SNPData objects. For positions present
#' in both objects, counts are summed element-wise. For positions present in
#' only one object (but retained by the join strategy), counts are kept with
#' zero-filling for the missing dimension.
#'
#' **Cell Merging Strategy:**
#' When both SNPData objects contain a \code{barcode} column in their
#' \code{barcode_info}, cells are matched and merged by their cell barcodes
#' rather than by \code{cell_id}. This allows proper merging of data from the
#' same physical cells across datasets, even if they have different internal
#' \code{cell_id} values. For cells with the same barcode, counts are summed.
#' If the \code{barcode} column is absent from either object, the function falls
#' back to merging by \code{cell_id}.
#'
#' Metadata is merged using the corresponding dplyr join function. For overlapping
#' SNPs or cells with conflicting metadata, values from x take priority. Auto-
#' computed columns (coverage, non_zero_samples, library_size, non_zero_snps)
#' are recalculated for the merged object.
#'
#' The independent join parameters enable fine-grained control:
#' \itemize{
#'   \item \code{snp_join="union", cell_join="union"}: Maximum data retention
#'   \item \code{snp_join="intersect", cell_join="intersect"}: Strict QC, only validated entries
#'   \item \code{snp_join="union", cell_join="intersect"}: Track more SNPs in same cells
#'   \item \code{snp_join="intersect", cell_join="union"}: Track consensus SNPs across more cells
#'   \item \code{snp_join="left", cell_join="left"}: Augment primary dataset
#' }
#'
#' @examples
#' \dontrun{
#' # Technical replicates - sum depth for validated entries only
#' rep1 <- import_cellsnp("replicate1/")
#' rep2 <- import_cellsnp("replicate2/")
#' combined <- merge_snpdata(rep1, rep2,
#'                          snp_join = "intersect",
#'                          cell_join = "intersect")
#'
#' # Batch integration - all cells, consensus SNPs only
#' batch1 <- import_cellsnp("donor1/")
#' batch2 <- import_cellsnp("donor2/")
#' integrated <- merge_snpdata(batch1, batch2,
#'                             snp_join = "intersect",
#'                             cell_join = "union")
#'
#' # SNP panel expansion - all SNPs, validated cells only
#' common <- import_cellsnp("common_vars/")
#' rare <- import_cellsnp("rare_vars/")
#' expanded <- merge_snpdata(common, rare,
#'                          snp_join = "union",
#'                          cell_join = "intersect")
#'
#' # Default: maximum data retention
#' dataset1 <- import_cellsnp("cohort_A/")
#' dataset2 <- import_cellsnp("cohort_B/")
#' combined <- merge_snpdata(dataset1, dataset2)
#' }
#'
#' @export
merge_snpdata <- function(
    x,
    y,
    snp_join = c("union", "intersect", "left", "right"),
    cell_join = c("union", "intersect", "left", "right")
) {
    # Match arguments
    snp_join <- match.arg(snp_join)
    cell_join <- match.arg(cell_join)

    # Validate both inputs are SNPData objects
    stopifnot(methods::is(x, "SNPData"))
    stopifnot(methods::is(y, "SNPData"))

    # Check for required identifier columns
    stopifnot("snp_id" %in% colnames(x@snp_info))
    stopifnot("snp_id" %in% colnames(y@snp_info))
    stopifnot("cell_id" %in% colnames(x@barcode_info))
    stopifnot("cell_id" %in% colnames(y@barcode_info))

    # Extract identifiers
    snp_ids_x <- rownames(x@ref_count)
    snp_ids_y <- rownames(y@ref_count)
    cell_ids_x <- colnames(x@ref_count)
    cell_ids_y <- colnames(y@ref_count)

    # Check if barcode column exists in both objects
    has_barcode_x <- "barcode" %in% colnames(x@barcode_info)
    has_barcode_y <- "barcode" %in% colnames(y@barcode_info)

    # Use barcodes for merging if available, otherwise fall back to cell_ids
    if (has_barcode_x && has_barcode_y) {
        # Extract barcodes from barcode_info
        barcodes_x <- x@barcode_info$barcode
        barcodes_y <- y@barcode_info$barcode

        # Create mapping from barcode to original cell_id (colnames)
        barcode_to_colname_x <- setNames(cell_ids_x, barcodes_x)
        barcode_to_colname_y <- setNames(cell_ids_y, barcodes_y)

        # Determine which barcodes to retain based on cell_join strategy
        barcodes_retained <- switch(
            cell_join,
            "union" = union(barcodes_x, barcodes_y),
            "intersect" = intersect(barcodes_x, barcodes_y),
            "left" = barcodes_x,
            "right" = barcodes_y
        )

        # Create mapping from original column names to retained barcode positions
        # This will be used to reorganize the matrices
        col_mapping_x <- match(barcodes_x, barcodes_retained)
        col_mapping_y <- match(barcodes_y, barcodes_retained)

        # Create new cell_ids for the retained barcodes
        cell_ids_retained <- paste0("cell_", seq_along(barcodes_retained))
    } else {
        # Fall back to using cell_ids from column names
        cell_ids_retained <- switch(
            cell_join,
            "union" = union(cell_ids_x, cell_ids_y),
            "intersect" = intersect(cell_ids_x, cell_ids_y),
            "left" = cell_ids_x,
            "right" = cell_ids_y
        )

        barcodes_retained <- cell_ids_retained
        col_mapping_x <- NULL
        col_mapping_y <- NULL
    }

    # Determine which SNPs to retain based on snp_join strategy
    snp_ids_retained <- switch(
        snp_join,
        "union" = union(snp_ids_x, snp_ids_y),
        "intersect" = intersect(snp_ids_x, snp_ids_y),
        "left" = snp_ids_x,
        "right" = snp_ids_y
    )

    # Handle edge case: empty result
    if (length(snp_ids_retained) == 0) {
        stop("No SNPs to retain after merge. Check your join strategy and input data.")
    }
    if (length(cell_ids_retained) == 0) {
        stop("No cells to retain after merge. Check your join strategy and input data.")
    }

    # Identify which entries need data from each object
    snp_from_x <- intersect(snp_ids_retained, snp_ids_x)
    snp_from_y <- intersect(snp_ids_retained, snp_ids_y)
    cell_from_x <- intersect(cell_ids_retained, cell_ids_x)
    cell_from_y <- intersect(cell_ids_retained, cell_ids_y)

    # Expand all matrices to retained dimensions
    ref_x_exp <- .expand_subset_matrix(x@ref_count, snp_ids_retained, cell_ids_retained, col_mapping_x)
    alt_x_exp <- .expand_subset_matrix(x@alt_count, snp_ids_retained, cell_ids_retained, col_mapping_x)
    oth_x_exp <- .expand_subset_matrix(x@oth_count, snp_ids_retained, cell_ids_retained, col_mapping_x)

    ref_y_exp <- .expand_subset_matrix(y@ref_count, snp_ids_retained, cell_ids_retained, col_mapping_y)
    alt_y_exp <- .expand_subset_matrix(y@alt_count, snp_ids_retained, cell_ids_retained, col_mapping_y)
    oth_y_exp <- .expand_subset_matrix(y@oth_count, snp_ids_retained, cell_ids_retained, col_mapping_y)

    # Merge by addition
    ref_merged <- ref_x_exp + ref_y_exp
    alt_merged <- alt_x_exp + alt_y_exp
    oth_merged <- oth_x_exp + oth_y_exp

    # Merge metadata using helper functions
    snp_info_merged <- .merge_snp_info(x, y, snp_ids_retained, snp_join)
    barcode_info_merged <- .merge_barcode_info(x, y, barcodes_retained, cell_join)

    # Create new SNPData object
    # The initialize method will automatically compute:
    # - coverage, non_zero_samples (in snp_info)
    # - library_size, non_zero_snps (in barcode_info)
    merged_obj <- SNPData(
        ref_count = ref_merged,
        alt_count = alt_merged,
        oth_count = oth_merged,
        snp_info = snp_info_merged,
        barcode_info = barcode_info_merged
    )

    logger::log_success(
        "Merged SNPData: {nrow(merged_obj)} SNPs x {ncol(merged_obj)} cells"
    )

    return(merged_obj)
}

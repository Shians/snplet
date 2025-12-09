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
#' @rdname barcode_count_df
setMethod(
    "barcode_count_df",
    signature(x = "SNPData"),
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
            out <- test_maf(out)
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
#' @rdname donor_count_df
setGeneric("donor_count_df", function(x, test_maf = TRUE) standardGeneric("donor_count_df"))
#' @rdname donor_count_df
setMethod(
    "donor_count_df",
    signature(x = "SNPData"),
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
            out <- test_maf(out)
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
#' @rdname clonotype_count_df
setGeneric("clonotype_count_df", function(x, test_maf = TRUE) standardGeneric("clonotype_count_df"))
#' @rdname clonotype_count_df
setMethod(
    "clonotype_count_df",
    signature(x = "SNPData"),
    function(x, test_maf = TRUE) {
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
            tidyr::pivot_longer(contains("clonotype"), names_to = "clonotype", values_to = "ref_count")

        logger::log_info("Extracting alternate counts")
        alt_count_grouped <- groupedRowSums(alt_count(x), barcode_info$clonotype)
        alt_count_df <- alt_count_grouped %>%
            tibble::as_tibble() %>%
            dplyr::mutate(snp_id = rownames(alt_count_grouped), .before = 1) %>%
            tidyr::pivot_longer(contains("clonotype"), names_to = "clonotype", values_to = "alt_count")

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
#' @rdname filter_snps
setMethod(
    "filter_snps",
    signature(.data = "SNPData"),
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
#' @rdname filter_barcodes
setMethod(
    "filter_barcodes",
    signature(.data = "SNPData"),
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
    "filter_samples",
    signature(.data = "SNPData"),
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
#' @rdname aggregate_count_df
setGeneric("aggregate_count_df", function(x, group_by, test_maf = TRUE) standardGeneric("aggregate_count_df"))
#' @rdname aggregate_count_df
setMethod(
    "aggregate_count_df",
    signature(x = "SNPData"),
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

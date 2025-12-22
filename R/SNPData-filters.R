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

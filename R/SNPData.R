setClass("SNPData",
    slots = c(
        alt_count = "Matrix",
        ref_count = "Matrix",
        snp_info = "data.frame",
        sample_info = "data.frame"
    )
)

setMethod("initialize", signature(.Object = "SNPData"),
    function(.Object, ref_count, alt_count, snp_info, sample_info) {
        # validate dimension compatibility
        stopifnot(nrow(alt_count) == nrow(ref_count))
        stopifnot(ncol(alt_count) == ncol(ref_count))

        stopifnot(ncol(alt_count) == nrow(sample_info))
        stopifnot(nrow(ref_count) == nrow(snp_info))

        # assign snp_ids if none exist
        if (!"snp_id" %in% colnames(snp_info)) {
            snp_info$snp_id <- paste0("snp_", seq_len(nrow(snp_info)))
        }

        # assign cell_ids if none exist
        if (!"cell_id" %in% colnames(sample_info)) {
            sample_info$cell_id <- paste0("cell_", seq_len(nrow(sample_info)))
        }

        colnames(ref_count) <- sample_info$cell_id
        colnames(alt_count) <- sample_info$cell_id

        rownames(ref_count) <- snp_info$snp_id
        rownames(alt_count) <- snp_info$snp_id

        .Object@ref_count <- ref_count
        .Object@alt_count <- alt_count
        .Object@snp_info <- snp_info
        .Object@sample_info <- sample_info

        .Object@snp_info$coverage <- Matrix::rowSums(alt_count + ref_count)
        .Object@snp_info$non_zero_samples <- Matrix::rowSums(alt_count + ref_count > 0)

        .Object@sample_info$library_size <- Matrix::colSums(alt_count + ref_count)
        .Object@sample_info$non_zero_snps <- Matrix::colSums(alt_count + ref_count > 0)

        .Object
    }
)

setMethod("[", signature(x = "SNPData", i = "ANY", j = "ANY"),
    function(x, i, j) {
        if (missing(i)) i <- seq_len(nrow(x@alt_count))
        if (missing(j)) j <- seq_len(ncol(x@alt_count))

        ref_count <- x@ref_count[i, j, drop = FALSE]
        alt_count <- x@alt_count[i, j, drop = FALSE]
        snp_info <- x@snp_info[i, ]
        sample_info <- x@sample_info[j, ]

        snp_info$coverage <- Matrix::rowSums(alt_count + ref_count)
        snp_info$non_zero_samples <- Matrix::rowSums(alt_count + ref_count > 0)

        sample_info$library_size <- Matrix::colSums(alt_count + ref_count)
        sample_info$non_zero_snps <- Matrix::colSums(alt_count + ref_count > 0)

        new("SNPData", alt_count = alt_count, ref_count = ref_count, snp_info = snp_info, sample_info = sample_info)
    }
)

# Constructor
setGeneric("SNPData", function(alt_count, ref_count, snp_info, sample_info) standardGeneric("SNPData"))
setMethod("SNPData", signature(alt_count = "Matrix", ref_count = "Matrix", snp_info = "data.frame", sample_info = "data.frame"),
    function(alt_count, ref_count, snp_info, sample_info) {
        new("SNPData", alt_count = alt_count, ref_count = ref_count, snp_info = snp_info, sample_info = sample_info)
    }
)

# Accessors
setGeneric("ref_count", function(object) standardGeneric("ref_count"))
setMethod("ref_count", signature(object = "SNPData"), function(object) object@ref_count)

setGeneric("alt_count", function(object) standardGeneric("alt_count"))
setMethod("alt_count", signature(object = "SNPData"), function(object) object@alt_count)

setGeneric("get_snp_info", function(object) standardGeneric("get_snp_info"))
setMethod("get_snp_info", signature(object = "SNPData"), function(object) object@snp_info)

setGeneric("get_sample_info", function(object) standardGeneric("get_sample_info"))
setMethod("get_sample_info", signature(object = "SNPData"), function(object) object@sample_info)

setMethod("nrow", signature(x = "SNPData"), function(x) nrow(x@ref_count))
setMethod("ncol", signature(x = "SNPData"), function(x) ncol(x@ref_count))

# Dimensions
setMethod("dim", signature(x = "SNPData"), function(x) c(nrow(x@alt_count), ncol(x@alt_count)))
setMethod("nrow", signature(x = "SNPData"), function(x) nrow(x@alt_count))
setMethod("ncol", signature(x = "SNPData"), function(x) ncol(x@alt_count))

setMethod("rownames", signature(x = "SNPData"), function(x) rownames(x@alt_count))
setMethod("colnames", signature(x = "SNPData"), function(x) colnames(x@alt_count))

# Show method
setMethod("show", signature(object = "SNPData"),
    function(object) {
        cat("Object of class 'SNPData'", "\n")
        cat("Dimensions: ", nrow(object), " SNPs x ", ncol(object), " samples", "\n")
        cat("SNP info:", "\n")
        print(object@snp_info)
        cat("Sample info:", "\n")
        print(object@sample_info)
    }
)

# Coverage method
setGeneric("coverage", function(x) standardGeneric("coverage"))
setMethod("coverage", signature(x = "SNPData"),
    function(x) {
        x@alt_count + x@ref_count
    }
)

# Reference fraction method
setGeneric("ref_fraction", function(x) standardGeneric("ref_fraction"))
setMethod("ref_fraction", signature(x = "SNPData"),
    function(x) {
        x@ref_count / (x@ref_count + x@alt_count)
    }
)

# Major allele fraction method
setGeneric("m_allele_frac", function(x) standardGeneric("m_allele_frac"))
setMethod("m_allele_frac", signature(x = "SNPData"),
    function(x) {
        abs(x@ref_count / (x@ref_count + x@alt_count) - 0.5) + 0.5
    }
)

setGeneric("cell_count_df", function(x) standardGeneric("cell_count_df"))
setMethod("cell_count_df", signature(x = "SNPData"),
    function(x) {
        ref_count_df <- x@ref_count %>%
            as.matrix() %>%
            tibble::as_tibble() %>%
            dplyr::mutate(snp_id = rownames(x@ref_count), .before = 1) %>%
            tidyr::pivot_longer(contains("cell"), names_to = "cell_id", values_to = "ref_count")

        alt_count_df <- x@alt_count %>%
            as.matrix() %>%
            tibble::as_tibble() %>%
            dplyr::mutate(snp_id = rownames(x@alt_count), .before = 1) %>%
            tidyr::pivot_longer(contains("cell"), names_to = "cell_id", values_to = "alt_count")

        dplyr::inner_join(ref_count_df, alt_count_df, by = c("snp_id", "cell_id")) %>%
            dplyr::mutate(
                total_count = ref_count + alt_count,
                ref_ratio = ref_count / total_count,
                maf = pmin(ref_count, alt_count) / total_count
            ) %>%
            dplyr::filter(total_count > 0)
    }
)

setGeneric("donor_count_df", function(x) standardGeneric("donor_count_df"))
setMethod("donor_count_df", signature(x = "SNPData"),
    function(x) {
        logger::log_info("Calculating donor level counts")

        logger::log_info("Extracting reference counts")
        ref_count_grouped <- groupedRowSums(ref_count(x), get_sample_info(x)$donor)
        ref_count_df <- ref_count_grouped %>%
            tibble::as_tibble() %>%
            dplyr::select(-any_of(c("unassigned", "doublet"))) %>%
            dplyr::mutate(snp_id = rownames(ref_count_grouped), .before = 1) %>%
            tidyr::pivot_longer(contains("donor"), names_to = "donor", values_to = "ref_count")

        logger::log_info("Extracting alternate counts")
        alt_count_grouped <- groupedRowSums(alt_count(x), get_sample_info(x)$donor)
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

        logger::log_success("Donor level counts calculated")
        out
    }
)

setGeneric("clonotype_count_df", function(x) standardGeneric("clonotype_count_df"))
setMethod("clonotype_count_df", signature(x = "SNPData"),
    function(x) {
        logger::log_info("Calculating clonotype level counts")

        logger::log_info("Extracting reference counts")
        ref_count_grouped <- groupedRowSums(ref_count(x), get_sample_info(x)$clonotype)
        ref_count_df <- ref_count_grouped %>%
            tibble::as_tibble() %>%
            dplyr::mutate(snp_id = rownames(ref_count_grouped), .before = 1) %>%
            tidyr::pivot_longer(contains("clonotype"), names_to = "clonotype", values_to = "ref_count")

        logger::log_info("Extracting alternate counts")
        alt_count_grouped <- groupedRowSums(alt_count(x), get_sample_info(x)$clonotype)
        alt_count_df <- alt_count_grouped %>%
            tibble::as_tibble() %>%
            dplyr::mutate(snp_id = rownames(alt_count_grouped), .before = 1) %>%
            tidyr::pivot_longer(contains("clonotype"), names_to = "clonotype", values_to = "alt_count")

        logger::log_info("Processing reference and alternate counts")
        most_likely_donor <- get_sample_info(x) %>%
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

        logger::log_success("Clonotype level counts calculated")
        out
    }
)

#' Filter SNPData object by SNP information
#'
#' @param .data A SNPData object
#' @param ... Logical expressions to filter by, based on snp_info columns
#' @return A new filtered SNPData object
#' @export
#'
#' @examples
#' \dontrun{
#' # Filter SNPs with coverage > 10
#' filtered_snps <- filter_snps(snp_data, coverage > 10)
#'
#' # Filter SNPs with multiple conditions
#' filtered_snps <- filter_snps(snp_data, coverage > 10, non_zero_samples > 5)
#' }
setGeneric("filter_snps", function(.data, ...) standardGeneric("filter_snps"))
setMethod("filter_snps", signature(.data = "SNPData"),
    function(.data, ...) {
        # Capture NSE expressions
        dots <- rlang::enquos(...)

        # Apply filter to snp_info
        selected_snps <- .data@snp_info %>%
            dplyr::filter(!!!dots)

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
#' @param ... Logical expressions to filter by, based on sample_info columns
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
setMethod("filter_barcodes", signature(.data = "SNPData"),
    function(.data, ...) {
        # Capture NSE expressions
        dots <- rlang::enquos(...)

        # Apply filter to sample_info
        selected_samples <- .data@sample_info %>%
            dplyr::filter(!!!dots)

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

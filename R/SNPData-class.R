#' SNPData: S4 class for single-cell SNP count data
#'
#' An S4 class to store and manipulate single-cell SNP count matrices with associated metadata.
#' This class integrates reference/alternate allele count matrices, SNP information, and
#' cell/barcode metadata to support allele-specific expression analysis workflows.
#'
#' @param ref_count A sparse Matrix containing reference allele counts (SNPs x cells)
#' @param alt_count A sparse Matrix containing alternate allele counts (SNPs x cells)
#' @param oth_count A sparse Matrix containing other allele counts (SNPs x cells), optional
#' @param snp_info A data.frame containing SNP metadata
#' @param barcode_info A data.frame containing cell/barcode metadata
#' @param object A SNPData object for show method
#' @param x A SNPData object
#' @param i Numeric or logical vector for subsetting SNPs (rows)
#' @param j Numeric or logical vector for subsetting samples (columns)
#' @param value A data.frame for replacement methods (barcode_info<- or snp_info<-)
#'
#' @slot ref_count A sparse Matrix containing reference allele counts (SNPs x cells)
#' @slot alt_count A sparse Matrix containing alternate allele counts (SNPs x cells)
#' @slot oth_count A sparse Matrix containing other allele counts (SNPs x cells)
#' @slot snp_info A data.frame containing SNP metadata with automatically computed coverage and non_zero_samples columns
#' @slot barcode_info A data.frame containing cell/barcode metadata with automatically computed library_size and non_zero_snps columns
#'
#' @section Accessors:
#' \describe{
#'   \item{\code{ref_count(x)}}{Get reference allele count matrix}
#'   \item{\code{alt_count(x)}}{Get alternate allele count matrix}
#'   \item{\code{oth_count(x)}}{Get other allele count matrix}
#'   \item{\code{get_snp_info(x)}}{Get SNP metadata data.frame}
#'   \item{\code{get_barcode_info(x)}}{Get cell/barcode metadata data.frame}
#'   \item{\code{get_sample_info(x)}}{Alias for get_barcode_info()}
#'   \item{\code{coverage(x)}}{Get total coverage matrix (ref + alt counts)}
#' }
#'
#' @section Dimensions:
#' \describe{
#'   \item{\code{nrow(x)}, \code{ncol(x)}, \code{dim(x)}}{Get object dimensions}
#'   \item{\code{rownames(x)}, \code{colnames(x)}}{Get row/column names (SNP/cell IDs)}
#' }
#'
#' @section Note:
#' The class automatically computes summary statistics in metadata:
#' \itemize{
#'   \item \code{snp_info$coverage}: Total counts per SNP across all cells
#'   \item \code{snp_info$non_zero_samples}: Number of cells with non-zero counts per SNP
#'   \item \code{barcode_info$library_size}: Total counts per cell across all SNPs
#'   \item \code{barcode_info$non_zero_snps}: Number of SNPs with non-zero counts per cell
#' }
#'
#' @examples
#' \dontrun{
#' # Create SNPData object
#' snp_data <- SNPData(
#'   ref_count = ref_matrix,
#'   alt_count = alt_matrix,
#'   snp_info = snp_df,
#'   barcode_info = cell_df
#' )
#'
#' # Access count matrices
#' ref_counts <- ref_count(snp_data)
#' alt_counts <- alt_count(snp_data)
#'
#' # Access metadata
#' snp_metadata <- get_snp_info(snp_data)
#' cell_metadata <- get_barcode_info(snp_data)
#'
#' # Calculate derived metrics
#' total_coverage <- coverage(snp_data)
#' }
#'
#' @exportClass SNPData
#' @include SNPData-helpers.R
#' @export
setClass(
    "SNPData",
    slots = c(
        ref_count = "Matrix",
        alt_count = "Matrix",
        oth_count = "Matrix",
        snp_info = "data.frame",
        barcode_info = "data.frame"
    )
)

setMethod(
    "initialize",
    signature(.Object = "SNPData"),
    function(
        .Object,
        ref_count,
        alt_count,
        oth_count = NULL,
        snp_info,
        barcode_info
    ) {
        oth_count <- .validate_count_dims(ref_count, alt_count, oth_count)
        .validate_info_dims(ref_count, alt_count, snp_info, barcode_info)

        snp_info <- .assign_snp_ids(snp_info)
        barcode_info <- .assign_cell_ids(barcode_info)

        deduped <- .dedupe_snps(ref_count, alt_count, oth_count, snp_info)
        ref_count <- deduped$ref_count
        alt_count <- deduped$alt_count
        oth_count <- deduped$oth_count
        snp_info <- deduped$snp_info

        # convert to tibble
        snp_info <- tibble::as_tibble(snp_info)
        barcode_info <- tibble::as_tibble(barcode_info)

        dimmed <- .set_dimnames(ref_count, alt_count, oth_count, snp_info, barcode_info)
        ref_count <- dimmed$ref_count
        alt_count <- dimmed$alt_count
        oth_count <- dimmed$oth_count

        .Object@ref_count <- ref_count
        .Object@alt_count <- alt_count
        .Object@oth_count <- oth_count
        metrics <- .recompute_metrics(snp_info, barcode_info, ref_count, alt_count)
        .Object@snp_info <- metrics$snp_info
        .Object@barcode_info <- metrics$barcode_info

        .Object
    }
)

#' Subset a SNPData object
#'
#' @param x A SNPData object
#' @param i Numeric or logical vector for subsetting SNPs (rows)
#' @param j Numeric or logical vector for subsetting samples (columns)
#' @return A subsetted SNPData object
#' @rdname SNPData-class
#' @export
setMethod(
    "[",
    signature(x = "SNPData", i = "ANY", j = "ANY"),
    function(x, i, j) {
        if (missing(i)) {
            i <- seq_len(nrow(x@alt_count))
        }
        if (missing(j)) {
            j <- seq_len(ncol(x@alt_count))
        }

        ref_count <- x@ref_count[i, j, drop = FALSE]
        alt_count <- x@alt_count[i, j, drop = FALSE]
        oth_count <- x@oth_count[i, j, drop = FALSE]
        snp_info <- x@snp_info[i, ]
        barcode_info <- x@barcode_info[j, ]

        new(
            "SNPData",
            alt_count = alt_count,
            ref_count = ref_count,
            oth_count = oth_count,
            snp_info = snp_info,
            barcode_info = barcode_info
        )
    }
)

# Constructor
#' @exportMethod SNPData
#' @rdname SNPData-class
setGeneric(
    "SNPData",
    function(ref_count, alt_count, snp_info, barcode_info, oth_count = NULL) {
        standardGeneric("SNPData")
    }
)
#' @exportMethod SNPData
#' @rdname SNPData-class
setMethod(
    "SNPData",
    signature(
        ref_count = "Matrix",
        alt_count = "Matrix",
        snp_info = "data.frame",
        barcode_info = "data.frame"
    ),
    function(ref_count, alt_count, snp_info, barcode_info, oth_count = NULL) {
        new(
            "SNPData",
            ref_count = ref_count,
            alt_count = alt_count,
            oth_count = oth_count,
            snp_info = snp_info,
            barcode_info = barcode_info
        )
    }
)

# Accessors
#' @exportMethod ref_count
#' @rdname SNPData-class
setGeneric("ref_count", function(x) standardGeneric("ref_count"))
#' @exportMethod ref_count
#' @rdname SNPData-class
setMethod("ref_count", signature(x = "SNPData"), function(x) x@ref_count)

#' @exportMethod alt_count
#' @rdname SNPData-class
setGeneric("alt_count", function(x) standardGeneric("alt_count"))
#' @exportMethod alt_count
#' @rdname SNPData-class
setMethod("alt_count", signature(x = "SNPData"), function(x) x@alt_count)

#' @exportMethod oth_count
#' @rdname SNPData-class
setGeneric("oth_count", function(x) standardGeneric("oth_count"))
#' @exportMethod oth_count
#' @rdname SNPData-class
setMethod("oth_count", signature(x = "SNPData"), function(x) x@oth_count)

#' @exportMethod get_snp_info
#' @rdname SNPData-class
setGeneric("get_snp_info", function(x) standardGeneric("get_snp_info"))
#' @exportMethod get_snp_info
#' @rdname SNPData-class
setMethod("get_snp_info", signature(x = "SNPData"), function(x) x@snp_info)

#' @exportMethod get_barcode_info
#' @rdname SNPData-class
setGeneric("get_barcode_info", function(x) standardGeneric("get_barcode_info"))
#' @exportMethod get_barcode_info
#' @rdname SNPData-class
setMethod("get_barcode_info", signature(x = "SNPData"), function(x) {
    x@barcode_info
})

#' @exportMethod get_sample_info
#' @rdname SNPData-class
setGeneric("get_sample_info", function(x) standardGeneric("get_sample_info"))
#' @exportMethod get_sample_info
#' @rdname SNPData-class
setMethod("get_sample_info", signature(x = "SNPData"), function(x) {
    get_barcode_info(x)
})

# Setters
#' @exportMethod barcode_info<-
#' @rdname SNPData-class
setGeneric("barcode_info<-", function(x, value) standardGeneric("barcode_info<-"))
#' @exportMethod barcode_info<-
#' @rdname SNPData-class
setReplaceMethod("barcode_info", signature(x = "SNPData", value = "data.frame"), function(x, value) {
    # Validate dimensions
    if (nrow(value) != ncol(x@ref_count)) {
        stop("Number of rows in barcode_info must match number of columns in count matrices")
    }
    # Validate required column exists
    if (!"cell_id" %in% colnames(value)) {
        stop("barcode_info must contain a 'cell_id' column")
    }
    # Validate cell_id matches column names of matrices
    if (!identical(value$cell_id, colnames(x@ref_count))) {
        stop("barcode_info$cell_id must match column names of count matrices")
    }
    x@barcode_info <- value
    x
})

#' @exportMethod snp_info<-
#' @rdname SNPData-class
setGeneric("snp_info<-", function(x, value) standardGeneric("snp_info<-"))
#' @exportMethod snp_info<-
#' @rdname SNPData-class
setReplaceMethod("snp_info", signature(x = "SNPData", value = "data.frame"), function(x, value) {
    # Validate dimensions
    if (nrow(value) != nrow(x@ref_count)) {
        stop("Number of rows in snp_info must match number of rows in count matrices")
    }
    # Validate required column exists
    if (!"snp_id" %in% colnames(value)) {
        stop("snp_info must contain a 'snp_id' column")
    }
    # Validate snp_id matches row names of matrices
    if (!identical(value$snp_id, rownames(x@ref_count))) {
        stop("snp_info$snp_id must match row names of count matrices")
    }
    x@snp_info <- value
    x
})

# Dimensions
#' Get dimensions of a SNPData object
#'
#' @param x A SNPData object
#' @return A numeric vector of length 2 giving the number of SNPs and samples
#' @rdname SNPData-class
#' @exportMethod dim
setMethod("dim", signature(x = "SNPData"), function(x) {
    c(nrow(x@alt_count), ncol(x@alt_count))
})
#' @rdname SNPData-class
#' @exportMethod nrow
setMethod("nrow", signature(x = "SNPData"), function(x) nrow(x@alt_count))
#' @rdname SNPData-class
#' @exportMethod ncol
setMethod("ncol", signature(x = "SNPData"), function(x) ncol(x@alt_count))

#' @exportMethod rownames
#' @rdname SNPData-class
setMethod("rownames", signature(x = "SNPData"), function(x) {
    rownames(x@alt_count)
})
#' @exportMethod colnames
#' @rdname SNPData-class
setMethod("colnames", signature(x = "SNPData"), function(x) {
    colnames(x@alt_count)
})

# Show method
#' @exportMethod show
#' @rdname SNPData-class
setMethod(
    "show",
    signature(object = "SNPData"),
    function(object) {
        cat("Object of class 'SNPData'", "\n")
        cat(
            "Dimensions: ",
            nrow(object),
            " SNPs x ",
            ncol(object),
            " samples",
            "\n"
        )
        cat("SNP info:", "\n")
        print(object@snp_info)
        cat("Sample info:", "\n")
        print(object@barcode_info)
    }
)

# Coverage method
#' @exportMethod coverage
#' @rdname SNPData-class
setGeneric("coverage", function(x) standardGeneric("coverage"))
#' @exportMethod coverage
#' @rdname SNPData-class
setMethod(
    "coverage",
    signature(x = "SNPData"),
    function(x) {
        x@alt_count + x@ref_count
    }
)

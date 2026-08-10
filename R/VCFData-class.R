#' VCFData: S4 class for VCF variant data
#'
#' An S4 class to store and manipulate VCF (Variant Call Format) data.
#' This class stores variant information, sample metadata, and VCF header information
#' to support genomic variant analysis workflows.
#'
#' @param header Character vector, required. VCF header lines.
#' @param samples Character vector, required. Sample names from the VCF file.
#' @param variants A data.frame, required, with standard VCF columns. Variant information.
#' @param object A VCFData object, required. Passed to the show method.
#' @param x A VCFData object, required.
#' @param i Numeric or logical vector, optional. Subsets variants (rows).
#' @param j Numeric or logical vector, optional. Subsets columns.
#'
#' @slot header A character vector containing VCF header lines
#' @slot samples A character vector of sample names from the VCF file
#' @slot variants A data.frame containing variant information with standard VCF columns
#'
#' @section Accessors:
#' \describe{
#'   \item{\code{header(x)}}{Get VCF header lines (alias: \code{get_header()})}
#'   \item{\code{samples(x)}}{Get sample names (alias: \code{get_samples()})}
#'   \item{\code{variants(x)}}{Get variant data.frame (alias: \code{get_variants()})}
#' }
#'
#' @section Dimensions:
#' \describe{
#'   \item{\code{nrow(x)}, \code{ncol(x)}, \code{dim(x)}}{Get object dimensions}
#'   \item{\code{rownames(x)}, \code{colnames(x)}}{Get row/column names}
#' }
#'
#' @examples
#' \dontrun{
#' # Create VCFData object
#' vcf_data <- VCFData(
#'   header = header_lines,
#'   samples = sample_names,
#'   variants = variant_df
#' )
#'
#' # Access components
#' header <- header(vcf_data)
#' samples <- samples(vcf_data)
#' variants <- variants(vcf_data)
#' }
#'
#' @exportClass VCFData
#' @export
setClass(
    "VCFData",
    slots = c(
        header = "character",
        samples = "character",
        variants = "data.frame"
    )
)

setMethod("initialize", signature(.Object = "VCFData"), function(.Object, header, samples, variants) {
    # Validate inputs
    stopifnot(is.character(header))
    stopifnot(is.character(samples))
    stopifnot(is.data.frame(variants))

    # Check that variants has required VCF columns
    required_cols <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")
    missing_cols <- setdiff(required_cols, colnames(variants))
    if (length(missing_cols) > 0) {
        stop("Missing required VCF columns: ", paste(missing_cols, collapse = ", "))
    }

    # Ensure POS is numeric
    if (!is.numeric(variants$POS)) {
        variants$POS <- as.numeric(variants$POS)
    }

    # Assign to slots
    .Object@header <- header
    .Object@samples <- samples
    .Object@variants <- variants

    .Object
})

# Constructor
#' @exportMethod VCFData
#' @rdname VCFData-class
setGeneric("VCFData", function(header, samples, variants) standardGeneric("VCFData"))
#' @exportMethod VCFData
#' @rdname VCFData-class
setMethod(
    "VCFData",
    signature(header = "character", samples = "character", variants = "data.frame"),
    function(header, samples, variants) {
        new("VCFData", header = header, samples = samples, variants = variants)
    }
)

# Accessors
#' @exportMethod header
#' @rdname VCFData-class
setGeneric("header", function(x) standardGeneric("header"))
#' @exportMethod header
#' @rdname VCFData-class
setMethod("header", signature(x = "VCFData"), function(x) x@header)

#' @export
#' @rdname VCFData-class
get_header <- function(x) header(x)

#' @exportMethod samples
#' @rdname VCFData-class
setGeneric("samples", function(x) standardGeneric("samples"))
#' @exportMethod samples
#' @rdname VCFData-class
setMethod("samples", signature(x = "VCFData"), function(x) x@samples)

#' @export
#' @rdname VCFData-class
get_samples <- function(x) samples(x)

#' @exportMethod variants
#' @rdname VCFData-class
setGeneric("variants", function(x) standardGeneric("variants"))
#' @exportMethod variants
#' @rdname VCFData-class
setMethod("variants", signature(x = "VCFData"), function(x) x@variants)

#' @export
#' @rdname VCFData-class
get_variants <- function(x) variants(x)

# Dimensions
#' @exportMethod nrow
#' @rdname VCFData-class
setMethod("nrow", signature(x = "VCFData"), function(x) nrow(x@variants))
#' @exportMethod ncol
#' @rdname VCFData-class
setMethod("ncol", signature(x = "VCFData"), function(x) ncol(x@variants))

#' Get dimensions of a VCFData object
#'
#' @param x A VCFData object, required.
#' @return A numeric vector of length 2 giving the number of variants and columns
#' @rdname VCFData-class
#' @exportMethod dim
setMethod("dim", signature(x = "VCFData"), function(x) c(nrow(x@variants), ncol(x@variants)))

#' @exportMethod rownames
#' @rdname VCFData-class
setMethod("rownames", signature(x = "VCFData"), function(x) rownames(x@variants))
#' @exportMethod colnames
#' @rdname VCFData-class
setMethod("colnames", signature(x = "VCFData"), function(x) colnames(x@variants))

# Show method
#' @exportMethod show
#' @rdname VCFData-class
setMethod("show", signature(object = "VCFData"), function(object) {
    cat("Object of class 'VCFData'", "\n")
    cat("Dimensions: ", nrow(object), " variants x ", ncol(object), " columns", "\n")
    cat("Samples: ", length(object@samples), " (", paste(head(object@samples, 3), collapse = ", "))
    if (length(object@samples) > 3) {
        cat("...)")
    } else {
        cat(")")
    }
    cat("\n")
    cat("Header lines: ", length(object@header), "\n")
    cat("Variants preview:", "\n")
    print(head(object@variants))
})

# Subset method
#' Subset a VCFData object
#'
#' @param x A VCFData object, required.
#' @param i Numeric or logical vector, optional. Subsets variants (rows).
#' @param j Numeric or logical vector, optional. Subsets columns.
#' @return A subsetted VCFData object
#' @rdname VCFData-class
#' @export
setMethod("[", signature(x = "VCFData", i = "ANY", j = "ANY"), function(x, i, j) {
    if (missing(i)) {
        i <- seq_len(nrow(x@variants))
    }
    if (missing(j)) {
        j <- seq_len(ncol(x@variants))
    }

    variants_subset <- x@variants[i, j, drop = FALSE]

    new("VCFData", header = x@header, samples = x@samples, variants = variants_subset)
})

#' Read VCF file and create VCFData object
#'
#' Reads a VCF (Variant Call Format) file and parses it into a VCFData object.
#' Handles header detection, sample extraction, and proper column naming.
#'
#' @param file Character scalar, required. Path to a VCF file (can be gzipped).
#' @return A VCFData object containing header, samples, and variant data
#' @export
#'
#' @examples
#' \dontrun{
#' vcf_data <- read_vcf("example.vcf")
#' vcf_data <- read_vcf("example.vcf.gz")
#' }
read_vcf <- function(file) {
    # Read header lines using "##" as comment to exclude #CHROM line
    if (grepl("\\.gz$", file)) {
        header_lines <- readLines(gzfile(file))
    } else {
        header_lines <- readLines(file)
    }

    # Extract only the ## header lines (not #CHROM)
    header_lines <- header_lines[grepl("^##", header_lines)]

    # Read the variant data using read_tsv with "##" comment
    variants <- readr::read_tsv(
        file,
        comment = "##",
        col_types = readr::cols(
            .default = readr::col_character(),
            POS = readr::col_integer()
        ),
        show_col_types = FALSE
    )

    # Clean up column names (remove # from CHROM)
    colnames(variants)[1] <- gsub("^#", "", colnames(variants)[1])

    # Extract sample names (columns after FORMAT, if present)
    col_names <- colnames(variants)
    format_idx <- which(col_names == "FORMAT")
    if (length(format_idx) > 0) {
        sample_names <- col_names[(format_idx + 1):length(col_names)]
    } else {
        # No FORMAT column, check if there are sample columns after INFO
        info_idx <- which(col_names == "INFO")
        if (length(info_idx) > 0 && length(col_names) > info_idx) {
            sample_names <- col_names[(info_idx + 1):length(col_names)]
        } else {
            sample_names <- character(0)
        }
    }

    # Create VCFData object
    vcf_data <- VCFData(
        header = header_lines,
        samples = sample_names,
        variants = variants
    )

    return(vcf_data)
}

#' VCFData: S4 class for VCF variant data
#'
#' An S4 class to store and manipulate VCF (Variant Call Format) data.
#' This class stores variant information, sample metadata, and VCF header information
#' to support genomic variant analysis workflows.
#'
#' @slot header A character vector containing VCF header lines
#' @slot samples A character vector of sample names from the VCF file
#' @slot variants A data.frame containing variant information with standard VCF columns
#'
#' @section Accessors:
#' \describe{
#'   \item{\code{get_header(x)}}{Get VCF header lines}
#'   \item{\code{get_samples(x)}}{Get sample names}
#'   \item{\code{get_variants(x)}}{Get variant data.frame}
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
#' header <- get_header(vcf_data)
#' samples <- get_samples(vcf_data)
#' variants <- get_variants(vcf_data)
#' }
#'
#' @exportClass VCFData
#' @export
setClass("VCFData",
    slots = c(
        header = "character",
        samples = "character",
        variants = "data.frame"
    )
)

setMethod("initialize", signature(.Object = "VCFData"),
    function(.Object, header, samples, variants) {
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
    }
)

# Constructor
#' @exportMethod VCFData
setGeneric("VCFData", function(header, samples, variants) standardGeneric("VCFData"))
#' @exportMethod VCFData
setMethod("VCFData", signature(header = "character", samples = "character", variants = "data.frame"),
    function(header, samples, variants) {
        new("VCFData", header = header, samples = samples, variants = variants)
    }
)

# Accessors
#' @exportMethod get_header
setGeneric("get_header", function(x) standardGeneric("get_header"))
#' @exportMethod get_header
setMethod("get_header", signature(x = "VCFData"), function(x) x@header)

#' @exportMethod get_samples
setGeneric("get_samples", function(x) standardGeneric("get_samples"))
#' @exportMethod get_samples
setMethod("get_samples", signature(x = "VCFData"), function(x) x@samples)

#' @exportMethod get_variants
setGeneric("get_variants", function(x) standardGeneric("get_variants"))
#' @exportMethod get_variants
setMethod("get_variants", signature(x = "VCFData"), function(x) x@variants)

# Dimensions
#' @exportMethod nrow
setMethod("nrow", signature(x = "VCFData"), function(x) nrow(x@variants))
#' @exportMethod ncol
setMethod("ncol", signature(x = "VCFData"), function(x) ncol(x@variants))

#' Get dimensions of a VCFData object
#'
#' @param x A VCFData object
#' @return A numeric vector of length 2 giving the number of variants and columns
#' @rdname VCFData-class
#' @exportMethod dim
setMethod("dim", signature(x = "VCFData"), function(x) c(nrow(x@variants), ncol(x@variants)))

#' @exportMethod rownames
setMethod("rownames", signature(x = "VCFData"), function(x) rownames(x@variants))
#' @exportMethod colnames
setMethod("colnames", signature(x = "VCFData"), function(x) colnames(x@variants))

# Show method
#' @exportMethod show
setMethod("show", signature(object = "VCFData"),
    function(object) {
        cat("Object of class 'VCFData'", "\n")
        cat("Dimensions: ", nrow(object), " variants x ", ncol(object), " columns", "\n")
        cat("Samples: ", length(object@samples), " (",
            paste(head(object@samples, 3), collapse = ", "),
            if (length(object@samples) > 3) "..." else "", ")", "\n")
        cat("Header lines: ", length(object@header), "\n")
        cat("Variants preview:", "\n")
        print(head(object@variants))
    }
)

# Subset method
#' Subset a VCFData object
#'
#' @param x A VCFData object
#' @param i Numeric or logical vector for subsetting variants (rows)
#' @param j Numeric or logical vector for subsetting columns
#' @return A subsetted VCFData object
#' @rdname VCFData-class
#' @export
setMethod("[", signature(x = "VCFData", i = "ANY", j = "ANY"),
    function(x, i, j) {
        if (missing(i)) i <- seq_len(nrow(x@variants))
        if (missing(j)) j <- seq_len(ncol(x@variants))

        variants_subset <- x@variants[i, j, drop = FALSE]

        new("VCFData",
            header = x@header,
            samples = x@samples,
            variants = variants_subset)
    }
)

#' Read VCF file and create VCFData object
#'
#' Reads a VCF (Variant Call Format) file and parses it into a VCFData object.
#' Handles header detection, sample extraction, and proper column naming.
#'
#' @param file Path to VCF file (can be gzipped)
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

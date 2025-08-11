setClass("SNPData",
    slots = c(
        ref_count = "Matrix",
        alt_count = "Matrix",
        oth_count = "Matrix",
        snp_info = "data.frame",
        barcode_info = "data.frame"
    )
)

setMethod("initialize", signature(.Object = "SNPData"),
    function(.Object, ref_count, alt_count, oth_count = NULL, snp_info, barcode_info) {
        # validate dimension compatibility
        stopifnot(nrow(alt_count) == nrow(ref_count))
        stopifnot(ncol(alt_count) == ncol(ref_count))
        # OTH matrix checks
        if (!is.null(oth_count)) {
            stopifnot(nrow(oth_count) == nrow(ref_count))
            stopifnot(ncol(oth_count) == ncol(ref_count))
        } else {
            oth_count <- Matrix::Matrix(0, nrow = nrow(ref_count), ncol = ncol(ref_count), sparse = TRUE)
        }

        stopifnot(ncol(alt_count) == nrow(barcode_info))
        stopifnot(nrow(ref_count) == nrow(snp_info))

        # assign snp_ids if none exist
        if (!"snp_id" %in% colnames(snp_info)) {
            snp_info$snp_id <- paste0("snp_", seq_len(nrow(snp_info)))
        }

        # assign cell_ids if none exist
        if (!"cell_id" %in% colnames(barcode_info)) {
            barcode_info$cell_id <- paste0("cell_", seq_len(nrow(barcode_info)))
        }

        colnames(ref_count) <- barcode_info$cell_id
        colnames(alt_count) <- barcode_info$cell_id
        colnames(oth_count) <- barcode_info$cell_id
        rownames(ref_count) <- snp_info$snp_id
        rownames(alt_count) <- snp_info$snp_id
        rownames(oth_count) <- snp_info$snp_id

        .Object@ref_count <- ref_count
        .Object@alt_count <- alt_count
        .Object@oth_count <- oth_count
        .Object@snp_info <- snp_info
        .Object@barcode_info <- barcode_info

        .Object@snp_info$coverage <- Matrix::rowSums(alt_count + ref_count)
        .Object@snp_info$non_zero_samples <- Matrix::rowSums(alt_count + ref_count > 0)

        .Object@barcode_info$library_size <- Matrix::colSums(alt_count + ref_count)
        .Object@barcode_info$non_zero_snps <- Matrix::colSums(alt_count + ref_count > 0)

        .Object
    }
)

setMethod("[", signature(x = "SNPData", i = "ANY", j = "ANY"),
    function(x, i, j) {
        if (missing(i)) i <- seq_len(nrow(x@alt_count))
        if (missing(j)) j <- seq_len(ncol(x@alt_count))

        ref_count <- x@ref_count[i, j, drop = FALSE]
        alt_count <- x@alt_count[i, j, drop = FALSE]
        oth_count <- x@oth_count[i, j, drop = FALSE]
        snp_info <- x@snp_info[i, ]
        barcode_info <- x@barcode_info[j, ]

        snp_info$coverage <- Matrix::rowSums(alt_count + ref_count)
        snp_info$non_zero_samples <- Matrix::rowSums(alt_count + ref_count > 0)

        barcode_info$library_size <- Matrix::colSums(alt_count + ref_count)
        barcode_info$non_zero_snps <- Matrix::colSums(alt_count + ref_count > 0)

        new("SNPData", alt_count = alt_count, ref_count = ref_count, oth_count = oth_count, snp_info = snp_info, barcode_info = barcode_info)
    }
)

# Constructor
setGeneric("SNPData", function(ref_count, alt_count, snp_info, barcode_info, oth_count = NULL) standardGeneric("SNPData"))
setMethod("SNPData", signature(ref_count = "Matrix", alt_count = "Matrix", snp_info = "data.frame", barcode_info = "data.frame"),
    function(ref_count, alt_count, snp_info, barcode_info, oth_count = NULL) {
        new("SNPData", ref_count = ref_count, alt_count = alt_count, oth_count = oth_count, snp_info = snp_info, barcode_info = barcode_info)
    }
)

# Accessors
setGeneric("ref_count", function(object) standardGeneric("ref_count"))
setMethod("ref_count", signature(object = "SNPData"), function(object) object@ref_count)

setGeneric("alt_count", function(object) standardGeneric("alt_count"))
setMethod("alt_count", signature(object = "SNPData"), function(object) object@alt_count)

setGeneric("oth_count", function(object) standardGeneric("oth_count"))
setMethod("oth_count", signature(object = "SNPData"), function(object) object@oth_count)

setGeneric("get_snp_info", function(object) standardGeneric("get_snp_info"))
setMethod("get_snp_info", signature(object = "SNPData"), function(object) object@snp_info)

setGeneric("get_barcode_info", function(object) standardGeneric("get_barcode_info"))
setMethod("get_barcode_info", signature(object = "SNPData"), function(object) object@barcode_info)

setGeneric("get_sample_info", function(object) standardGeneric("get_sample_info"))
setMethod("get_sample_info", signature(object = "SNPData"), function(object) get_barcode_info(object))

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
        print(object@barcode_info)
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
setGeneric("major_allele_frac", function(x) standardGeneric("major_allele_frac"))
setMethod("major_allele_frac", signature(x = "SNPData"),
    function(x) {
        abs(x@ref_count / (x@ref_count + x@alt_count) - 0.5) + 0.5
    }
)

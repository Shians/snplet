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
#' @param donor_info A data.frame with one row per donor, optional. Defaults to one row per
#'   distinct value of \code{barcode_info$donor} (empty if \code{barcode_info} has no
#'   \code{donor} column). Must contain a \code{donor} column; \code{n_cells} is computed
#'   automatically.
#' @param donor_snp_info A data.frame with one row per (\code{snp_id}, \code{donor},
#'   \code{zygosity_source}) triple, optional -- a pair may carry more than one row when
#'   more than one source has called it (e.g. a \code{"vireo_gt"} row and a
#'   \code{"binomial"} row for the same pair). Defaults to an empty table. Must contain
#'   \code{snp_id} and \code{donor} columns, both of which must already appear in
#'   \code{snp_info}/\code{donor_info}; a row with a non-\code{NA} \code{zygosity} must
#'   also carry a non-\code{NA} \code{zygosity_source}.
#' @param source Character vector of \code{zygosity_source} value(s) to filter
#'   \code{donor_snp_info(x)} to, or \code{"all"} to return every source unfiltered.
#'   Defaults to \code{NULL}, meaning the object's active source (see
#'   \code{\link{zygosity_source}}).
#' @param donor_map A named character vector, \code{c(new_label = old_label, ...)} (the same
#'   \code{new = old} convention as \code{dplyr::rename()}), optional. Applied to
#'   \code{barcode_info$donor} and, if supplied, \code{donor_info$donor}/
#'   \code{donor_snp_info$donor} before the object is built -- useful since Vireo assigns
#'   arbitrary labels (\code{donor0}, \code{donor1}, ...) that this can relabel at import
#'   time rather than after the fact via \code{\link{rename_donor}}.
#' @param object A SNPData object for show method
#' @param x A SNPData object
#' @param i Numeric or logical vector for subsetting SNPs (rows)
#' @param j Numeric or logical vector for subsetting samples (columns)
#' @param value A data.frame for replacement methods (barcode_info<- or snp_info<-)
#' @param ... Additional arguments (unused; required by the \code{updateObject} generic)
#' @param verbose Logical, whether \code{updateObject} logs which legacy slots it migrated
#'   (default \code{FALSE})
#'
#' @slot ref_count A sparse Matrix containing reference allele counts (SNPs x cells)
#' @slot alt_count A sparse Matrix containing alternate allele counts (SNPs x cells)
#' @slot oth_count A sparse Matrix containing other allele counts (SNPs x cells)
#' @slot snp_info A data.frame containing SNP metadata with automatically computed coverage and non_zero_samples columns
#' @slot barcode_info A data.frame containing cell/barcode metadata with automatically computed library_size and non_zero_snps columns
#' @slot chr_style Character string indicating the chromosome naming style. One of: "numeric", "ucsc", "refseq_mouse", "genbank_mouse", "refseq_human", "genbank_human", or "unknown"
#' @slot donor_info A tibble with one row per donor and an automatically computed
#'   \code{n_cells} column. Rows are dropped when a donor loses all of its cells, via
#'   subsetting or \code{merge_snpdata()}.
#' @slot donor_snp_info A tibble with one row per (\code{snp_id}, \code{donor},
#'   \code{zygosity_source}) triple, carrying zygosity calls (\code{zygosity},
#'   \code{zygosity_source} and per-source confidence columns) and XCI fit diagnostics
#'   (\code{xci_informative}, \code{allele_on_x1}, \code{xci_escape_fraction}) written by
#'   \code{\link{assign_xci}} for whichever source was active when it ran. Rows are dropped
#'   along with their donor, same as \code{donor_info}.
#' @slot zygosity_source Character string naming the \emph{active} zygosity-call source
#'   (a value of \code{donor_snp_info$zygosity_source}), or \code{NA_character_} if none
#'   has been established yet. \code{\link{donor_snp_info}}, \code{\link{assign_xci}}, and
#'   other zygosity-dependent functions default to this source; switch it with
#'   \code{zygosity_source<-()}. Set automatically at import (e.g. \code{"vireo_gt"} after
#'   importing Vireo genotypes) or by \code{\link{infer_zygosity}}.
#'
#' @section Accessors:
#' \describe{
#'   \item{\code{ref_count(x)}}{Get reference allele count matrix}
#'   \item{\code{alt_count(x)}}{Get alternate allele count matrix}
#'   \item{\code{oth_count(x)}}{Get other allele count matrix}
#'   \item{\code{snp_info(x)}}{Get SNP metadata data.frame (alias: \code{get_snp_info()})}
#'   \item{\code{barcode_info(x)}}{Get cell/barcode metadata data.frame (alias: \code{get_barcode_info()})}
#'   \item{\code{get_sample_info(x)}}{Alias for barcode_info()}
#'   \item{\code{donor_info(x)}}{Get per-donor metadata tibble (alias: \code{get_donor_info()})}
#'   \item{\code{donor_snp_info(x, source = NULL)}}{Get per-(SNP, donor) metadata tibble,
#'     filtered to the active zygosity source by default (alias: \code{get_donor_snp_info()})}
#'   \item{\code{chr_style(x)}}{Get chromosome naming style}
#'   \item{\code{zygosity_source(x)}}{Get the active zygosity-call source; set with
#'     \code{zygosity_source<-()}}
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
#' snp_metadata <- snp_info(snp_data)
#' cell_metadata <- barcode_info(snp_data)
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
        barcode_info = "data.frame",
        chr_style = "character",
        donor_info = "tbl_df",
        donor_snp_info = "tbl_df",
        zygosity_source = "character"
    )
)

setValidity("SNPData", function(object) {
    if (!methods::.hasSlot(object, "donor_snp_info")) {
        return(TRUE)
    }
    problems <- .validate_donor_dims(object@donor_info, object@donor_snp_info, object@snp_info)
    if (length(problems) == 0) TRUE else problems
})

setMethod(
    "initialize",
    signature(.Object = "SNPData"),
    function(
        .Object,
        ref_count,
        alt_count,
        oth_count = NULL,
        snp_info,
        barcode_info,
        donor_info = NULL,
        donor_snp_info = NULL,
        donor_map = NULL
    ) {
        oth_count <- .validate_count_dims(ref_count, alt_count, oth_count)
        .validate_info_dims(ref_count, alt_count, snp_info, barcode_info)

        snp_info <- .assign_snp_ids(snp_info)
        barcode_info <- .assign_cell_ids(barcode_info)

        # Relabel donors (e.g. Vireo's arbitrary donor0..n) before donor_info
        # is derived, so barcode_info, donor_info, and any caller-supplied
        # donor_snp_info (e.g. from a Vireo GT VCF) all key on the same labels.
        if (!is.null(donor_map)) {
            if ("donor" %in% colnames(barcode_info)) {
                barcode_info$donor <- .apply_donor_map(barcode_info$donor, donor_map)
            }
            if (!is.null(donor_info) && "donor" %in% colnames(donor_info)) {
                donor_info$donor <- .apply_donor_map(donor_info$donor, donor_map)
            }
            if (!is.null(donor_snp_info) && "donor" %in% colnames(donor_snp_info)) {
                donor_snp_info$donor <- .apply_donor_map(donor_snp_info$donor, donor_map)
            }
        }

        donor_info <- donor_info %||% .default_donor_info(barcode_info)
        donor_snp_info <- donor_snp_info %||% .empty_donor_snp_info()

        deduped <- .dedupe_snps(ref_count, alt_count, oth_count, snp_info, donor_snp_info)
        ref_count <- deduped$ref_count
        alt_count <- deduped$alt_count
        oth_count <- deduped$oth_count
        snp_info <- deduped$snp_info
        donor_snp_info <- deduped$donor_snp_info

        # convert to tibble
        snp_info <- tibble::as_tibble(snp_info)
        barcode_info <- tibble::as_tibble(barcode_info)
        donor_info <- tibble::as_tibble(donor_info)
        donor_snp_info <- tibble::as_tibble(donor_snp_info)

        # Detect chromosome style and add canonical chromosome column
        if ("chrom" %in% colnames(snp_info)) {
            chr_style <- detect_chr_style(snp_info$chrom)
            snp_info$chrom_canonical <- normalise_chr_names(snp_info$chrom, from_style = chr_style)
        } else {
            chr_style <- "unknown"
        }

        dimmed <- .set_dimnames(ref_count, alt_count, oth_count, snp_info, barcode_info)
        ref_count <- dimmed$ref_count
        alt_count <- dimmed$alt_count
        oth_count <- dimmed$oth_count

        .Object@ref_count <- ref_count
        .Object@alt_count <- alt_count
        .Object@oth_count <- oth_count
        .Object@chr_style <- chr_style
        metrics <- .recompute_metrics(snp_info, barcode_info, donor_info, ref_count, alt_count)
        .Object@snp_info <- metrics$snp_info
        .Object@barcode_info <- metrics$barcode_info
        .Object@donor_info <- metrics$donor_info
        .Object@donor_snp_info <- donor_snp_info
        .Object@zygosity_source <- .derive_zygosity_source(donor_snp_info)

        methods::validObject(.Object)
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

        # Donor genotypes are a property of the donor, not of the retained
        # cells: a donor with no cells left after subsetting has its rows
        # dropped from both donor tables rather than carried over stale.
        if (methods::.hasSlot(x, "donor_snp_info")) {
            surviving_donors <- if ("donor" %in% colnames(barcode_info)) {
                unique(stats::na.omit(barcode_info$donor))
            } else {
                character(0)
            }
            donor_info <- x@donor_info[x@donor_info$donor %in% surviving_donors, , drop = FALSE]
            donor_snp_info <- x@donor_snp_info[
                x@donor_snp_info$snp_id %in% snp_info$snp_id & x@donor_snp_info$donor %in% surviving_donors,
                ,
                drop = FALSE
            ]
        } else {
            donor_info <- NULL
            donor_snp_info <- NULL
        }

        obj <- new(
            "SNPData",
            alt_count = alt_count,
            ref_count = ref_count,
            oth_count = oth_count,
            snp_info = snp_info,
            barcode_info = barcode_info,
            donor_info = donor_info,
            donor_snp_info = donor_snp_info
        )
        # Handle backwards compatibility with older SNPData objects
        if (methods::.hasSlot(x, "chr_style")) {
            obj@chr_style <- x@chr_style
        }
        obj <- .propagate_zygosity_source(obj, x)
        obj
    }
)

# Constructor
#' @exportMethod SNPData
#' @rdname SNPData-class
setGeneric(
    "SNPData",
    function(
        ref_count,
        alt_count,
        snp_info,
        barcode_info,
        oth_count = NULL,
        donor_info = NULL,
        donor_snp_info = NULL,
        donor_map = NULL
    ) {
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
    function(
        ref_count,
        alt_count,
        snp_info,
        barcode_info,
        oth_count = NULL,
        donor_info = NULL,
        donor_snp_info = NULL,
        donor_map = NULL
    ) {
        new(
            "SNPData",
            ref_count = ref_count,
            alt_count = alt_count,
            oth_count = oth_count,
            snp_info = snp_info,
            barcode_info = barcode_info,
            donor_info = donor_info,
            donor_snp_info = donor_snp_info,
            donor_map = donor_map
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

#' @exportMethod snp_info
#' @rdname SNPData-class
setGeneric("snp_info", function(x) standardGeneric("snp_info"))
#' @exportMethod snp_info
#' @rdname SNPData-class
setMethod("snp_info", signature(x = "SNPData"), function(x) x@snp_info)

#' @export
#' @rdname SNPData-class
get_snp_info <- function(x) snp_info(x)

#' @exportMethod barcode_info
#' @rdname SNPData-class
setGeneric("barcode_info", function(x) standardGeneric("barcode_info"))
#' @exportMethod barcode_info
#' @rdname SNPData-class
setMethod("barcode_info", signature(x = "SNPData"), function(x) {
    x@barcode_info
})

#' @export
#' @rdname SNPData-class
get_barcode_info <- function(x) barcode_info(x)

#' @exportMethod get_sample_info
#' @rdname SNPData-class
setGeneric("get_sample_info", function(x) standardGeneric("get_sample_info"))
#' @exportMethod get_sample_info
#' @rdname SNPData-class
setMethod("get_sample_info", signature(x = "SNPData"), function(x) {
    barcode_info(x)
})

#' @exportMethod donor_info
#' @rdname SNPData-class
setGeneric("donor_info", function(x) standardGeneric("donor_info"))
#' @exportMethod donor_info
#' @rdname SNPData-class
setMethod("donor_info", signature(x = "SNPData"), function(x) {
    # Handle backwards compatibility with older SNPData objects
    if (!methods::.hasSlot(x, "donor_info")) {
        return(.default_donor_info(x@barcode_info))
    }
    x@donor_info
})

#' @export
#' @rdname SNPData-class
get_donor_info <- function(x) donor_info(x)

#' @exportMethod donor_snp_info
#' @rdname SNPData-class
setGeneric("donor_snp_info", function(x, source = NULL) standardGeneric("donor_snp_info"))
#' @exportMethod donor_snp_info
#' @rdname SNPData-class
setMethod("donor_snp_info", signature(x = "SNPData"), function(x, source = NULL) {
    # Handle backwards compatibility with older SNPData objects
    if (!methods::.hasSlot(x, "donor_snp_info")) {
        return(.empty_donor_snp_info())
    }
    raw <- x@donor_snp_info
    if (identical(source, "all")) {
        return(raw)
    }
    active_source <- if (is.null(source)) zygosity_source(x) else source
    dplyr::filter(raw, .data$zygosity_source %in% active_source)
})

#' @export
#' @rdname SNPData-class
get_donor_snp_info <- function(x, source = NULL) donor_snp_info(x, source = source)

#' @exportMethod chr_style
#' @rdname SNPData-class
setGeneric("chr_style", function(x) standardGeneric("chr_style"))
#' @exportMethod chr_style
#' @rdname SNPData-class
setMethod("chr_style", signature(x = "SNPData"), function(x) {
    # Handle backwards compatibility with older SNPData objects
    if (!methods::.hasSlot(x, "chr_style")) {
        return("unknown")
    }
    x@chr_style
})

#' @exportMethod zygosity_source
#' @rdname SNPData-class
setGeneric("zygosity_source", function(x) standardGeneric("zygosity_source"))
#' @exportMethod zygosity_source
#' @rdname SNPData-class
setMethod("zygosity_source", signature(x = "SNPData"), function(x) {
    # Handle backwards compatibility with older SNPData objects
    if (!methods::.hasSlot(x, "zygosity_source")) {
        return(NA_character_)
    }
    x@zygosity_source
})

#' @exportMethod zygosity_source<-
#' @rdname SNPData-class
setGeneric("zygosity_source<-", function(x, value) standardGeneric("zygosity_source<-"))
#' @exportMethod zygosity_source<-
#' @rdname SNPData-class
setReplaceMethod("zygosity_source", signature(x = "SNPData", value = "character"), function(x, value) {
    if (length(value) != 1 || is.na(value)) {
        stop("zygosity_source<- requires a single, non-NA character value")
    }
    available <- unique(stats::na.omit(donor_snp_info(x, source = "all")$zygosity_source))
    if (!value %in% available) {
        stop(sprintf(
            "zygosity_source<-: '%s' not found in donor_snp_info$zygosity_source. Available source(s): %s",
            value,
            if (length(available) > 0) paste(available, collapse = ", ") else "(none)"
        ))
    }
    x@zygosity_source <- value
    x
})

# Propagates the active zygosity_source slot from an existing object to one
# just rebuilt via new("SNPData", ...) -- the constructor's own
# auto-derivation (see .derive_zygosity_source()) collapses to NA once more
# than one source is present, so anything that reconstructs from an existing
# object (rather than building fresh) must carry an already-active source
# over explicitly rather than relying on re-derivation. If the old object had
# no active source yet (e.g. a direct add_donor_snp_metadata() call just
# introduced the first one, as infer_zygosity() does), it derives one from
# the *new* object's donor_snp_info instead of leaving it stuck at NA.
.propagate_zygosity_source <- function(new_x, old_x) {
    if (!methods::.hasSlot(old_x, "zygosity_source")) {
        return(new_x)
    }
    new_x@zygosity_source <- if (!is.na(old_x@zygosity_source)) {
        old_x@zygosity_source
    } else {
        .derive_zygosity_source(donor_snp_info(new_x, source = "all"))
    }
    new_x
}

#' @exportMethod updateObject
#' @rdname SNPData-class
setMethod("updateObject", signature(object = "SNPData"), function(object, ..., verbose = FALSE) {
    if (!methods::.hasSlot(object, "chr_style")) {
        snp_info <- object@snp_info
        if ("chrom" %in% colnames(snp_info)) {
            chr_style <- detect_chr_style(snp_info$chrom)
            if (!"chrom_canonical" %in% colnames(snp_info)) {
                snp_info$chrom_canonical <- normalise_chr_names(snp_info$chrom, from_style = chr_style)
            }
        } else {
            chr_style <- "unknown"
        }

        if (verbose) {
            log_info("updateObject(SNPData): detected chromosome style '{chr_style}'")
        }

        object@snp_info <- snp_info
        object@chr_style <- chr_style
    }

    if (!methods::.hasSlot(object, "donor_snp_info")) {
        snp_info <- object@snp_info
        packed_cols <- c(
            "xci_informative",
            "xci_informative_donor",
            "xci_allele_on_x1_by_donor",
            "xci_escape_fraction_by_donor"
        )

        if (all(packed_cols %in% colnames(snp_info))) {
            # Unpack the three CSV-packed *_by_donor columns into the long
            # donor_snp_info shape. Legacy objects never called zygosity
            # through this path, so those columns come back NA.
            donor_snp_info <- snp_info %>%
                dplyr::filter(xci_informative) %>%
                dplyr::select(
                    snp_id,
                    donor = xci_informative_donor,
                    allele_on_x1 = xci_allele_on_x1_by_donor,
                    xci_escape_fraction = xci_escape_fraction_by_donor
                ) %>%
                tidyr::separate_longer_delim(
                    c(donor, allele_on_x1, xci_escape_fraction),
                    delim = ","
                ) %>%
                dplyr::mutate(
                    xci_escape_fraction = as.numeric(xci_escape_fraction),
                    xci_informative = TRUE,
                    zygosity = NA_character_,
                    zygosity_source = NA_character_,
                    zygosity_p_val = NA_real_,
                    zygosity_adj_p_val = NA_real_,
                    zygosity_gt_prob = NA_real_
                )
            snp_info <- snp_info %>% dplyr::select(-dplyr::all_of(packed_cols))
        } else {
            donor_snp_info <- .empty_donor_snp_info()
        }

        if (verbose) {
            log_info(
                "updateObject(SNPData): migrated {nrow(donor_snp_info)} donor_snp_info row(s) from packed snp_info columns"
            )
        }

        object@snp_info <- snp_info
        object@donor_info <- .default_donor_info(object@barcode_info)
        object@donor_snp_info <- donor_snp_info
    }

    if (!methods::.hasSlot(object, "zygosity_source")) {
        zygosity_source <- .derive_zygosity_source(object@donor_snp_info)
        if (verbose) {
            log_info("updateObject(SNPData): derived zygosity_source '{zygosity_source}'")
        }
        object@zygosity_source <- zygosity_source
    }

    methods::validObject(object)
    object
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
    # Once donor_info carries data, the 'donor' column can only change via
    # rename_donor() -- a wholesale relabel, not an arbitrary per-cell
    # reassignment that would silently desync donor_info/donor_snp_info.
    if (methods::.hasSlot(x, "donor_info") && nrow(x@donor_info) > 0) {
        old_donor <- if ("donor" %in% colnames(x@barcode_info)) x@barcode_info$donor else NULL
        new_donor <- if ("donor" %in% colnames(value)) value$donor else NULL
        if (!identical(old_donor, new_donor)) {
            stop(
                "barcode_info<- cannot change the 'donor' column once donor_info carries data. ",
                "Use rename_donor() to relabel donors."
            )
        }
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

#' Relabel donors in a SNPData object
#'
#' Wholesale-relabels one or more donors, updating \code{barcode_info$donor},
#' \code{donor_info$donor}, and \code{donor_snp_info$donor} together so the three stay
#' consistent. Every cell currently labelled with an old donor moves to its new label;
#' individual cells cannot be reassigned between donors this way (\code{barcode_info<-}
#' rejects any 'donor' column change once \code{donor_info} carries data, precisely to
#' force donor-label changes through this single, atomic path).
#'
#' @param x A SNPData object
#' @param donor_map A named character vector, \code{c(new_label = old_label, ...)},
#'   following the same \code{new = old} convention as \code{dplyr::rename()}. Every value
#'   must be an existing donor (see \code{donor_info(x)}); donors not mentioned keep
#'   their current label. The resulting donor labels (renamed and unrenamed together) must
#'   stay unique -- \code{rename_donor()} relabels, it does not merge donors.
#'
#' @return A SNPData object with the relabelled donors
#' @export
#'
#' @examples
#' \dontrun{
#' # Vireo names donors arbitrarily (donor0, donor1, ...); give them real identities
#' snp_data <- rename_donor(snp_data, c(PatientA = "donor0", PatientB = "donor1"))
#' }
rename_donor <- function(x, donor_map) {
    if (!methods::is(x, "SNPData")) {
        stop("Input must be a SNPData object")
    }

    donor_info <- donor_info(x)
    unknown_old <- setdiff(donor_map, donor_info$donor)
    if (length(unknown_old) > 0) {
        stop(paste0(
            "donor_map references donor(s) not present in this object: ",
            paste(unknown_old, collapse = ", ")
        ))
    }

    new_labels <- .apply_donor_map(donor_info$donor, donor_map)
    if (any(duplicated(new_labels))) {
        stop(paste0(
            "donor_map produces duplicate donor label(s): ",
            paste(unique(new_labels[duplicated(new_labels)]), collapse = ", ")
        ))
    }

    barcode_info <- barcode_info(x)
    if ("donor" %in% colnames(barcode_info)) {
        barcode_info$donor <- .apply_donor_map(barcode_info$donor, donor_map)
    }
    donor_info$donor <- new_labels
    donor_snp_info <- donor_snp_info(x, source = "all")
    donor_snp_info$donor <- .apply_donor_map(donor_snp_info$donor, donor_map)

    result <- new(
        "SNPData",
        ref_count = x@ref_count,
        alt_count = x@alt_count,
        oth_count = x@oth_count,
        snp_info = x@snp_info,
        barcode_info = barcode_info,
        donor_info = donor_info,
        donor_snp_info = donor_snp_info
    )
    .propagate_zygosity_source(result, x)
}

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
        # Handle backwards compatibility with older SNPData objects
        if (methods::.hasSlot(object, "chr_style")) {
            cat("Chromosome style:", object@chr_style, "\n")
        }
        if (methods::.hasSlot(object, "zygosity_source")) {
            cat("Active zygosity source:", object@zygosity_source, "\n")
        }
        print(object@snp_info)
        cat("Barcode info (barcode_info()):", "\n")
        print(object@barcode_info)
        cat("Donor info (donor_info()):", "\n")
        print(donor_info(object))
        cat("Donor SNP info (donor_snp_info()):", "\n")
        print(donor_snp_info(object))
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

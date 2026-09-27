#' MoleculeCalls: S4 class for read-backed molecule allele calls
#'
#' Holds the per-molecule allele calls \code{\link{phase_from_molecules}}
#' extracts from a SNPData object's BAM files, together with what is needed to
#' count them by gene and to trace them back to their source. A SNPData object
#' carries one in its \code{molecules} slot, empty until
#' \code{phase_from_molecules} fills it.
#'
#' Calls are keyed on (\code{library_id}, \code{barcode}) rather than on donor,
#' so relabelling donors with \code{\link{rename_donor}} leaves them valid;
#' functions that count them look each cell's donor up from
#' \code{barcode_info}. A SNPData object passes its MoleculeCalls through
#' subsetting unchanged, since those functions only count calls for cells and
#' SNPs the object still holds.
#'
#' @param calls A data.frame, optional (default: an empty table). One row per
#'   (molecule, SNP), with columns \code{library_id} (character, \code{NA} for
#'   an unlabelled library), \code{barcode}, \code{umi}, \code{snp_id},
#'   \code{allele} (\code{"REF"}/\code{"ALT"}), and \code{transcript_strand}
#'   (\code{"+"}/\code{"-"}/\code{NA}).
#' @param snp_gene_map A data.frame, optional (default: an empty table). One
#'   row per (\code{snp_id}, candidate \code{gene_name}), with columns
#'   \code{snp_id}, \code{gene_name}, \code{gene_strand}, \code{ambiguous}, as
#'   returned by \code{\link{assign_snp_genes}}.
#' @param bam_calibration A data.frame, optional (default: an empty table). One
#'   row per BAM file scanned, with columns \code{bam_file}, \code{orientation},
#'   \code{n_ts_reads}, \code{concordance}, \code{n_scanned}.
#' @param bam_files A named list of character vectors, optional (default: an
#'   empty list). \code{library_id = path(s)}: the BAM files the calls were
#'   extracted from.
#' @param params A named list, optional (default: an empty list). The
#'   extraction settings used, recorded for provenance.
#' @param object A MoleculeCalls object, required. Passed to the show method.
#' @param x A MoleculeCalls object, required.
#'
#' @slot calls A tibble of per-molecule allele calls (see the \code{calls} argument).
#' @slot snp_gene_map A tibble mapping SNPs to candidate genes (see the
#'   \code{snp_gene_map} argument). Distinct from \code{snp_info$gene_name},
#'   which comma-joins overlapping genes into one display label: this keeps each
#'   candidate as its own row with its strand, which is what attributing a
#'   molecule at a multi-gene SNP needs.
#' @slot bam_calibration A tibble of per-BAM strand-orientation diagnostics.
#' @slot bam_files A named list of the BAM paths the calls came from.
#' @slot params A named list of the extraction settings used.
#'
#' @section Accessors:
#' \describe{
#'   \item{\code{molecule_calls(x)}}{Get the per-molecule allele calls}
#'   \item{\code{snp_gene_map(x)}}{Get the SNP-to-gene map}
#'   \item{\code{bam_calibration(x)}}{Get the per-BAM strand calibration}
#' }
#'
#' @seealso \code{\link{molecules}} to get a SNPData object's MoleculeCalls.
#'
#' @examples
#' # An empty container, as held by a SNPData object before phasing
#' MoleculeCalls()
#'
#' @family molecule-level allele counting functions
#' @exportClass MoleculeCalls
#' @export
setClass(
    "MoleculeCalls",
    slots = c(
        calls = "tbl_df",
        snp_gene_map = "tbl_df",
        bam_calibration = "tbl_df",
        bam_files = "list",
        params = "list"
    )
)

.MOLECULE_CALL_COLUMNS <- c("library_id", "barcode", "umi", "snp_id", "allele", "transcript_strand")
.SNP_GENE_MAP_COLUMNS <- c("snp_id", "gene_name", "gene_strand", "ambiguous")

setValidity("MoleculeCalls", function(object) {
    problems <- c(
        .missing_columns_problem(object@calls, .MOLECULE_CALL_COLUMNS, "calls"),
        .missing_columns_problem(object@snp_gene_map, .SNP_GENE_MAP_COLUMNS, "snp_gene_map")
    )
    if (length(problems) == 0) TRUE else problems
})

.missing_columns_problem <- function(df, required, df_name) {
    missing_cols <- setdiff(required, colnames(df))
    if (length(missing_cols) == 0) {
        return(character(0))
    }
    paste0(df_name, " is missing required column(s): ", paste(missing_cols, collapse = ", "))
}

.empty_molecule_calls <- function() {
    tibble::tibble(
        library_id = character(0),
        barcode = character(0),
        umi = character(0),
        snp_id = character(0),
        allele = character(0),
        transcript_strand = character(0)
    )
}

.empty_snp_gene_map <- function() {
    tibble::tibble(
        snp_id = character(0),
        gene_name = character(0),
        gene_strand = character(0),
        ambiguous = logical(0)
    )
}

.empty_bam_calibration <- function() {
    tibble::tibble(
        bam_file = character(0),
        orientation = character(0),
        n_ts_reads = integer(0),
        concordance = double(0),
        n_scanned = integer(0)
    )
}

#' @rdname MoleculeCalls-class
#' @export
MoleculeCalls <- function(
    calls = NULL,
    snp_gene_map = NULL,
    bam_calibration = NULL,
    bam_files = list(),
    params = list()
) {
    new(
        "MoleculeCalls",
        calls = .as_tibble_or(calls, .empty_molecule_calls()),
        snp_gene_map = .as_tibble_or(snp_gene_map, .empty_snp_gene_map()),
        bam_calibration = .as_tibble_or(bam_calibration, .empty_bam_calibration()),
        bam_files = as.list(bam_files),
        params = params
    )
}

# A table argument as a tibble, or its empty template when not supplied.
.as_tibble_or <- function(value, empty) {
    if (is.null(value)) {
        return(empty)
    }
    tibble::as_tibble(value)
}

# TRUE once phase_from_molecules() has stored calls; an empty container is the
# default every SNPData object starts with.
.has_molecule_calls <- function(molecules) {
    nrow(molecules@calls) > 0
}

#' @exportMethod molecule_calls
#' @rdname MoleculeCalls-class
setGeneric("molecule_calls", function(x) standardGeneric("molecule_calls"))
#' @exportMethod molecule_calls
#' @rdname MoleculeCalls-class
setMethod("molecule_calls", signature(x = "MoleculeCalls"), function(x) x@calls)

#' @exportMethod snp_gene_map
#' @rdname MoleculeCalls-class
setGeneric("snp_gene_map", function(x) standardGeneric("snp_gene_map"))
#' @exportMethod snp_gene_map
#' @rdname MoleculeCalls-class
setMethod("snp_gene_map", signature(x = "MoleculeCalls"), function(x) x@snp_gene_map)

#' @exportMethod bam_calibration
#' @rdname MoleculeCalls-class
setGeneric("bam_calibration", function(x) standardGeneric("bam_calibration"))
#' @exportMethod bam_calibration
#' @rdname MoleculeCalls-class
setMethod("bam_calibration", signature(x = "MoleculeCalls"), function(x) x@bam_calibration)

#' @exportMethod show
#' @rdname MoleculeCalls-class
setMethod("show", signature(object = "MoleculeCalls"), function(object) {
    cat("Object of class 'MoleculeCalls'", "\n")
    if (!.has_molecule_calls(object)) {
        cat("Empty: no molecule calls extracted (see phase_from_molecules())", "\n")
        return(invisible(NULL))
    }
    calls <- object@calls
    cat(
        nrow(calls),
        " calls from ",
        nrow(dplyr::distinct(calls, library_id, barcode, umi)),
        " molecules at ",
        dplyr::n_distinct(calls$snp_id),
        " SNPs",
        "\n",
        sep = ""
    )
    cat("SNP-to-gene map:", nrow(object@snp_gene_map), "rows", "\n")
    for (library_id in names(object@bam_files)) {
        cat("BAM files [", library_id, "]: ", paste(object@bam_files[[library_id]], collapse = ", "), "\n", sep = "")
    }
    invisible(NULL)
})

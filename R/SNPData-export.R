#' Export SNPData object to cellSNP-compatible files
#'
#' Writes the count matrices and cell annotations of a SNPData object to an
#' output folder laid out like cellSNP-lite output, so that the data can be
#' handed to external tools that read that format (and re-read by
#' \code{\link{import_cellsnp}}).
#'
#' @param snpdata A SNPData object, required. Must hold cells from a single
#'   sequencing library; see Details.
#' @param out_dir Character scalar, required. Output directory to write files to.
#'
#' @details
#' The cellSNP format describes one library: it keys cells by barcode alone,
#' and carries no slot for anything the analysis has since derived. Exporting
#' is therefore lossy, and deliberately so — this writes an interchange format,
#' not a saved object. Use \code{saveRDS()} to round-trip a SNPData object
#' without loss.
#'
#' Written out:
#' \itemize{
#'   \item \code{cellSNP.tag.AD.mtx}, \code{cellSNP.tag.DP.mtx},
#'     \code{cellSNP.tag.OTH.mtx} — ALT, total (REF + ALT), and OTH counts.
#'   \item \code{cellSNP.base.vcf.gz} — \code{chrom}, \code{pos}, \code{snp_id},
#'     \code{ref}, \code{alt} from \code{snp_info}.
#'   \item \code{cellSNP.samples.tsv} — one barcode per line, in matrix column
#'     order. Kept single-column, as external readers expect.
#'   \item \code{donor_ids.tsv} — one row per cell, in matrix column order,
#'     carrying \code{cell} and \code{donor_id} as Vireo writes them, plus this
#'     object's \code{cell_id} and \code{library_id} so the exported directory
#'     records which library it came from. \code{import_cellsnp} ignores those
#'     two extra columns and takes \code{library_id} from its own argument.
#'   \item \code{filtered_contig_annotations.csv} — \code{barcode} and
#'     \code{raw_clonotype_id}, written only when clonotypes are present.
#' }
#'
#' Dropped, with no cellSNP-format equivalent: \code{donor_snp_info} (every
#' zygosity call and every \code{\link{assign_xci}} result), the active
#' \code{\link{zygosity_source}}, \code{library_info} (including recorded BAM
#' paths), \code{donor_info}, and any columns of \code{snp_info} or
#' \code{barcode_info} beyond those listed above — notably \code{gene_name},
#' which \code{import_cellsnp} regenerates from its \code{gene_annotation}
#' argument.
#'
#' A multi-library object cannot be represented: a 10x barcode is unique only
#' within its library, so writing two libraries into one barcode-keyed
#' directory would fuse or mis-assign cells that happen to share a barcode.
#' Such an object is rejected rather than exported; split it first, e.g.
#' \code{filter_barcodes(x, library_id == "run1")}.
#'
#' @family import and export functions
#' @export
#'
#' @examples
#' \dontrun{
#' snp_data <- get_example_snpdata()
#' export_cellsnp(snp_data, "exported_cellsnp")
#'
#' # A merged object must be split into its libraries first
#' export_cellsnp(filter_barcodes(merged, library_id == "run1"), "run1_cellsnp")
#' }
export_cellsnp <- function(snpdata, out_dir) {
    barcode_info <- barcode_info(snpdata)
    .check_single_library(barcode_info)
    .check_unique_barcodes(barcode_info)

    if (!dir.exists(out_dir)) {
        dir.create(out_dir, recursive = TRUE)
    }
    logger::log_info("Exporting SNPData to {out_dir}")

    # Write ALT, REF, and OTH matrices as Matrix Market files
    ad_file <- file.path(out_dir, "cellSNP.tag.AD.mtx")
    dp_file <- file.path(out_dir, "cellSNP.tag.DP.mtx")
    oth_file <- file.path(out_dir, "cellSNP.tag.OTH.mtx")
    Matrix::writeMM(alt_count(snpdata), ad_file)
    logger::log_info("ALT count matrix written to: {ad_file}")
    Matrix::writeMM(ref_count(snpdata) + alt_count(snpdata), dp_file)
    logger::log_info("DP (total) count matrix written to: {dp_file}")
    Matrix::writeMM(oth_count(snpdata), oth_file)
    logger::log_info("OTH count matrix written to: {oth_file}")

    # Write SNP info as a VCF-like file (minimal, for import_cellsnp)
    snp_info <- snpdata@snp_info
    vcf_file <- file.path(out_dir, "cellSNP.base.vcf")
    vcf_df <- snp_info %>%
        dplyr::transmute(
            chrom = chrom,
            pos = pos,
            id = snp_id,
            ref = ref,
            alt = alt,
            qual = ".",
            filter = ".",
            info = "."
        )
    readr::write_tsv(vcf_df, vcf_file, col_names = FALSE)
    logger::log_info("SNP info written to: {vcf_file} (gzipping...)")
    R.utils::gzip(vcf_file, overwrite = TRUE)
    logger::log_info("Gzipped VCF file written to: {vcf_file}.gz")

    # One row per cell, in matrix column order: the count matrices are the
    # authority on which cell is which, and de-duplicating here would let the
    # annotation tables fall out of step with them.
    donor_file <- file.path(out_dir, "donor_ids.tsv")
    donor_df <- barcode_info %>%
        dplyr::transmute(
            cell = barcode,
            donor_id = donor,
            cell_id = cell_id,
            library_id = library_id
        )
    readr::write_tsv(donor_df, donor_file)
    logger::log_info("Donor info written to: {donor_file}")

    # Write VDJ info as filtered_contig_annotations.csv (if clonotype data available)
    if ("clonotype" %in% colnames(barcode_info) && !all(is.na(barcode_info$clonotype))) {
        vdj_file <- file.path(out_dir, "filtered_contig_annotations.csv")
        vdj_df <- barcode_info %>%
            dplyr::transmute(barcode = barcode, raw_clonotype_id = clonotype)
        readr::write_csv(vdj_df, vdj_file)
        logger::log_info("VDJ info written to: {vdj_file}")
    } else {
        logger::log_info("Skipping VDJ export (no clonotype information available)")
    }

    # Write barcodes into cellSNP.samples.tsv. Single column, unlike
    # donor_ids.tsv, because external cellSNP readers expect a bare list.
    samples_file <- file.path(out_dir, "cellSNP.samples.tsv")
    readr::write_tsv(
        dplyr::select(barcode_info, barcode),
        samples_file,
        col_names = FALSE
    )

    .log_dropped_slots(snpdata)
    logger::log_success("SNPData exported to {out_dir}")
}

# A cellSNP directory describes one library, and keys cells by barcode alone.
# Two libraries in one directory would silently fuse or mis-assign the cells
# that happen to share a barcode, so refuse instead of writing a wrong answer.
.check_single_library <- function(barcode_info) {
    if (!"library_id" %in% colnames(barcode_info)) {
        return(invisible(NULL))
    }
    libraries <- unique(stats::na.omit(barcode_info$library_id))
    if (length(libraries) <= 1) {
        return(invisible(NULL))
    }
    stop(
        "export_cellsnp() cannot write a multi-library object: the cellSNP format keys cells by ",
        "barcode, which is unique only within a library. Export each library separately, e.g. ",
        sprintf(
            "filter_barcodes(x, library_id == \"%s\"). Libraries present: %s.",
            libraries[1],
            paste(libraries, collapse = ", ")
        )
    )
}

# Even within one library the annotation files are barcode-keyed, so a repeated
# barcode cannot be written unambiguously.
.check_unique_barcodes <- function(barcode_info) {
    duplicates <- unique(barcode_info$barcode[duplicated(barcode_info$barcode)])
    if (length(duplicates) == 0) {
        return(invisible(NULL))
    }
    stop(
        "export_cellsnp() cannot write repeated barcodes, which the cellSNP format cannot tell ",
        sprintf(
            "apart: %s%s.",
            paste(utils::head(duplicates, 5), collapse = ", "),
            if (length(duplicates) > 5) ", ..." else ""
        )
    )
}

# Exporting is lossy by design; say which analysis results are being left
# behind rather than letting them disappear quietly.
.log_dropped_slots <- function(snpdata) {
    dropped <- character(0)
    if (nrow(donor_snp_info(snpdata, source = "all")) > 0) {
        dropped <- c(dropped, "donor_snp_info (zygosity calls and XCI results)")
    }
    if (any(lengths(library_info(snpdata)$bam_files) > 0)) {
        dropped <- c(dropped, "library_info BAM paths")
    }
    if (nrow(snp_gene_map(snpdata)) > 0) {
        dropped <- c(dropped, "the SNP-to-gene map")
    }
    if (length(dropped) > 0) {
        logger::log_warn(
            "The cellSNP format cannot carry {paste(dropped, collapse = ' and ')}; ",
            "use saveRDS() to keep them."
        )
    }
}

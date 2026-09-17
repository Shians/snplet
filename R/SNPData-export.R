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

#' Export a SNPData object's cell-level counts as a SingleCellExperiment
#'
#' Wraps the REF/ALT (and OTH, if present) count matrices as assays of a
#' \code{SingleCellExperiment}, with \code{barcode_info} as \code{colData} and
#' \code{snp_info} as \code{rowData}, so cell-level allelic calls (e.g.
#' \code{\link{assign_xci}}'s \code{active_x}) can be carried into
#' Seurat (\code{Seurat::as.Seurat()}) or other SingleCellExperiment-based
#' tools alongside a cell's existing cluster/UMAP annotation.
#'
#' @details
#' This is the cell-level counterpart to \code{\link{as_escape_experiment}},
#' which exports donor-level escape counts instead; the two are not
#' interchangeable, since a \code{SingleCellExperiment} is cell-keyed and an
#' XCI escape comparison is donor-keyed (see \code{\link{as_escape_experiment}}
#' for why cell-level counts cannot stand in for that).
#'
#' Assays: \code{"ref"} and \code{"alt"} always; \code{"oth"} only when
#' \code{oth_count(snpdata)} has any non-zero entry, since most objects carry
#' an all-zero \code{oth_count} and an all-zero assay would be dead weight IN
#' every downstream operation. \code{rowData} and \code{colData} carry every
#' column of \code{snp_info}/\code{barcode_info} respectively, other than the
#' identifiers already used as dimnames (\code{snp_id}, \code{cell_id}).
#'
#' Dropped, with no \code{SingleCellExperiment} slot to hold them:
#' \code{donor_info}, \code{donor_snp_info} (zygosity calls and XCI gene-level
#' diagnostics), \code{library_info}, and \code{snp_gene_map}. Use
#' \code{saveRDS()} to keep the full SNPData object, or re-derive per-cell
#' summaries of these (e.g. \code{donor_info}'s \code{xci_skew} is one value
#' per donor, not per cell, so it has no natural column here; join it onto
#' \code{colData} afterwards by \code{donor} if needed).
#'
#' @param snpdata A SNPData object, required.
#'
#' @return A \code{SingleCellExperiment} with SNPs as rows and cells as
#'   columns, dimnamed by \code{snp_id} and \code{cell_id}.
#'
#' @family import and export functions
#' @export
#'
#' @examples
#' \dontrun{
#' snp_data <- assign_xci(get_example_snpdata())
#' sce <- as_singlecellexperiment(snp_data)
#'
#' # active_x rides along in colData, ready to plot against an existing
#' # Seurat object's clustering once merged
#' SummarizedExperiment::colData(sce)$active_x
#' }
as_singlecellexperiment <- function(snpdata) {
    if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
        stop(
            "as_singlecellexperiment() requires the SingleCellExperiment package. ",
            "Install it with BiocManager::install(\"SingleCellExperiment\")."
        )
    }

    snp_info <- snp_info(snpdata)
    barcode_info <- barcode_info(snpdata)

    assays <- list(ref = ref_count(snpdata), alt = alt_count(snpdata))
    oth <- oth_count(snpdata)
    if (any(oth != 0)) {
        assays$oth <- oth
    }

    row_data <- DataFrame(
        dplyr::select(snp_info, -snp_id),
        row.names = snp_info$snp_id
    )
    col_data <- DataFrame(
        dplyr::select(barcode_info, -cell_id),
        row.names = barcode_info$cell_id
    )

    SingleCellExperiment::SingleCellExperiment(
        assays = assays,
        rowData = row_data,
        colData = col_data
    )
}

#' Export donor-level XCI escape counts as a SummarizedExperiment
#'
#' Pivots \code{\link{test_escape}}'s per-(donor, gene) escape counts into a
#' gene x donor \code{SummarizedExperiment} with \code{"active"} and
#' \code{"inactive"} assays, so escape (or skew) can be compared between
#' experimental groups with external tools built for paired-count data, most
#' directly edgeR's differential methylation workflow
#' (\code{edgeR::modelMatrixMeth()}, treating \code{"active"}/\code{"inactive"}
#' the way that workflow treats methylated/unmethylated read counts per CpG).
#'
#' @details
#' This is the donor-level counterpart to
#' \code{\link{as_singlecellexperiment}}. The two cannot substitute for each
#' other: an escape comparison needs donors as replicates (independent units a
#' group label can attach to and a test can be powered over), and a donor's
#' cells cannot serve as those replicates themselves, since they are not
#' independent — they share one donor's XCI skew and one donor's genotype, the
#' entire premise \code{\link{assign_xci}} fits per donor rather than
#' per cell.
#'
#' \code{\link{haplotype_expression}} (or
#' \code{\link{haplotype_expression_by_molecule}}, used automatically once
#' \code{\link{add_molecule_phase}} has run, matching \code{\link{test_escape}}'s
#' own source selection) excludes a (donor, gene) pair with no qualifying SNP
#' entirely, rather than reporting zero coverage for it — so the gene x donor
#' grid is not dense before this function pads it. A missing pair becomes
#' \code{NA} in both assays here, not a dropped row: a
#' \code{SummarizedExperiment}'s assays must share one set of dimnames, so a
#' gene present for some donors and absent for others cannot be represented by
#' varying which rows exist per column, only by \code{NA} in the columns where
#' it is missing.
#'
#' \code{colData} is \code{donor_info(snpdata)}, so an experimental-group
#' label attached beforehand with
#' \code{add_donor_metadata(x, data.frame(donor = ..., group = ...))} carries
#' straight into the result, alongside \code{xci_skew} and the other stored
#' per-donor diagnostics. \code{rowData} has one row per \code{gene_name}.
#'
#' edgeR has no native concept of two paired assays: its differential
#' methylation path expects one matrix with \code{<sample>-Me}/\code{<sample>-Un}
#' column pairs (see \code{edgeR::modelMatrixMeth()}), so combine this
#' object's two assays into that shape before constructing a
#' \code{edgeR::DGEList()} (see Examples); this function stops at the
#' \code{SummarizedExperiment}, which is the natural, tool-agnostic
#' Bioconductor container for the gene x donor active/inactive counts
#' themselves, and does not depend on edgeR.
#'
#' NA-containing rows/columns are not filtered here; edgeR's own count-based
#' functions generally do not tolerate NA counts, so filter to genes with
#' complete coverage across the donors being compared before fitting.
#'
#' @param snpdata A SNPData object, required, that had XCI diagnostics stored
#'   by \code{\link{assign_xci}} or \code{\link{assign_xci_by_clonotype}}.
#'
#' @return A \code{SummarizedExperiment} with genes as rows and donors as
#'   columns, assays \code{"active"} and \code{"inactive"} (integer, \code{NA}
#'   where that gene had no qualifying SNP in that donor).
#'
#' @family X-chromosome inactivation functions
#' @export
#'
#' @examples
#' \dontrun{
#' snp_data <- assign_xci(snp_data)
#' snp_data <- add_donor_metadata(
#'   snp_data,
#'   data.frame(donor = c("donor0", "donor1"), group = c("control", "treated"))
#' )
#' se <- as_escape_experiment(snp_data)
#'
#' # Reshape into edgeR's methylation-style Me/Un column-pair matrix and fit
#' active <- SummarizedExperiment::assay(se, "active")
#' inactive <- SummarizedExperiment::assay(se, "inactive")
#' donors <- colnames(se)
#' counts <- do.call(cbind, lapply(donors, function(d) {
#'   m <- cbind(active[, d], inactive[, d])
#'   colnames(m) <- paste0(d, c("-Me", "-Un"))
#'   m
#' }))
#' keep <- stats::complete.cases(counts)
#' y <- edgeR::DGEList(counts[keep, ])
#' design_samples <- model.matrix(~0 + group, data = as.data.frame(SummarizedExperiment::colData(se)))
#' design <- edgeR::modelMatrixMeth(design_samples)
#' y <- edgeR::estimateDisp(y, design, trend = "none")
#' fit <- edgeR::glmFit(y, design)
#' }
as_escape_experiment <- function(snpdata) {
    if (!.has_xci_diagnostics(snpdata)) {
        stop("No stored XCI diagnostics found. Run assign_xci(snpdata) first.")
    }

    counts <- .escape_counts(snpdata)
    if (!"donor" %in% colnames(counts)) {
        stop("as_escape_experiment() needs donor assignments to build a gene x donor matrix.")
    }

    donors <- sort(unique(counts$donor))
    genes <- sort(unique(counts$gene_name))

    # haplotype_expression()/haplotype_expression_by_molecule() drop a (donor,
    # gene) pair with no qualifying SNP rather than reporting zero coverage
    # for it, so the grid is padded here rather than already dense: a missing
    # pair becomes NA, not a false zero, since zero coverage and "never
    # measured" are not the same claim.
    dense <- tidyr::expand_grid(donor = donors, gene_name = genes) %>%
        dplyr::left_join(counts, by = c("donor", "gene_name"))

    to_matrix <- function(value_col) {
        m <- dense %>%
            dplyr::select(gene_name, donor, dplyr::all_of(value_col)) %>%
            tidyr::pivot_wider(names_from = donor, values_from = dplyr::all_of(value_col)) %>%
            tibble::column_to_rownames("gene_name") %>%
            as.matrix()
        m[genes, donors, drop = FALSE]
    }

    assays <- list(active = to_matrix("active_count"), inactive = to_matrix("inactive_count"))

    donor_info <- donor_info(snpdata)
    col_data <- DataFrame(
        dplyr::select(donor_info, -donor)[match(donors, donor_info$donor), , drop = FALSE],
        row.names = donors
    )
    row_data <- DataFrame(gene_name = genes, row.names = genes)

    SummarizedExperiment(
        assays = assays,
        rowData = row_data,
        colData = col_data
    )
}

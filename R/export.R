#' Export SNPData object to cellSNP-compatible files
#'
#' Writes the components of a SNPData object to an output folder in a format
#' compatible with import_cellsnp.
#'
#' @param snpdata A SNPData object
#' @param out_dir Output directory to write files
#' @export
#'
#' @examples
#' \dontrun{
#' export_cellsnp(snp_data, "exported_cellsnp")
#' }
export_cellsnp <- function(snpdata, out_dir) {
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
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

    # Write donor info as donors.tsv (if available)
    barcode_info <- get_barcode_info(snpdata)
    donor_file <- file.path(out_dir, "donor_ids.tsv")
    donor_df <- barcode_info %>%
        dplyr::select(barcode, donor) %>%
        dplyr::rename(donor_id = donor) %>%
        dplyr::distinct()
    readr::write_tsv(donor_df, donor_file)
    logger::log_info("Donor info written to: {donor_file}")

    # Write VDJ info as filtered_contig_annotations.csv
    vdj_file <- file.path(out_dir, "filtered_contig_annotations.csv")
    vdj_df <- barcode_info %>%
        dplyr::select(barcode, clonotype) %>%
        dplyr::rename(raw_clonotype_id = clonotype) %>%
        dplyr::distinct()
    readr::write_csv(vdj_df, vdj_file)
    logger::log_info("VDJ info written to: {vdj_file}")

    # Write barcodes into cellSNP.samples.tsv
    samples_file <- file.path(out_dir, "cellSNP.samples.tsv")
    cells <- get_barcode_info(snpdata)
    readr::write_tsv(
        dplyr::select(cells, barcode), samples_file, col_names = FALSE
    )

    logger::log_success("SNPData exported to {out_dir}")
}

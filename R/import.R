#' Import cellSNP data and create a SNPData object
#'
#' This function imports data from cellSNP-lite output, VDJ annotations from cellranger,
#' and donor information from Vireo to create a SNPData object.
#'
#' @param cellsnp_dir Directory containing cellSNP-lite output files
#' @param vdj_file Path to filtered_contig_annotations.csv from cellranger VDJ
#' @param gene_annotation Data frame with gene annotations for SNPs (must have same number of rows as SNP matrices)
#' @param vireo_file Path to donors.tsv file from Vireo (optional, default: NA)
#' @param barcode_column Name of the column in vdj_file that contains cell barcodes
#' @param clonotype_column Name of the column in vdj_file that contains clonotype information
#'
#' @return A SNPData object
#' @export
#'
#' @examples
#' \dontrun{
#' snp_data <- import_cellsnp(
#'   cellsnp_dir = "path/to/cellsnp_output",
#'   vdj_file = "path/to/filtered_contig_annotations.csv",
#'   gene_annotation = gene_anno_df
#' )
#' }
import_cellsnp <- function(
    cellsnp_dir,
    vdj_file,
    gene_annotation,
    vireo_file = NA,
    barcode_column = "barcode",
    clonotype_column = "raw_clonotype_id"
) {
    # Validate gene_annotation columns
    required_gene_cols <- c("chrom", "gene_name")
    missing_cols <- setdiff(required_gene_cols, colnames(gene_annotation))
    if (length(missing_cols) > 0) {
        stop(
            sprintf(
                "gene_annotation is missing required columns: %s",
                paste(missing_cols, collapse = ", ")
            )
        )
    }

    # Check if required files exist
    dp_file <- fs::path(cellsnp_dir, "cellSNP.tag.DP.mtx")
    ad_file <- fs::path(cellsnp_dir, "cellSNP.tag.AD.mtx")
    oth_file <- fs::path(cellsnp_dir, "cellSNP.tag.OTH.mtx")
    base_file <- fs::path(cellsnp_dir, "cellSNP.base.vcf.gz")
    samples_file <- fs::path(cellsnp_dir, "cellSNP.samples.tsv")

    for (file in c(dp_file, ad_file, oth_file, base_file, vdj_file)) {
        check_file(file)
    }
    # Only check vireo_file if not NA
    if (!is.na(vireo_file)) {
        check_file(vireo_file)
    }

    # Read cellSNP matrices
    coverage <- Matrix::readMM(dp_file)
    alt_count <- Matrix::readMM(ad_file)
    oth_count <- Matrix::readMM(oth_file)
    ref_count <- coverage - alt_count # Only subtract alt_count

    # Read SNP information from VCF file
    snp_vcf_data <- read_vcf_base(base_file)

    # Read cell barcodes from cells file
    cells <- readr::read_tsv(
        samples_file,
        col_names = "barcode",
        col_types = readr::cols(barcode = readr::col_character())
    )

    # Merge SNP info with gene annotation
    # Create a GRanges object for SNPs
    snps_gr <- plyranges::as_granges(
        snp_vcf_data,
        seqnames = chrom,
        start = pos,
        end = pos,
    )

    # Create a GRanges object for gene annotations
    gene_anno_gr <- plyranges::as_granges(
        gene_annotation,
        seqnames = chrom
    )

    # Merge SNPs with gene annotations
    snp_info <- plyranges::join_overlap_left(snps_gr, gene_anno_gr) %>%
        tibble::as_tibble() %>%
        dplyr::rename(chrom = seqnames, pos = start) %>%
        dplyr::select(snp_id, chrom, pos, ref, alt, gene_name) %>%
        dplyr::summarise(
            chrom = dplyr::first(chrom),
            pos = dplyr::first(pos),
            ref = dplyr::first(ref),
            alt = dplyr::first(alt),
            gene_name = paste(unique(gene_name), collapse = ", "),
            .by = snp_id
        ) %>%
        dplyr::mutate(
            gene_name = dplyr::if_else(gene_name == "NA", NA_character_, gene_name)
        )

    # Read donor information if provided, else create dummy donor info
    if (!is.na(vireo_file)) {
        donor_info <- readr::read_tsv(
            vireo_file,
            col_types = readr::cols(.default = readr::col_character())
        )
    } else {
        # get barcodes from cell_snp output
        donor_info <- cells %>%
            dplyr::mutate(donor = "donor0")
    }

    # Read VDJ clonotype information
    vdj_info <- readr::read_csv(
        vdj_file,
        col_types = readr::cols(.default = readr::col_character())
    ) %>%
        dplyr::mutate(
            barcode = stringr::str_remove(barcode, "-[0-9]+$") # Remove suffix if present
        )

    # Merge donor and clonotype information
    barcode_info <- merge_cell_annotations(
        donor_info,
        vdj_info,
        barcode_column,
        clonotype_column
    )

    # Create SNPData object
    logger::log_info("Creating SNPData object with {nrow(barcode_info)} barcodes and {nrow(snp_info)} SNPs")
    snp_data <- SNPData(
        alt_count = alt_count,
        ref_count = ref_count,
        oth_count = oth_count,
        snp_info = snp_info,
        barcode_info = barcode_info
    )

    return(snp_data)
}

#' Read the base VCF file from cellSNP output
#'
#' @param vcf_file Path to cellSNP.base.vcf.gz file
#'
#' @return Data frame with SNP information
#' @keywords internal
read_vcf_base <- function(vcf_file) {
    vcf_data <- readr::read_tsv(
        vcf_file,
        comment = "#",
        col_names = c(
            "chrom",
            "pos",
            "id",
            "ref",
            "alt",
            "qual",
            "filter",
            "info"
        ),
        col_types = readr::cols(
            chrom = readr::col_character(),
            pos = readr::col_integer(),
            id = readr::col_character(),
            ref = readr::col_character(),
            alt = readr::col_character(),
            qual = readr::col_character(),
            filter = readr::col_character(),
            info = readr::col_character()
        )
    )

    # Generate SNP IDs
    vcf_data$snp_id <- paste0("snp_", seq_len(nrow(vcf_data)))

    # Reorder columns to have snp_id first
    vcf_data <- vcf_data[, c("snp_id", names(vcf_data)[names(vcf_data) != "snp_id"])]

    return(vcf_data)
}

#' Merge donor and clonotype information
#'
#' @param donor_info Data frame with donor information from Vireo
#' @param vdj_info Data frame with VDJ information from cellranger
#' @param barcode_column Name of the column in vdj_info containing cell barcodes
#' @param clonotype_column Name of the column in vdj_info containing clonotype information
#'
#' @return Data frame with merged cell annotations
#' @keywords internal
merge_cell_annotations <- function(donor_info, vdj_info, barcode_column, clonotype_column) {
    # Standardize column names
    if ("cell" %in% colnames(donor_info)) {
        donor_info <- donor_info %>%
            dplyr::rename(
                barcode = cell,
                donor = donor_id
            )
    }

    # Ensure barcode_column exists in vdj_info
    if (!barcode_column %in% colnames(vdj_info)) {
        stop(paste("Column", barcode_column, "not found in VDJ annotation file"))
    }

    # Ensure clonotype_column exists in vdj_info
    if (!clonotype_column %in% colnames(vdj_info)) {
        stop(paste("Column", clonotype_column, "not found in VDJ annotation file"))
    }

    # Extract relevant columns from VDJ data
    vdj_subset <- vdj_info %>%
        dplyr::select(!!barcode_column, !!clonotype_column) %>%
        dplyr::rename(
            barcode = !!barcode_column,
            clonotype = !!clonotype_column
        ) %>%
        dplyr::distinct()

    # Merge donor info with VDJ info
    barcode_info <- donor_info %>%
        dplyr::left_join(vdj_subset, by = "barcode") %>%
        dplyr::mutate(
            cell_id = paste0("cell_", seq_len(dplyr::n()))
        )

    # Ensure required columns exist
    if (!"donor" %in% colnames(barcode_info)) {
        barcode_info$donor <- NA_character_
    }

    if (!"clonotype" %in% colnames(barcode_info)) {
        barcode_info$clonotype <- NA_character_
    }

    barcode_info <- barcode_info %>%
        dplyr::select(cell_id, barcode, donor, clonotype, everything())

    return(barcode_info)
}

#' Load example SNPData object for demonstration and testing
#'
#' This function loads a small example SNPData object using bundled example files
#' from the snplet package. It demonstrates the import workflow and is useful for
#' testing and vignettes.
#'
#' @return A SNPData object constructed from example data included with the package.
#' @examples
#' snp_data <- get_example_snpdata()
get_example_snpdata <- function() {
    import_cellsnp(
        cellsnp_dir = system.file("extdata/example_snpdata", package = "snplet"),
        vdj_file = system.file("extdata/example_snpdata/filtered_contig_annotations.csv", package = "snplet"),
        gene_annotation = readr::read_tsv(system.file("extdata/example_gene_anno.tsv", package = "snplet")),
        vireo_file = system.file("extdata/example_snpdata/donor_ids.tsv", package = "snplet")
    )
}

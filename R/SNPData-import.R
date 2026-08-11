#' Import cellSNP data and create a SNPData object
#'
#' This function imports data from cellSNP-lite output, with optional VDJ annotations from cellranger
#' and donor information from Vireo to create a SNPData object.
#'
#' @param cellsnp_dir Character scalar, required. Directory containing
#'   cellSNP-lite output files.
#' @param gene_annotation A data.frame, required, with columns \code{chrom},
#'   \code{start}, \code{end}, \code{gene_name}. Gene annotations.
#' @param library_id Character scalar, required. Name of the sequencing
#'   library this cellSNP-lite run came from, stored on every cell. A 10x
#'   barcode is unique only within its library, so \code{\link{merge_snpdata}}
#'   needs this label to tell a repeated cell from two different cells that
#'   happened to draw the same barcode.
#' @param vdj_file Character scalar, optional (default \code{NULL}). Path to
#'   \code{filtered_contig_annotations.csv} from cellranger VDJ.
#' @param vireo_folder Character scalar, optional (default \code{NULL}).
#'   Path to a Vireo output directory. Donor assignments are read from
#'   \code{donor_ids.tsv} inside it. If it also contains a genotype VCF
#'   (\code{GT_donors.vireo.vcf.gz}, one sample column per donor), that is
#'   used to populate per-(SNP, donor) zygosity calls at construction time;
#'   if the genotype VCF is absent, only donor assignments are read.
#' @param donor_map A named character vector, \code{c(new_label = old_label, ...)}
#'   (the same \code{new = old} convention as \code{dplyr::rename()}),
#'   optional (default \code{NULL}). Relabels donors at import time, useful
#'   since Vireo assigns arbitrary labels (\code{donor0}, \code{donor1}, ...),
#'   so applying the map here keeps every donor-keyed table consistent from
#'   the start.
#' @param barcode_column Character scalar (default \code{"barcode"}). Name
#'   of the column in \code{vdj_file} containing cell barcodes (only used if
#'   \code{vdj_file} is provided).
#' @param clonotype_column Character scalar (default \code{"raw_clonotype_id"}).
#'   Name of the column in \code{vdj_file} containing clonotype information
#'   (only used if \code{vdj_file} is provided).
#' @param bam_files Character vector, optional (default \code{NULL}, recording
#'   no paths). The BAM file or files this cellSNP run was made from, stored
#'   against \code{library_id} in \code{library_info} so that
#'   \code{\link{add_molecule_phase}} can find them without being told again.
#'   Paths are kept as given and only checked when read.
#'
#' @return A SNPData object
#' @family import and export functions
#' @export
#'
#' @examples
#' \dontrun{
#' # Import with VDJ and Vireo data
#' snp_data <- import_cellsnp(
#'   cellsnp_dir = "path/to/cellsnp_output",
#'   gene_annotation = gene_anno_df,
#'   vdj_file = "path/to/filtered_contig_annotations.csv",
#'   vireo_folder = "path/to/vireo_output",
#'   library_id = "run1"
#' )
#'
#' # Import without VDJ data (no clonotype information)
#' snp_data <- import_cellsnp(
#'   cellsnp_dir = "path/to/cellsnp_output",
#'   gene_annotation = gene_anno_df,
#'   library_id = "run1"
#' )
#'
#' # A Vireo output folder containing GT_donors.vireo.vcf.gz alongside
#' # donor_ids.tsv also populates per-donor zygosity
#' snp_data <- import_cellsnp(
#'   cellsnp_dir = "path/to/cellsnp_output",
#'   gene_annotation = gene_anno_df,
#'   vireo_folder = "path/to/vireo_output",
#'   library_id = "run1"
#' )
#'
#' # Relabel Vireo's arbitrary donor0/donor1 to real identities at import time
#' snp_data <- import_cellsnp(
#'   cellsnp_dir = "path/to/cellsnp_output",
#'   gene_annotation = gene_anno_df,
#'   vireo_folder = "path/to/vireo_output",
#'   donor_map = c(PatientA = "donor0", PatientB = "donor1"),
#'   library_id = "run1"
#' )
#'
#' # Two libraries labelled distinctly, so merge_snpdata() keeps cells that
#' # share a barcode by chance apart instead of fusing them
#' run1 <- import_cellsnp("run1/", gene_anno_df, library_id = "run1")
#' run2 <- import_cellsnp("run2/", gene_anno_df, library_id = "run2")
#' combined <- merge_snpdata(run1, run2)
#' }
import_cellsnp <- function(
    cellsnp_dir,
    gene_annotation,
    library_id,
    vdj_file = NULL,
    vireo_folder = NULL,
    donor_map = NULL,
    barcode_column = "barcode",
    clonotype_column = "raw_clonotype_id",
    bam_files = NULL
) {
    # Validate gene_annotation columns
    required_gene_cols <- c("chrom", "start", "end", "gene_name")
    missing_cols <- setdiff(required_gene_cols, colnames(gene_annotation))
    if (length(missing_cols) > 0) {
        stop(
            sprintf(
                "gene_annotation is missing required columns: %s",
                paste(missing_cols, collapse = ", ")
            )
        )
    }

    # Required rather than defaulted: merge_snpdata() distinguishes a barcode
    # shared within a library from one shared across libraries, and nothing in
    # the cellSNP output records which library a run came from, so a default
    # would silently make every object unmergeable.
    if (missing(library_id)) {
        stop(
            "library_id is required: name the sequencing library this cellSNP run came from, ",
            "e.g. import_cellsnp(..., library_id = \"run1\"). ",
            "merge_snpdata() needs it to tell a repeated cell from two cells sharing a barcode."
        )
    }
    if (length(library_id) != 1 || is.na(library_id)) {
        stop("library_id must be a single non-NA string naming the library this cellSNP run came from.")
    }

    # Check if required files exist
    dp_file <- fs::path(cellsnp_dir, "cellSNP.tag.DP.mtx")
    ad_file <- fs::path(cellsnp_dir, "cellSNP.tag.AD.mtx")
    oth_file <- fs::path(cellsnp_dir, "cellSNP.tag.OTH.mtx")
    base_file <- fs::path(cellsnp_dir, "cellSNP.base.vcf.gz")
    samples_file <- fs::path(cellsnp_dir, "cellSNP.samples.tsv")

    for (file in c(dp_file, ad_file, oth_file, base_file)) {
        check_file(file)
    }
    # Check optional files if provided
    if (!is.null(vdj_file)) {
        check_file(vdj_file)
    }
    # Check BAM files if provided
    if (!is.null(bam_files)) {
        for (bam_file in bam_files) {
            check_file(bam_file)
        }
    }

    # donor_ids.tsv is the point of pointing at a Vireo folder, so it must exist;
    # the genotype VCF is a bonus feature of that same run and is silently skipped
    # if the folder doesn't have one (e.g. an older or genotype-free Vireo run).
    vireo_file <- NULL
    gt_file <- NULL
    if (!is.null(vireo_folder)) {
        vireo_file <- fs::path(vireo_folder, "donor_ids.tsv")
        check_file(vireo_file)
        candidate_gt_file <- fs::path(vireo_folder, "GT_donors.vireo.vcf.gz")
        if (fs::file_exists(candidate_gt_file)) {
            check_file(candidate_gt_file)
            gt_file <- candidate_gt_file
        } else {
            logger::log_warn("Vireo folder does not contain GT_donors.vireo.vcf.gz; skipping genotype import")
        }
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
    snp_info_full <- add_snp_gene_names(snp_vcf_data, gene_annotation) %>%
        dplyr::select(snp_id, chrom, pos, ref, alt, gene_name)

    # Identify first occurrence of each unique SNP in the original VCF
    # This handles duplicates from both the VCF file and gene annotation overlaps
    keep_rows <- !duplicated(snp_vcf_data$snp_id)

    # Subset matrices to match deduplicated SNP info
    coverage <- coverage[keep_rows, , drop = FALSE]
    alt_count <- alt_count[keep_rows, , drop = FALSE]
    oth_count <- oth_count[keep_rows, , drop = FALSE]
    ref_count <- ref_count[keep_rows, , drop = FALSE]

    snp_info <- snp_info_full[keep_rows, , drop = FALSE]

    # Read donor information if provided, else create dummy donor info
    if (!is.null(vireo_file)) {
        donor_info <- readr::read_tsv(
            vireo_file,
            col_types = readr::cols(.default = readr::col_character())
        )
    } else {
        # get barcodes from cell_snp output
        donor_info <- cells %>%
            dplyr::mutate(donor = "donor0")
    }

    # Read VDJ clonotype information if provided
    if (!is.null(vdj_file)) {
        vdj_info <- readr::read_csv(
            vdj_file,
            col_types = readr::cols(.default = readr::col_character())
        ) %>%
            dplyr::mutate(
                barcode = stringr::str_remove(barcode, "-[0-9]+$") # Remove suffix if present
            )
    } else {
        vdj_info <- NULL
    }

    # Merge donor and clonotype information
    barcode_info <- merge_cell_annotations(
        donor_info,
        vdj_info,
        barcode_column,
        clonotype_column
    )

    # One cellSNP-lite run covers one library, so the label is constant across cells.
    barcode_info$library_id <- as.character(library_id)

    # Read Vireo genotype calls, if provided, to populate per-(SNP, donor)
    # zygosity at construction time
    donor_snp_info <- if (!is.null(gt_file)) {
        .read_vireo_gt(gt_file, snp_info)
    } else {
        NULL
    }

    # Create SNPData object
    logger::log_info("Creating SNPData object with {nrow(barcode_info)} barcodes and {nrow(snp_info)} SNPs")
    snp_data <- SNPData(
        alt_count = alt_count,
        ref_count = ref_count,
        oth_count = oth_count,
        snp_info = snp_info,
        barcode_info = barcode_info,
        donor_snp_info = donor_snp_info,
        donor_map = donor_map
    )

    # Import is when a BAM path is actually known -- this cellSNP run was made
    # from it -- so recording it here means add_molecule_phase() never has to
    # be told again, and the path survives every later merge.
    if (!is.null(bam_files)) {
        snp_data <- add_library_bams(snp_data, stats::setNames(list(bam_files), library_id))
    }

    return(snp_data)
}

#' Read Vireo genotype calls and classify per-donor zygosity
#'
#' Parses a Vireo genotype VCF (\code{GT_donors.vireo.vcf.gz}, one sample column per
#' donor) into the long \code{donor_snp_info} shape: GT calls of \code{"0/1"}/\code{"1/0"}
#' are classified \code{"het"}, \code{"0/0"}/\code{"1/1"} are classified \code{"hom"}.
#' When a \code{PL} (phred-scaled genotype likelihood) field is present, the posterior
#' probability of the called genotype is derived from it
#' (\code{10^(-PL/10)} normalised across the three genotypes) and stored as
#' \code{zygosity_gt_prob}.
#'
#' @param gt_file Path to a Vireo genotype VCF
#' @param snp_info The snp_info being constructed for this import, used to restrict
#'   calls to SNPs actually present (matched on the \code{chrom:pos:ref:alt} snp_id key)
#'
#' @return A tibble with columns \code{snp_id}, \code{donor}, \code{zygosity},
#'   \code{zygosity_source}, \code{zygosity_gt_prob}
#' @keywords internal
.read_vireo_gt <- function(gt_file, snp_info) {
    vcf <- read_vcf(gt_file)
    variants <- variants(vcf)
    donors <- samples(vcf)

    variants$snp_id <- make_snp_id(variants$CHROM, variants$POS, variants$REF, variants$ALT)

    format_fields <- strsplit(variants$FORMAT[1], ":")[[1]]
    gt_idx <- match("GT", format_fields)
    pl_idx <- match("PL", format_fields)
    if (is.na(gt_idx)) {
        stop("Vireo GT VCF has no GT field in its FORMAT column")
    }

    donor_calls <- purrr::map(donors, function(d) {
        fields <- stringr::str_split_fixed(variants[[d]], ":", length(format_fields))
        gt <- fields[, gt_idx]
        zygosity <- dplyr::case_when(
            gt %in% c("0/1", "1/0") ~ "het",
            gt %in% c("0/0", "1/1") ~ "hom",
            TRUE ~ NA_character_
        )

        zygosity_gt_prob <- rep(NA_real_, length(gt))
        if (!is.na(pl_idx)) {
            pl <- stringr::str_split_fixed(fields[, pl_idx], ",", 3)
            # A missing PL is written as "." by Vireo; coerces to NA as intended.
            suppressWarnings(storage.mode(pl) <- "double")
            gt_index <- dplyr::case_when(
                gt == "0/0" ~ 1L,
                gt %in% c("0/1", "1/0") ~ 2L,
                gt == "1/1" ~ 3L,
                TRUE ~ NA_integer_
            )
            has_pl <- !is.na(gt_index) & !purrr::map_lgl(asplit(pl, 1), anyNA)
            likelihood <- 10^(-pl[has_pl, , drop = FALSE] / 10)
            posterior <- likelihood / rowSums(likelihood)
            zygosity_gt_prob[has_pl] <- posterior[cbind(seq_len(sum(has_pl)), gt_index[has_pl])]
        }

        tibble::tibble(
            snp_id = variants$snp_id,
            donor = d,
            zygosity = zygosity,
            zygosity_source = ifelse(is.na(zygosity), NA_character_, "vireo_gt"),
            zygosity_gt_prob = zygosity_gt_prob
        )
    }) %>%
        dplyr::bind_rows()

    donor_calls %>%
        dplyr::filter(snp_id %in% snp_info$snp_id, !is.na(zygosity))
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

    # Generate standardised SNP IDs
    vcf_data$snp_id <- make_snp_id(
        vcf_data$chrom,
        vcf_data$pos,
        vcf_data$ref,
        vcf_data$alt
    )

    # Reorder columns to have snp_id first
    vcf_data <- vcf_data[, c("snp_id", names(vcf_data)[names(vcf_data) != "snp_id"])]

    return(vcf_data)
}

#' Merge donor and clonotype information
#'
#' @param donor_info Data frame with donor information from Vireo
#' @param vdj_info Data frame with VDJ information from cellranger (NULL if not provided)
#' @param barcode_column Name of the column in vdj_info containing cell barcodes (only used if vdj_info provided)
#' @param clonotype_column Name of the column in vdj_info containing clonotype information (only used if vdj_info provided)
#'
#' @return Data frame with merged cell annotations
#' @keywords internal
merge_cell_annotations <- function(donor_info, vdj_info = NULL, barcode_column = NULL, clonotype_column = NULL) {
    # Standardise column names
    if ("cell" %in% colnames(donor_info)) {
        donor_info <- donor_info %>%
            dplyr::rename(barcode = cell)
    } else if ("cell_id" %in% colnames(donor_info)) {
        donor_info <- donor_info %>%
            dplyr::rename(barcode = cell_id)
    }

    if ("donor_id" %in% colnames(donor_info)) {
        donor_info <- donor_info %>%
            dplyr::rename(donor = donor_id)
    }

    # If VDJ info not provided, create barcode_info without clonotype
    if (is.null(vdj_info)) {
        barcode_info <- donor_info %>%
            dplyr::mutate(
                cell_id = paste0("cell_", seq_len(dplyr::n())),
                clonotype = NA_character_
            ) %>%
            dplyr::select(cell_id, barcode, donor, clonotype, everything())

        return(barcode_info)
    }

    # Ensure barcode_column exists in vdj_info
    if (!barcode_column %in% colnames(vdj_info)) {
        stop(paste0(
            "Column ",
            barcode_column,
            " not found in VDJ annotation file"
        ))
    }

    # Ensure clonotype_column exists in vdj_info
    if (!clonotype_column %in% colnames(vdj_info)) {
        stop(paste0(
            "Column ",
            clonotype_column,
            " not found in VDJ annotation file"
        ))
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
#' @family import and export functions
#' @export
#' @examples
#' snp_data <- get_example_snpdata()
get_example_snpdata <- function() {
    import_cellsnp(
        cellsnp_dir = system.file("extdata/example_snpdata", package = "snplet"),
        gene_annotation = readr::read_tsv(
            system.file("extdata/example_gene_anno.tsv", package = "snplet"),
            show_col_types = FALSE
        ),
        vdj_file = system.file("extdata/example_snpdata/filtered_contig_annotations.csv", package = "snplet"),
        vireo_folder = system.file("extdata/example_snpdata", package = "snplet"),
        library_id = "example"
    )
}

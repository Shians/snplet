#' Import cellSNP data and create a SNPData object
#'
#' This function imports data from cellSNP-lite output, with optional VDJ annotations from cellranger
#' and donor information from Vireo to create a SNPData object.
#'
#' @param cellsnp_dir Character scalar, required. Directory containing
#'   cellSNP-lite output files.
#' @param gene_annotation A data.frame, required, with columns \code{chrom},
#'   \code{start}, \code{end}, \code{gene_name}. Gene annotations.
#' @param library_id Character scalar, optional (default \code{NA}). Name of
#'   the sequencing library this cellSNP-lite run came from, stored on every
#'   cell. See \sQuote{Merging libraries} below.
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
#'
#' @return A SNPData object
#'
#' @section Merging libraries:
#' To import several runs or libraries into one object, use
#' \code{\link{import_cellsnp_libraries}}. A 10x barcode is unique only within
#' its library, so that function matches cells on (\code{library_id},
#' \code{barcode}) and requires a \code{library_id} for every run. A
#' single-library import can safely leave \code{library_id} unset.
#'
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
#' }
import_cellsnp <- function(
    cellsnp_dir,
    gene_annotation,
    library_id = NA_character_,
    vdj_file = NULL,
    vireo_folder = NULL,
    donor_map = NULL,
    barcode_column = "barcode",
    clonotype_column = "raw_clonotype_id"
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

    # Left NA rather than defaulted to a guessed label: nothing in the cellSNP
    # output records which library a run came from, and a guessed default
    # (e.g. the directory name) risks two different libraries colliding
    # silently at merge time. import_cellsnp_libraries() requires a
    # library_id on every row, so a single-library import can safely leave
    # this unset, and a multi-library import is stopped there instead.
    if (length(library_id) != 1) {
        stop("library_id must be a single string naming the library this cellSNP run came from, or NA.")
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

    # Read cellSNP matrices. read_mtx() parses the coordinate triples with
    # readr rather than Matrix::readMM()'s scan()-based reader, which measurably
    # speeds up the multi-hundred-megabyte matrices cellSNP-lite produces.
    coverage <- read_mtx(dp_file)
    alt_count <- read_mtx(ad_file)
    oth_count <- read_mtx(oth_file)
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
        donor_map = donor_map,
        total_count = coverage
    )

    return(snp_data)
}

#' Import several cellSNP-lite runs into one SNPData object
#'
#' Imports each cellSNP-lite run listed in a sample sheet with
#' \code{\link{import_cellsnp}} and combines them into a single SNPData object.
#' Runs from different libraries contribute separate cells; runs that share a
#' \code{library_id} are treated as repeat runs over the same cells, and their
#' counts are summed (see \sQuote{Repeat runs of a library}).
#'
#' @param sheet A data.frame, required, with one row per cellSNP-lite run. Must
#'   contain the columns \code{cellsnp_dir} and \code{library_id} (character,
#'   non-\code{NA}), and may contain \code{vdj_file}, \code{vireo_folder},
#'   \code{bam_files}, and \code{donor_map}. Each column except
#'   \code{bam_files} is passed to the \code{\link{import_cellsnp}} argument of
#'   the same name for that row; see \sQuote{Sample sheet} below.
#' @param gene_annotation A data.frame, required, with columns \code{chrom},
#'   \code{start}, \code{end}, \code{gene_name}. Gene annotations, shared by
#'   every run.
#' @param barcode_column Character scalar (default \code{"barcode"}). Name of
#'   the barcode column in every run's \code{vdj_file}.
#' @param clonotype_column Character scalar (default
#'   \code{"raw_clonotype_id"}). Name of the clonotype column in every run's
#'   \code{vdj_file}.
#'
#' @return A SNPData object holding the cells of every run. Its SNPs are the
#'   union of all runs' SNPs, zero-filled where a run did not measure a SNP.
#'
#' @section Sample sheet:
#' \describe{
#'   \item{\code{cellsnp_dir}}{Required. Directory of cellSNP-lite output. Each
#'     directory may appear only once, since importing it twice would double
#'     its counts.}
#'   \item{\code{library_id}}{Required. Name of the sequencing library the run
#'     came from. A 10x barcode is unique only within its library, so cells
#'     are matched on (\code{library_id}, \code{barcode}): a barcode repeated
#'     across rows with the same \code{library_id} is one cell, and one
#'     repeated across different libraries is two.}
#'   \item{\code{vdj_file}, \code{vireo_folder}}{Optional character columns;
#'     \code{NA} means none for that row.}
#'   \item{\code{bam_files}}{Optional. A character column (one BAM per row) or a
#'     list column of character vectors (several BAMs per row); \code{NA} or
#'     \code{NULL} means none. Used only to check that repeat runs of a
#'     library do not count the same reads twice (see \sQuote{Repeat runs of a
#'     library}); the paths are not stored on the returned object. Each path is
#'     checked to exist.}
#'   \item{\code{donor_map}}{Optional list column of named character vectors,
#'     \code{c(new_label = old_label, ...)}, relabelling that row's donors;
#'     \code{NULL} means no relabelling.}
#' }
#'
#' @section Repeat runs of a library:
#' Summing the counts of runs that share a \code{library_id} is correct only
#' when no read is counted by both: either the runs come from different BAM
#' files (the library sequenced twice), or they come from the same BAM but
#' cover disjoint SNPs (e.g. one run per chromosome). Two runs over the same
#' BAM that share SNPs count those reads twice, which breaks the independence
#' assumed by the binomial and beta-binomial tests. The function therefore:
#' \itemize{
#'   \item errors when two runs of a library share a \code{bam_files} path
#'     and at least one SNP;
#'   \item warns when two runs of a library share SNPs but either lacks
#'     \code{bam_files}, since the two cases cannot then be told apart;
#'   \item errors when two runs of a library assign the same cell to
#'     different donors, as can happen when each run has its own
#'     \code{vireo_folder}, since Vireo labels donors independently per run.
#' }
#'
#' @section Donor labels:
#' Donors must not span libraries. Vireo numbers donors from \code{donor0} in
#' every run, and a run imported without \code{vireo_folder} labels all its
#' cells \code{donor0}, so the same label in two libraries almost always names
#' two different people. The function therefore errors when a donor label
#' (other than \code{"doublet"} or \code{"unassigned"}) appears in more than
#' one library; give each library distinct labels through \code{donor_map}.
#'
#' @family import and export functions
#' @export
#'
#' @examples
#' \dontrun{
#' # Two libraries, the first sequenced twice (so two BAMs holding different
#' # reads). Vireo labels donors from donor0 in each library, so donor_map
#' # gives every library its own labels.
#' sheet <- tibble::tibble(
#'   cellsnp_dir = c("lib1_run1/", "lib1_run2/", "lib2/"),
#'   library_id = c("lib1", "lib1", "lib2"),
#'   vireo_folder = c("vireo_lib1/", "vireo_lib1/", "vireo_lib2/"),
#'   bam_files = c("lib1_seq1.bam", "lib1_seq2.bam", "lib2.bam"),
#'   donor_map = list(
#'     c(PatientA = "donor0", PatientB = "donor1"),
#'     c(PatientA = "donor0", PatientB = "donor1"),
#'     c(PatientC = "donor0", PatientD = "donor1")
#'   )
#' )
#' snp_data <- import_cellsnp_libraries(sheet, gene_anno_df)
#' }
import_cellsnp_libraries <- function(
    sheet,
    gene_annotation,
    barcode_column = "barcode",
    clonotype_column = "raw_clonotype_id"
) {
    sheet <- .validate_import_sheet(sheet)

    logger::log_info(
        "Importing {nrow(sheet)} cellSNP run(s) from {dplyr::n_distinct(sheet$library_id)} library(ies)"
    )
    runs <- purrr::pmap(sheet, function(cellsnp_dir, library_id, vdj_file, vireo_folder, bam_files, donor_map) {
        import_cellsnp(
            cellsnp_dir = cellsnp_dir,
            gene_annotation = gene_annotation,
            library_id = library_id,
            vdj_file = vdj_file,
            vireo_folder = vireo_folder,
            donor_map = donor_map,
            barcode_column = barcode_column,
            clonotype_column = clonotype_column
        )
    })

    .check_repeat_run_overlap(sheet, runs)
    .check_repeat_run_donors(runs)
    .check_donors_within_library(runs)

    # merge_snpdata() defaults to union joins on both axes, which is what a
    # sample sheet means: keep every run's SNPs and cells, summing only where
    # two runs share a (library_id, barcode) cell.
    Reduce(merge_snpdata, runs)
}

# Columns a sample sheet may carry, in the order .validate_import_sheet()
# returns them and import_cellsnp_libraries() unpacks them.
.IMPORT_SHEET_COLUMNS <- c("cellsnp_dir", "library_id", "vdj_file", "vireo_folder", "bam_files", "donor_map")

# Checks a sample sheet and normalises it to a tibble with every column in
# .IMPORT_SHEET_COLUMNS. Optional columns become list columns holding NULL
# where a row has no value, so each row unpacks straight into import_cellsnp()
# arguments. Unknown columns are rejected rather than ignored, so a misspelt
# optional column cannot silently drop its data.
.validate_import_sheet <- function(sheet) {
    if (!is.data.frame(sheet) || nrow(sheet) == 0) {
        stop("sheet must be a data.frame with one row per cellSNP-lite run.")
    }

    missing_cols <- setdiff(c("cellsnp_dir", "library_id"), colnames(sheet))
    if (length(missing_cols) > 0) {
        stop("sheet is missing required column(s): ", paste(missing_cols, collapse = ", "))
    }

    unknown_cols <- setdiff(colnames(sheet), .IMPORT_SHEET_COLUMNS)
    if (length(unknown_cols) > 0) {
        stop(
            "sheet has unrecognised column(s): ",
            paste(unknown_cols, collapse = ", "),
            ". Allowed columns are: ",
            paste(.IMPORT_SHEET_COLUMNS, collapse = ", ")
        )
    }

    sheet <- tibble::as_tibble(sheet)
    sheet$cellsnp_dir <- as.character(sheet$cellsnp_dir)
    sheet$library_id <- as.character(sheet$library_id)

    if (anyNA(sheet$library_id) || any(sheet$library_id == "")) {
        stop(
            "Every row of sheet needs a library_id: cells from different runs are matched on ",
            "(library_id, barcode), so an unlabelled run cannot be combined safely."
        )
    }

    dir_keys <- as.character(fs::path_abs(sheet$cellsnp_dir))
    if (anyDuplicated(dir_keys)) {
        stop(
            "sheet lists the same cellsnp_dir more than once, which would double its counts: ",
            paste(unique(sheet$cellsnp_dir[duplicated(dir_keys)]), collapse = ", ")
        )
    }

    for (col in setdiff(.IMPORT_SHEET_COLUMNS, c("cellsnp_dir", "library_id"))) {
        sheet[[col]] <- .as_sheet_list_column(sheet[[col]], nrow(sheet))
    }
    for (bam_file in unlist(sheet$bam_files)) {
        check_file(bam_file)
    }

    sheet[.IMPORT_SHEET_COLUMNS]
}

# An optional sheet column as a list, one element per row, with NULL for an
# absent column or an NA entry: import_cellsnp() reads NULL, not NA, as "none".
.as_sheet_list_column <- function(values, n_rows) {
    if (is.null(values)) {
        return(vector("list", n_rows))
    }
    purrr::map(as.list(values), .na_to_null)
}

.na_to_null <- function(value) {
    if (length(value) == 1 && is.na(value)) {
        return(NULL)
    }
    value
}

# Runs of one library are summed cell by cell, which is only correct when no
# read is counted by both. Two runs over different BAMs (the library sequenced
# twice) hold different reads; two runs over the same BAM hold the same reads, so
# they may only be summed if they pile up different SNPs (e.g. split by
# chromosome). A shared BAM with shared SNPs double counts every read at those
# SNPs, breaking the independence the binomial and beta-binomial tests assume,
# so it errors. Without BAM paths the two cases look identical, so shared SNPs
# only warn.
.check_repeat_run_overlap <- function(sheet, runs) {
    run_pairs <- .same_library_run_pairs(sheet$library_id)
    if (nrow(run_pairs) == 0) {
        return(invisible(NULL))
    }

    run_pairs$n_shared_snps <- purrr::map2_int(run_pairs$first, run_pairs$second, function(first, second) {
        length(intersect(snp_info(runs[[first]])$snp_id, snp_info(runs[[second]])$snp_id))
    })
    run_pairs$bams_known <- purrr::map2_lgl(run_pairs$first, run_pairs$second, function(first, second) {
        !is.null(sheet$bam_files[[first]]) && !is.null(sheet$bam_files[[second]])
    })
    run_pairs$shares_bam <- purrr::map2_lgl(run_pairs$first, run_pairs$second, function(first, second) {
        shared <- intersect(.bam_keys(sheet$bam_files[[first]]), .bam_keys(sheet$bam_files[[second]]))
        length(shared) > 0
    })
    run_pairs$label <- paste0(
        sheet$cellsnp_dir[run_pairs$first],
        " and ",
        sheet$cellsnp_dir[run_pairs$second],
        " (",
        run_pairs$n_shared_snps,
        " shared SNPs)"
    )

    double_counted <- run_pairs[run_pairs$shares_bam & run_pairs$n_shared_snps > 0, ]
    if (nrow(double_counted) > 0) {
        stop(
            "Runs of the same library share a BAM file and SNPs, so summing them would count ",
            "the same reads twice: ",
            paste(double_counted$label, collapse = "; "),
            ". Runs over one BAM must cover disjoint SNPs, e.g. one run per chromosome."
        )
    }

    unverifiable <- run_pairs[!run_pairs$bams_known & run_pairs$n_shared_snps > 0, ]
    if (nrow(unverifiable) > 0) {
        warning(
            "Runs of the same library share SNPs, and their counts will be summed: ",
            paste(unverifiable$label, collapse = "; "),
            ". This is correct only if the runs were made from different BAM files; if they ",
            "share a BAM, the same reads are counted twice. Add a bam_files column to sheet ",
            "to have this checked.",
            call. = FALSE
        )
    }

    invisible(NULL)
}

# Every unordered pair of sheet rows that share a library_id, as row indices.
.same_library_run_pairs <- function(library_id) {
    same_library <- outer(library_id, library_id, "==") & upper.tri(diag(length(library_id)))
    pairs <- which(same_library, arr.ind = TRUE)
    tibble::tibble(first = unname(pairs[, 1]), second = unname(pairs[, 2]))
}

# A BAM path in a form comparable across rows, so the same file written two
# ways (relative vs absolute, via a symlink) is still recognised.
# .validate_import_sheet() has already checked that each path exists.
.bam_keys <- function(bam_files) {
    if (is.null(bam_files)) {
        return(character(0))
    }
    as.character(fs::path_real(bam_files))
}

# A cell seen by more than one run of its library is one cell, so every run must
# give it the same donor. Vireo labels each run independently, so two runs with
# their own vireo_folder can call the same person donor0 in one and donor1 in the
# other; merging would then keep the first run's label and attach the second
# run's genotype calls to the wrong person. Checked after each row's donor_map,
# since it is the final labels that must agree.
.check_repeat_run_donors <- function(runs) {
    cell_donors <- purrr::map(runs, function(run) {
        dplyr::select(barcode_info(run), library_id, barcode, donor)
    }) %>%
        dplyr::bind_rows() %>%
        dplyr::distinct()

    # After distinct(), a (library_id, barcode) key repeats only where its runs
    # disagree on the donor.
    conflicting <- duplicated(cell_donors[c("library_id", "barcode")])
    if (!any(conflicting)) {
        return(invisible(NULL))
    }

    conflict_keys <- unique(cell_donors[conflicting, c("library_id", "barcode")])
    examples <- utils::head(conflict_keys, 5)
    example_labels <- purrr::map2_chr(examples$library_id, examples$barcode, function(lib, bc) {
        donors <- cell_donors$donor[cell_donors$library_id == lib & cell_donors$barcode == bc]
        paste0(lib, "/", bc, " (", paste(donors, collapse = " vs "), ")")
    })

    stop(
        nrow(conflict_keys),
        " cell(s) are assigned different donors by different runs of the same library, e.g. ",
        paste(example_labels, collapse = ", "),
        ". Vireo labels donors independently in each run; use one vireo_folder for every run ",
        "of a library, or align the labels with donor_map."
    )
}

# Donors must not span libraries (a library's BAM holds all of its donors'
# cells, and molecule phasing looks donors up by library). Checked on the
# imported runs, after each row's donor_map, because it is the final labels that
# must be distinct. "doublet" and "unassigned" are Vireo's shared status labels
# rather than donors, so they may appear in every library.
.check_donors_within_library <- function(runs) {
    donor_libraries <- purrr::map(runs, function(run) {
        dplyr::distinct(barcode_info(run), library_id, donor)
    }) %>%
        dplyr::bind_rows() %>%
        dplyr::distinct() %>%
        dplyr::filter(!is.na(donor), !donor %in% c("doublet", "unassigned"))

    shared_donors <- unique(donor_libraries$donor[duplicated(donor_libraries$donor)])
    if (length(shared_donors) == 0) {
        return(invisible(NULL))
    }

    stop(
        "Donor label(s) found in more than one library: ",
        paste(shared_donors, collapse = ", "),
        ". Vireo numbers donors from donor0 in every library, and a run without vireo_folder labels ",
        "all its cells donor0, so a shared label usually names different people. Give each library ",
        "distinct labels with a donor_map list column in sheet."
    )
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

#' Read a MatrixMarket file into a sparse Matrix
#'
#' A faster drop-in for \code{Matrix::readMM()} on the plain-text
#' \code{.mtx} files cellSNP-lite produces: parses the coordinate entries
#' with \code{readr::read_table()} (a vectorised C++ parser) rather than the
#' \code{scan()}-based reader that \code{Matrix::readMM()} uses, which
#' measurably speeds up import on multi-hundred-megabyte matrices.
#'
#' Accepts the variations the MatrixMarket format permits: fields separated
#' by any run of spaces or tabs, leading and trailing whitespace, CRLF line
#' endings, blank lines, any number of comment lines, gzip compression, and
#' values in scientific notation. Both the sparse \code{coordinate} and dense
#' \code{array} formats are read, with the \code{real}, \code{integer}, or
#' (coordinate only) \code{pattern} field and the \code{general},
#' \code{symmetric}, or \code{skew-symmetric} symmetry; the stored triangle
#' of a symmetric matrix is mirrored into a full general matrix. A file
#' without a banner line is read as \code{coordinate real general}. Complex
#' and Hermitian matrices are rejected.
#'
#' @param mtx_file Path to a \code{.mtx} or \code{.mtx.gz} file
#'
#' @return A sparse \code{dgCMatrix}
#' @keywords internal
read_mtx <- function(mtx_file) {
    header <- read_mtx_header(mtx_file)

    if (header$format == "array") {
        entries <- read_mtx_array_entries(mtx_file, header)
    } else {
        entries <- read_mtx_coordinate_entries(mtx_file, header)
    }

    row_index <- entries$i
    col_index <- entries$j
    values <- entries$x

    # Symmetric files store only one triangle, so mirror the off-diagonal
    # entries (negated for skew-symmetric) to recover the full matrix.
    if (header$symmetry != "general") {
        is_off_diagonal <- row_index != col_index
        mirror_sign <- 1
        if (header$symmetry == "skew-symmetric") {
            mirror_sign <- -1
        }
        mirrored_rows <- col_index[is_off_diagonal]
        col_index <- c(col_index, row_index[is_off_diagonal])
        row_index <- c(row_index, mirrored_rows)
        values <- c(values, mirror_sign * values[is_off_diagonal])
    }

    Matrix::sparseMatrix(
        i = row_index,
        j = col_index,
        x = values,
        dims = header$dims,
        index1 = TRUE
    )
}

#' Read the entries of a coordinate-format MatrixMarket file
#'
#' @param mtx_file Path to a \code{.mtx} or \code{.mtx.gz} file
#' @param header The list returned by \code{read_mtx_header()}
#'
#' @return A list of 1-based row indices \code{i}, column indices \code{j},
#'   and values \code{x} (all ones for a pattern file), as stored in the file.
#' @keywords internal
read_mtx_coordinate_entries <- function(mtx_file, header) {
    value_cols <- c("i", "j", "x")
    if (header$field == "pattern") {
        value_cols <- c("i", "j")
    }

    # A line with a missing or surplus field parses to NA in some column (the
    # read_table() fallback puts surplus fields in its extra column), so
    # readr's parsing warnings are silenced and malformed lines found here.
    find_malformed <- function(entries) {
        # read_delim() sizes the table from the first line, so a malformed
        # first line can leave whole columns missing.
        if (!all(value_cols %in% names(entries))) {
            return(rep(TRUE, nrow(entries)))
        }
        is_malformed <- Reduce(`|`, lapply(entries[value_cols], is.na))
        if ("extra" %in% names(entries)) {
            is_malformed <- is_malformed | !is.na(entries$extra)
        }
        is_malformed
    }

    # Fast path: read_delim() is multithreaded but needs exactly one
    # separator character between fields, so guess it from the first entry.
    # Indices are parsed as integers to save memory, so an index written in
    # scientific notation also falls through to the slow path.
    delim <- " "
    if (stringr::str_detect(header$first_entry, "\t")) {
        delim <- "\t"
    }
    entries <- suppressWarnings(readr::read_delim(
        mtx_file,
        delim = delim,
        skip = header$n_lines,
        col_names = value_cols,
        col_types = substr("iid", 1, length(value_cols)),
        progress = FALSE
    ))

    # Slow path: read_table() splits on any run of whitespace, which handles
    # mixed separators, repeated separators, and padded lines.
    is_malformed <- find_malformed(entries)
    if (any(is_malformed)) {
        entries <- suppressWarnings(readr::read_table(
            mtx_file,
            skip = header$n_lines,
            col_names = c(value_cols, "extra"),
            col_types = paste0(strrep("d", length(value_cols)), "c"),
            progress = FALSE
        ))
        is_malformed <- find_malformed(entries)
    }

    if (nrow(entries) != header$nnz) {
        stop(
            "MatrixMarket file ",
            mtx_file,
            " declares ",
            header$nnz,
            " entries but contains ",
            nrow(entries),
            "; the file may be truncated."
        )
    }

    if (any(is_malformed)) {
        stop(
            "MatrixMarket file ",
            mtx_file,
            " has a malformed entry at entry ",
            which(is_malformed)[1],
            "; expected ",
            length(value_cols),
            " whitespace-separated numbers per line."
        )
    }

    # Check the index range cheaply first; only locate the offending entry
    # when there is one to report.
    has_invalid_index <- function(index, n) {
        index_range <- range(index, 1)
        index_range[1] < 1 || index_range[2] > n || (!is.integer(index) && any(index != round(index)))
    }
    if (has_invalid_index(entries$i, header$dims[1]) || has_invalid_index(entries$j, header$dims[2])) {
        is_valid_index <- function(index, n) {
            index == round(index) & index >= 1 & index <= n
        }
        is_out_of_range <- !is_valid_index(entries$i, header$dims[1]) | !is_valid_index(entries$j, header$dims[2])
        stop(
            "MatrixMarket file ",
            mtx_file,
            " has an invalid index at entry ",
            which(is_out_of_range)[1],
            "; indices must be whole numbers within the declared ",
            header$dims[1],
            " x ",
            header$dims[2],
            " dimensions."
        )
    }

    values <- rep(1, nrow(entries))
    if (header$field != "pattern") {
        values <- entries$x
    }
    list(i = entries$i, j = entries$j, x = values)
}

#' Read the entries of an array-format MatrixMarket file
#'
#' Array files list values in column-major order with no indices: every
#' element of a general matrix, the lower triangle including the diagonal of
#' a symmetric matrix, or the strict lower triangle of a skew-symmetric one.
#'
#' @inheritParams read_mtx_coordinate_entries
#'
#' @return A list of 1-based row indices \code{i}, column indices \code{j},
#'   and values \code{x} for the non-zero stored elements.
#' @keywords internal
read_mtx_array_entries <- function(mtx_file, header) {
    if (header$field == "pattern") {
        stop("MatrixMarket file ", mtx_file, " uses the pattern field, which the array format does not allow.")
    }

    values <- suppressWarnings(as.numeric(scan(
        mtx_file,
        what = "character",
        skip = header$n_lines,
        quiet = TRUE
    )))

    stored_positions <- matrix(TRUE, nrow = header$dims[1], ncol = header$dims[2])
    if (header$symmetry != "general") {
        stored_positions <- lower.tri(stored_positions, diag = header$symmetry == "symmetric")
    }
    positions <- which(stored_positions, arr.ind = TRUE)

    if (length(values) != nrow(positions) || anyNA(values)) {
        stop(
            "MatrixMarket file ",
            mtx_file,
            " should hold ",
            nrow(positions),
            " numeric array values but has ",
            length(values),
            " fields, ",
            sum(is.na(values)),
            " of them non-numeric."
        )
    }

    is_non_zero <- values != 0
    list(i = positions[is_non_zero, 1], j = positions[is_non_zero, 2], x = values[is_non_zero])
}

#' Parse the banner and size line of a MatrixMarket file
#'
#' @param mtx_file Path to a \code{.mtx} or \code{.mtx.gz} file
#'
#' @return A list with the banner's \code{format}, \code{field}, and
#'   \code{symmetry} (lower-cased), the matrix \code{dims}, the declared
#'   \code{nnz} (\code{NA} for array files), \code{n_lines}, the number of
#'   lines up to and including the size line, and \code{first_entry}, the
#'   first non-blank line after it (\code{""} if there is none).
#' @keywords internal
read_mtx_header <- function(mtx_file) {
    # The banner, comment lines, and blank lines precede the size line, and
    # the format allows any number of them, so read further until it is found.
    n_max <- 64L
    repeat {
        lines <- readr::read_lines(mtx_file, n_max = n_max, skip_empty_rows = FALSE)
        is_preamble <- stringr::str_detect(lines, "^\\s*(%|$)")
        size_line <- match(FALSE, is_preamble)
        if (!is.na(size_line) || length(lines) < n_max) {
            break
        }
        n_max <- n_max * 4L
    }
    if (is.na(size_line)) {
        stop("MatrixMarket file ", mtx_file, " has no size line.")
    }

    header <- list(format = "coordinate", field = "real", symmetry = "general")
    if (stringr::str_detect(lines[1], stringr::regex("^%%MatrixMarket", ignore_case = TRUE))) {
        banner <- tolower(stringr::str_split(stringr::str_trim(lines[1]), "\\s+")[[1]])
        if (length(banner) != 5 || banner[2] != "matrix") {
            stop("MatrixMarket file ", mtx_file, " has a malformed banner: ", lines[1])
        }
        header <- list(format = banner[3], field = banner[4], symmetry = banner[5])
    }

    if (!header$format %in% c("coordinate", "array")) {
        stop("MatrixMarket file ", mtx_file, " has unsupported format '", header$format, "'.")
    }
    if (!header$field %in% c("real", "integer", "pattern", "double")) {
        stop("MatrixMarket file ", mtx_file, " has unsupported field '", header$field, "'.")
    }
    if (!header$symmetry %in% c("general", "symmetric", "skew-symmetric")) {
        stop("MatrixMarket file ", mtx_file, " has unsupported symmetry '", header$symmetry, "'.")
    }

    n_size_fields <- 3L
    if (header$format == "array") {
        n_size_fields <- 2L
    }
    size <- suppressWarnings(as.numeric(stringr::str_split(stringr::str_trim(lines[size_line]), "\\s+")[[1]]))
    if (length(size) != n_size_fields || anyNA(size)) {
        stop(
            "MatrixMarket file ",
            mtx_file,
            " has a malformed size line: '",
            lines[size_line],
            "'; expected ",
            n_size_fields,
            " numbers."
        )
    }

    header$dims <- as.integer(size[1:2])
    header$nnz <- size[3]
    header$n_lines <- size_line

    following_lines <- readr::read_lines(mtx_file, skip = size_line, n_max = 16L)
    header$first_entry <- c(following_lines[following_lines != ""], "")[1]
    header
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

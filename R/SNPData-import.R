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
#' @param bam_files An unnamed character vector, optional (default \code{NULL},
#'   recording no paths). The BAM file or files this cellSNP run was made from.
#'   Requires \code{library_id} (this whole vector is stored against it as one
#'   library's paths in \code{library_info}, via \code{\link{add_library_bams}}),
#'   so that \code{\link{phase_from_molecules}} can find them without being told
#'   again. Must not be named: this call covers a single library, so any names
#'   on \code{bam_files} itself would be silently discarded rather than used as
#'   per-library keys; each path is checked to exist.
#'
#' @return A SNPData object
#'
#' @section Merging libraries:
#' A 10x barcode is unique only within its library, so
#' \code{\link{merge_snpdata}} uses \code{library_id} to tell a repeated cell
#' from two different cells that happened to draw the same barcode, and
#' refuses to merge any object carrying a \code{NA} \code{library_id}. A
#' single-library workflow that never calls \code{merge_snpdata} can safely
#' leave \code{library_id} unset.
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
    library_id = NA_character_,
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

    # Left NA rather than defaulted to a guessed label: nothing in the cellSNP
    # output records which library a run came from, and a guessed default
    # (e.g. the directory name) risks two different libraries colliding
    # silently at merge time. merge_snpdata() already refuses to merge any
    # object with an NA library_id (see .check_library_ids()), so a
    # single-library workflow that never merges can safely leave this unset,
    # and a workflow that does merge is stopped there instead.
    if (length(library_id) != 1) {
        stop("library_id must be a single string naming the library this cellSNP run came from, or NA.")
    }

    # bam_files is recorded against library_id in library_info, so there is
    # nothing to key it against when library_id was left NA. Checked ahead of
    # file existence so this points at the real cause rather than a confusing
    # add_library_bams() error once import has otherwise succeeded.
    if (!is.null(bam_files) && is.na(library_id)) {
        stop(
            "bam_files was supplied but library_id was not: BAM paths are recorded against ",
            "library_id in library_info, so set library_id = to use bam_files."
        )
    }

    # One import call covers one library, so bam_files is keyed by library_id
    # automatically below; any names on bam_files itself would be silently
    # discarded by that wrapping (add_library_bams()'s per-library keys come
    # from the outer list this constructs, not from bam_files' own names),
    # so a named vector is rejected here rather than left to fail silently.
    if (!is.null(bam_files) && !is.null(names(bam_files))) {
        stop(
            "bam_files must be an unnamed character vector of path(s) for the ",
            "single library named by library_id; names on bam_files are ignored ",
            "and would be silently dropped."
        )
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

    # Import is the one point where the gene annotation is in hand, so the
    # molecule-level SNP-to-gene map is derived here rather than asked for again
    # by haplotype_expression_by_molecule(). It needs strand, which the display
    # labels in snp_info$gene_name do not, so an unstranded annotation leaves the
    # map empty; assign_snp_genes() can supply it later via snp_gene_map<-().
    snp_gene_map <- if ("strand" %in% colnames(gene_annotation)) {
        assign_snp_genes(snp_info, gene_annotation)
    } else {
        logger::log_info(
            "gene_annotation has no strand column, so no SNP-to-gene map was built; ",
            "haplotype_expression_by_molecule() needs one, set later with snp_gene_map(x) <- ."
        )
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
        snp_gene_map = snp_gene_map,
        total_count = coverage
    )

    # Import is when a BAM path is actually known -- this cellSNP run was made
    # from it -- so recording it here means phase_from_molecules() never has to
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

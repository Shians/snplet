# Internal helpers for SNPData construction and validation

.validate_count_dims <- function(ref_count, alt_count, oth_count) {
    stopifnot(
        "ref_count and alt_count must have the same number of rows (SNPs)" = nrow(alt_count) == nrow(ref_count),
        "ref_count and alt_count must have the same number of columns (cells)" = ncol(alt_count) == ncol(ref_count)
    )

    if (is.null(oth_count)) {
        oth_count <- Matrix::Matrix(
            0,
            nrow = nrow(ref_count),
            ncol = ncol(ref_count),
            sparse = TRUE
        )
    } else {
        stopifnot(
            "oth_count must have the same number of rows (SNPs) as ref_count" = nrow(oth_count) == nrow(ref_count),
            "oth_count must have the same number of columns (cells) as ref_count" = ncol(oth_count) == ncol(ref_count)
        )
    }

    oth_count
}

.validate_total_count <- function(total_count, ref_count, alt_count) {
    stopifnot(
        "total_count must have the same number of rows (SNPs) as ref_count" = nrow(total_count) == nrow(ref_count),
        "total_count must have the same number of columns (cells) as ref_count" = ncol(total_count) == ncol(ref_count)
    )
    # A full total_count == alt_count + ref_count check would cost as much as
    # the addition it's meant to save, so this checks the row/column margins
    # instead -- catching a mismatched or stale total_count in practice
    # without materialising a second full sparse matrix.
    row_diff <- Matrix::rowSums(total_count) - Matrix::rowSums(alt_count) - Matrix::rowSums(ref_count)
    if (any(row_diff != 0)) {
        stop("total_count does not match alt_count + ref_count (row sums differ)")
    }
    total_count
}

.validate_info_dims <- function(ref_count, alt_count, snp_info, barcode_info) {
    stopifnot(
        "ncol(alt_count) must equal nrow(barcode_info)" = ncol(alt_count) == nrow(barcode_info),
        "nrow(ref_count) must equal nrow(snp_info)" = nrow(ref_count) == nrow(snp_info)
    )
}

.validate_donor_dims <- function(donor_info, donor_snp_info, snp_info) {
    problems <- character(0)

    if (!"donor" %in% colnames(donor_info)) {
        problems <- c(problems, "donor_info must contain a 'donor' column")
    } else if (any(duplicated(donor_info$donor))) {
        problems <- c(problems, "donor_info contains duplicate donor values")
    }

    missing_cols <- setdiff(c("snp_id", "donor"), colnames(donor_snp_info))
    if (length(missing_cols) > 0) {
        problems <- c(
            problems,
            paste0(
                "donor_snp_info is missing required column(s): ",
                paste(missing_cols, collapse = ", ")
            )
        )
        return(problems)
    }

    # A (snp_id, donor) pair may carry one row per zygosity_source (e.g. a
    # "vireo_gt" row and a "binomial" row coexisting) -- so the key only
    # widens to include zygosity_source when that column is actually present.
    dup_key <- if ("zygosity_source" %in% colnames(donor_snp_info)) {
        c("snp_id", "donor", "zygosity_source")
    } else {
        c("snp_id", "donor")
    }
    if (any(duplicated(donor_snp_info[dup_key]))) {
        problems <- c(problems, sprintf("donor_snp_info contains duplicate (%s) rows", paste(dup_key, collapse = ", ")))
    }

    unknown_snps <- setdiff(donor_snp_info$snp_id, snp_info$snp_id)
    if (length(unknown_snps) > 0) {
        problems <- c(
            problems,
            sprintf(
                "donor_snp_info references %d snp_id(s) not present in snp_info",
                length(unknown_snps)
            )
        )
    }

    unknown_donors <- setdiff(donor_snp_info$donor, donor_info$donor)
    if (length(unknown_donors) > 0) {
        problems <- c(
            problems,
            sprintf(
                "donor_snp_info references %d donor(s) not present in donor_info",
                length(unknown_donors)
            )
        )
    }

    if (all(c("zygosity", "zygosity_source") %in% colnames(donor_snp_info))) {
        missing_source <- !is.na(donor_snp_info$zygosity) & is.na(donor_snp_info$zygosity_source)
        if (any(missing_source)) {
            problems <- c(
                problems,
                sprintf(
                    "donor_snp_info has %d row(s) with a zygosity call but no zygosity_source",
                    sum(missing_source)
                )
            )
        }
    }

    problems
}

.empty_donor_snp_info <- function() {
    tibble::tibble(
        snp_id = character(0),
        donor = character(0),
        zygosity = character(0),
        zygosity_source = character(0),
        zygosity_p_val = double(0),
        zygosity_adj_p_val = double(0),
        zygosity_gt_prob = double(0),
        xci_informative = logical(0),
        allele_on_x1 = character(0),
        xci_escape_fraction = double(0)
    )
}

.apply_donor_map <- function(donor, donor_map) {
    if (is.null(donor_map)) {
        return(donor)
    }
    if (
        !is.character(donor_map) ||
            is.null(names(donor_map)) ||
            any(names(donor_map) == "") ||
            any(duplicated(names(donor_map)))
    ) {
        stop(
            "donor_map must be a named character vector (new donor label = old donor label) with unique, non-empty names"
        )
    }
    if (any(duplicated(donor_map))) {
        stop("donor_map must not map more than one new label to the same old donor")
    }

    idx <- match(donor, donor_map)
    matched <- !is.na(idx)
    # Not `ifelse()`: on zero-length input it returns `logical(0)` regardless
    # of `donor`'s type (a base R quirk), which breaks a downstream join
    # expecting `donor` to stay character.
    result <- donor
    result[matched] <- names(donor_map)[idx[matched]]
    result
}

.derive_zygosity_source <- function(donor_snp_info) {
    if (!"zygosity_source" %in% colnames(donor_snp_info)) {
        return(NA_character_)
    }
    sources <- unique(stats::na.omit(donor_snp_info$zygosity_source))
    if (length(sources) == 1) sources else NA_character_
}

.default_donor_info <- function(barcode_info) {
    if ("donor" %in% colnames(barcode_info)) {
        tibble::tibble(donor = sort(unique(stats::na.omit(barcode_info$donor))))
    } else {
        tibble::tibble(donor = character(0))
    }
}

# One row per library present in barcode_info. Derived rather than accepted
# from the caller: which libraries an object holds is a property of its cells,
# so there is nothing for a constructor argument to say that barcode_info does
# not already.
.default_library_info <- function(barcode_info) {
    if (!"library_id" %in% colnames(barcode_info)) {
        return(.empty_library_info())
    }
    libraries <- sort(unique(stats::na.omit(barcode_info$library_id)))
    if (length(libraries) == 0) {
        return(.empty_library_info())
    }
    tibble::tibble(
        library_id = libraries,
        n_cells = as.integer(tabulate(match(barcode_info$library_id, libraries), nbins = length(libraries)))
    )
}

.empty_library_info <- function() {
    tibble::tibble(
        library_id = character(0),
        n_cells = integer(0)
    )
}

# Re-derives library_info after barcode_info has been written to directly.
# Editing barcode_info$library_id is the one way to change which libraries an
# object holds without going through the constructor, so without this the two
# tables drift apart.
.resync_library_info <- function(x, barcode_info) {
    if (!methods::.hasSlot(x, "library_info")) {
        return(x)
    }
    x@library_info <- .default_library_info(barcode_info)
    x
}

# Carries the molecule calls from one object onto another rebuilt from it.
# Passed through whole rather than cut down to the surviving cells and SNPs:
# every function that counts them joins them to the object's own cells and
# phased SNPs first, so calls outside the object are never counted, and
# SNPData need not know how the calls are structured.
.propagate_molecules <- function(object, from) {
    object@molecules <- molecules(from)
    object
}

# The gene annotation columns an SNPData object keeps. strand is needed by the
# molecule functions, which build their SNP-to-gene map from this copy.
.GENE_ANNO_COLUMNS <- c("chrom", "start", "end", "gene_name", "strand")

.empty_gene_anno <- function() {
    tibble::tibble(
        chrom = character(0),
        start = integer(0),
        end = integer(0),
        gene_name = character(0),
        strand = character(0)
    )
}

# Trims a gene annotation to the columns stored on SNPData, dropping the rest
# (a GFF's attribute strings would otherwise bloat every saved object) and
# any attributes such as readr's column spec, so two copies of one annotation
# compare identical.
.as_stored_gene_anno <- function(gene_anno) {
    if (is.null(gene_anno)) {
        return(.empty_gene_anno())
    }
    missing_cols <- setdiff(.GENE_ANNO_COLUMNS, colnames(gene_anno))
    if (length(missing_cols) > 0) {
        stop("gene_anno is missing required column(s): ", paste(missing_cols, collapse = ", "))
    }
    .check_gene_strand(gene_anno$strand)
    stored <- tibble::as_tibble(lapply(as.list(gene_anno)[.GENE_ANNO_COLUMNS], as.vector))
    stored$chrom <- as.character(stored$chrom)
    dplyr::distinct(stored)
}

# Rejects gene strands other than "+"/"-". A molecule is credited to one of
# two overlapping genes by matching its own transcript strand to the gene's,
# so a gene with "." or "*" (common in GFFs) or NA would silently never
# receive a molecule wherever it overlaps another gene.
.check_gene_strand <- function(strand, arg_name = "gene_anno") {
    invalid <- unique(as.character(strand)[!as.character(strand) %in% c("+", "-")])
    if (length(invalid) > 0) {
        stop(
            arg_name,
            "$strand must be \"+\" or \"-\" for every gene; found: ",
            paste(ifelse(is.na(invalid), "NA", paste0("\"", invalid, "\"")), collapse = ", "),
            ". Drop or assign a strand to unstranded genes."
        )
    }
    invisible(NULL)
}

# Carries the gene annotation onto an object rebuilt from another. Kept whole
# rather than cut to the surviving SNPs: it describes the genome, not the
# object's contents.
.propagate_gene_anno <- function(object, from) {
    object@gene_anno <- gene_anno(from)
    object
}

# Re-keys stored molecule calls when barcode_info is about to be replaced.
# Calls are matched to cells on (library_id, barcode), so editing either column
# after phase_from_molecules() would otherwise leave them matching nothing, or
# the wrong cell. Each cell's old key is mapped to its new one; calls for cells
# no longer in the object have no new key and are dropped, since they could
# only ever be counted against a cell that later took over their old key.
.rekey_molecule_calls <- function(x, new_barcode_info) {
    molecules <- molecules(x)
    if (!.has_molecule_calls(molecules)) {
        return(x)
    }
    old_keys <- .cell_key_columns(x@barcode_info)
    new_keys <- .cell_key_columns(new_barcode_info)[match(x@barcode_info$cell_id, new_barcode_info$cell_id), ]
    if (identical(old_keys, new_keys)) {
        return(x)
    }
    if (anyDuplicated(new_keys) > 0) {
        stop(
            "barcode_info would give two cells the same (library_id, barcode) pair, so the molecule calls ",
            "phase_from_molecules() stored could not be told apart. Keep each pair unique."
        )
    }

    key_map <- dplyr::bind_cols(old_keys, dplyr::rename_with(new_keys, ~ paste0("new_", .x)))
    molecules@calls <- molecules@calls %>%
        dplyr::inner_join(key_map, by = c("library_id", "barcode")) %>%
        dplyr::mutate(library_id = new_library_id, barcode = new_barcode) %>%
        dplyr::select(-new_library_id, -new_barcode)
    x@molecules <- molecules
    x
}

# The (library_id, barcode) key molecule calls are matched on, as character
# columns, with NA standing in for a column barcode_info does not carry.
.cell_key_columns <- function(barcode_info) {
    key_column <- function(name) {
        if (!name %in% colnames(barcode_info)) {
            return(rep(NA_character_, nrow(barcode_info)))
        }
        as.character(barcode_info[[name]])
    }
    tibble::tibble(library_id = key_column("library_id"), barcode = key_column("barcode"))
}

# Moves the molecule data older versions kept elsewhere into the molecules
# slot: the "molecule_calls" and "bam_calibration" attributes
# phase_from_molecules() used to attach, the retired snp_gene_map slot, and the
# BAM paths library_info used to record. Old calls were keyed on donor, so each
# is given the library of the cell it names.
.migrate_molecule_attributes <- function(object) {
    old_calls <- attr(object, "molecule_calls")
    calls <- NULL
    if (!is.null(old_calls)) {
        barcode_info <- object@barcode_info
        if (!"library_id" %in% colnames(barcode_info)) {
            barcode_info$library_id <- NA_character_
        }
        cell_library <- dplyr::distinct(barcode_info, donor, barcode, library_id)
        calls <- old_calls %>%
            dplyr::left_join(cell_library, by = c("donor", "barcode")) %>%
            dplyr::select(dplyr::all_of(.MOLECULE_CALL_COLUMNS))
    }

    bam_files <- list()
    if (methods::.hasSlot(object, "library_info") && "bam_files" %in% colnames(object@library_info)) {
        stored <- object@library_info[lengths(object@library_info$bam_files) > 0, , drop = FALSE]
        bam_files <- stats::setNames(stored$bam_files, stored$library_id)
        object@library_info$bam_files <- NULL
    }

    molecules <- MoleculeCalls(
        calls = calls,
        snp_gene_map = attr(object, "snp_gene_map"),
        bam_calibration = attr(object, "bam_calibration"),
        bam_files = bam_files
    )
    attr(object, "molecule_calls") <- NULL
    attr(object, "bam_calibration") <- NULL
    attr(object, "snp_gene_map") <- NULL
    object@molecules <- molecules
    object
}

# Normalises the `bam_files` argument shared by phase_from_molecules() and
# molecule_haplotype_counts() into a named list of character vectors, one
# element per library.
.as_library_bam_list <- function(bam_files, arg_name = "bam_files") {
    if (length(bam_files) == 0) {
        stop(arg_name, " is empty; supply at least one library's BAM file(s).")
    }
    nms <- names(bam_files)
    if (is.null(nms) || anyNA(nms) || any(!nzchar(nms))) {
        stop(arg_name, " must be named, library_id = path(s).")
    }
    if (anyDuplicated(nms) > 0) {
        stop(
            arg_name,
            " has repeated library_id name(s): ",
            paste(unique(nms[duplicated(nms)]), collapse = ", "),
            ". Give each library one entry listing all of its BAM files."
        )
    }
    bam_files <- as.list(bam_files)
    for (lib in nms) {
        paths <- bam_files[[lib]]
        if (!is.character(paths) || length(paths) == 0 || anyNA(paths)) {
            stop(arg_name, "[['", lib, "']] must be a non-empty character vector of BAM paths.")
        }
    }
    bam_files
}

.assign_snp_ids <- function(snp_info) {
    if (!"snp_id" %in% colnames(snp_info)) {
        if (all(c("chrom", "pos", "ref", "alt") %in% colnames(snp_info))) {
            snp_info$snp_id <- make_snp_id(
                snp_info$chrom,
                snp_info$pos,
                snp_info$ref,
                snp_info$alt
            )
        } else {
            snp_info$snp_id <- paste0("snp_", seq_len(nrow(snp_info)))
        }
    }

    snp_info
}

.assign_cell_ids <- function(barcode_info) {
    if (!"cell_id" %in% colnames(barcode_info)) {
        barcode_info$cell_id <- paste0("cell_", seq_len(nrow(barcode_info)))
    }

    barcode_info
}

# Guarantees a `library_id` column on barcode_info, filled with NA where the
# caller did not supply one.
#
# The column records which sequencing library each cell's barcode was drawn
# from, which is what lets merge_snpdata() tell a barcode shared between two
# libraries (two different cells that collided on a 737k-barcode whitelist)
# from a barcode shared within one library (the same physical cell sequenced
# twice). Always present rather than optional so that merging code can key on
# (library_id, barcode) without first testing whether either object has the
# column; an all-NA column simply means the libraries are unlabelled.
.assign_library_id <- function(barcode_info) {
    if (!"library_id" %in% colnames(barcode_info)) {
        barcode_info$library_id <- NA_character_
    }

    barcode_info$library_id <- as.character(barcode_info$library_id)

    anchor <- if ("barcode" %in% colnames(barcode_info)) "barcode" else "cell_id"
    dplyr::relocate(barcode_info, "library_id", .after = dplyr::all_of(anchor))
}

.dedupe_snps <- function(ref_count, alt_count, oth_count, snp_info, donor_snp_info, total_count = NULL) {
    if (!any(duplicated(snp_info$snp_id))) {
        return(list(
            ref_count = ref_count,
            alt_count = alt_count,
            oth_count = oth_count,
            snp_info = snp_info,
            donor_snp_info = donor_snp_info,
            total_count = total_count
        ))
    }

    dup_positions <- which(duplicated(snp_info$snp_id))
    dup_labels <- paste0(
        snp_info$snp_id[dup_positions],
        " (row ",
        dup_positions,
        ")"
    )
    dup_snps_msg <- paste(head(dup_labels, 5), collapse = ", ")
    if (length(dup_positions) > 5) {
        dup_snps_msg <- paste0(
            dup_snps_msg,
            ", ... (",
            length(dup_positions),
            " duplicates)"
        )
    }

    warning(
        sprintf(
            "Duplicate SNP IDs detected (%s). Keeping first occurrence and dropping duplicates.",
            dup_snps_msg
        ),
        call. = FALSE
    )

    keep_snps <- !duplicated(snp_info$snp_id)
    kept_snp_info <- snp_info[keep_snps, , drop = FALSE]
    list(
        ref_count = ref_count[keep_snps, , drop = FALSE],
        alt_count = alt_count[keep_snps, , drop = FALSE],
        oth_count = oth_count[keep_snps, , drop = FALSE],
        snp_info = kept_snp_info,
        donor_snp_info = donor_snp_info[donor_snp_info$snp_id %in% kept_snp_info$snp_id, , drop = FALSE],
        total_count = if (is.null(total_count)) NULL else total_count[keep_snps, , drop = FALSE]
    )
}

.set_dimnames <- function(ref_count, alt_count, oth_count, snp_info, barcode_info) {
    colnames(ref_count) <- barcode_info$cell_id
    colnames(alt_count) <- barcode_info$cell_id
    colnames(oth_count) <- barcode_info$cell_id
    rownames(ref_count) <- snp_info$snp_id
    rownames(alt_count) <- snp_info$snp_id
    rownames(oth_count) <- snp_info$snp_id

    list(
        ref_count = ref_count,
        alt_count = alt_count,
        oth_count = oth_count
    )
}

# Converts a `[` index (numeric, negative, logical, or character) into
# positions along one dimension, erroring on any index that does not name or
# point at an existing SNP/cell rather than letting it become an NA row.
.as_index_positions <- function(index, dim_names, what) {
    positions <- stats::setNames(seq_along(dim_names), dim_names)[index]
    if (anyNA(positions)) {
        bad <- if (is.logical(index)) which(is.na(positions)) else index[is.na(positions)]
        stop(sprintf(
            "Cannot subset SNPData: %s index not found: %s",
            what,
            paste(utils::head(bad, 5), collapse = ", ")
        ))
    }
    unname(positions)
}

.recompute_snp_stats <- function(snp_info, total_count) {
    snp_info$coverage <- Matrix::rowSums(total_count)
    snp_info$non_zero_samples <- Matrix::rowSums(total_count > 0)
    snp_info
}

.recompute_barcode_stats <- function(barcode_info, total_count) {
    barcode_info$library_size <- Matrix::colSums(total_count)
    barcode_info$non_zero_snps <- Matrix::colSums(total_count > 0)
    barcode_info
}

.recompute_donor_stats <- function(donor_info, barcode_info) {
    if (nrow(donor_info) == 0) {
        donor_info$n_cells <- integer(0)
        return(donor_info)
    }

    if ("donor" %in% colnames(barcode_info)) {
        cell_counts <- table(barcode_info$donor)
        donor_info$n_cells <- as.integer(cell_counts[donor_info$donor])
        donor_info$n_cells[is.na(donor_info$n_cells)] <- 0L
    } else {
        donor_info$n_cells <- 0L
    }

    donor_info
}

.recompute_metrics <- function(snp_info, barcode_info, donor_info, ref_count, alt_count, total_count = NULL) {
    # total_count = alt_count + ref_count by construction; a caller that
    # already has it (e.g. import_cellsnp()'s DP matrix) passes it through to
    # skip re-deriving it from two large sparse matrices.
    if (is.null(total_count)) {
        total_count <- alt_count + ref_count
    }
    list(
        snp_info = .recompute_snp_stats(snp_info, total_count),
        barcode_info = .recompute_barcode_stats(barcode_info, total_count),
        donor_info = .recompute_donor_stats(donor_info, barcode_info)
    )
}

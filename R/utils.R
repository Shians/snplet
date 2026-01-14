# Suppress R CMD check notes about global variables used in data.frame operations
utils::globalVariables(c(
    ".snp_row",
    "adj_p_val",
    "alt",
    "barcode",
    "binom_p_value",
    "cell",
    "cell_id",
    "chrom",
    "chrom_canonical",
    "clonotype",
    "coverage",
    "donor",
    "donor_id",
    "end",
    "gene_name",
    "library_size",
    "maf",
    "minor_count",
    "minor_allele_count",
    "n",
    "n_het_donors",
    "p_val",
    "pos",
    "ref",
    "seqnames",
    "signif_snps_clonotype",
    "snp_id",
    "start",
    "tested",
    "total_count",
    "total_library_size",
    "zygosity"
))

percentile_summary <- function(x, percentiles = c(0.1, 0.25, 0.75, 0.9, 0.95, 0.99)) {
    out <- c(min = min(x))
    past_median <- FALSE
    for (p in percentiles) {
        if (p > 0.5 && !past_median) {
            out <- c(out, median(x))
            names(out)[length(out)] <- "median"
            past_median <- TRUE
        }
        out <- c(out, as.numeric(quantile(x, p)))
        names(out)[length(out)] <- paste0("p", as.character(round(p * 100)))
    }

    out <- c(out, max = max(x))
    names(out)[length(out)] <- "max"

    out
}

groupedRowMeans <- function(x, groups) {
    # calculate the mean of each group
    out <- matrix(NA, ncol = length(unique(groups)), nrow = nrow(x))
    colnames(out) <- sort(unique(groups))
    rownames(out) <- rownames(x)

    for (i in unique(groups)) {
        out[, i] <- Matrix::rowMeans(
            x[, groups == i, drop = FALSE],
            na.rm = TRUE
        )
    }

    out
}

groupedRowSums <- function(x, groups) {
    # calculate the sum of each group
    out <- matrix(NA, ncol = length(unique(groups)), nrow = nrow(x))
    colnames(out) <- sort(unique(groups))
    rownames(out) <- rownames(x)

    for (i in unique(groups)) {
        out[, i] <- Matrix::rowSums(
            x[, groups == i, drop = FALSE],
            na.rm = TRUE
        )
    }

    out
}

#' Generate standardized SNP IDs from genomic coordinates
#'
#' Creates SNP identifiers in the format "chr:pos:ref:alt" for consistent
#' identification across datasets. This format is deterministic and contains
#' all information needed for variant identification.
#'
#' @param chrom Character vector of chromosome names
#' @param pos Integer vector of genomic positions
#' @param ref Character vector of reference alleles
#' @param alt Character vector of alternate alleles
#'
#' @return Character vector of SNP IDs in format "chr:pos:ref:alt"
#'
#' @examples
#' \dontrun{
#' snp_ids <- make_snp_id(
#'     chrom = c("chr1", "chr2"),
#'     pos = c(12345, 67890),
#'     ref = c("A", "G"),
#'     alt = c("G", "T")
#' )
#' # Returns: c("chr1:12345:A:G", "chr2:67890:G:T")
#' }
#'
#' @keywords internal
make_snp_id <- function(chrom, pos, ref, alt) {
    paste(chrom, pos, ref, alt, sep = ":")
}

#' Check if a file exists
#' @param path Path to the file to check
check_file <- function(path) {
    if (file.exists(path)) {
        logger::log_info("File found: {path}")
    } else {
        stop(glue::glue("Required file not found: {path}"))
    }
}

# ==============================================================================
# Chromosome Name Utilities
# ==============================================================================

#' Load chromosome name lookup table
#'
#' @return A data.frame with chromosome name mappings
#' @keywords internal
.load_chr_table <- function() {
    table_path <- system.file("extdata", "chr_name_table.tsv", package = "snplet")
    if (!file.exists(table_path)) {
        stop("Chromosome name lookup table not found in package installation")
    }
    readr::read_tsv(
        table_path,
        col_types = readr::cols(.default = readr::col_character())
    )
}

#' Detect chromosome naming style
#'
#' Automatically detects the chromosome naming convention used in a character vector.
#' Returns one of: "numeric", "ucsc", "refseq_mouse", "genbank_mouse", "refseq_human",
#' "genbank_human", or "unknown".
#'
#' @param chr_names Character vector of chromosome names
#' @return Character string indicating the detected style
#'
#' @examples
#' \dontrun{
#' detect_chr_style(c("1", "2", "X"))  # Returns "numeric"
#' detect_chr_style(c("chr1", "chr2", "chrX"))  # Returns "ucsc"
#' detect_chr_style(c("NC_000067.6", "NC_000068.7"))  # Returns "refseq_mouse"
#' }
#'
#' @export
detect_chr_style <- function(chr_names) {
    if (length(chr_names) == 0) {
        return("unknown")
    }

    chr_table <- .load_chr_table()

    # Remove NA values for matching
    chr_names_clean <- chr_names[!is.na(chr_names)]
    if (length(chr_names_clean) == 0) {
        return("unknown")
    }

    # Get unique chromosome names to avoid bias from duplicates
    unique_chrs <- unique(chr_names_clean)

    # Check each style column for matches
    style_cols <- c("numeric", "ucsc", "refseq_mouse", "genbank_mouse", "refseq_human", "genbank_human")

    detected_styles <- character(0)

    for (style in style_cols) {
        table_values <- chr_table[[style]]
        table_values <- table_values[!is.na(table_values)]

        # Count how many unique input chromosomes match this style
        matching_chrs <- unique_chrs[unique_chrs %in% table_values]
        n_matches <- length(matching_chrs)

        if (n_matches == 0) {
            next
        }

        # Calculate fraction of input chromosomes that match
        input_fraction <- n_matches / length(unique_chrs)

        # A style is valid if it matches any of the input chromosomes
        is_valid <- n_matches > 0

        if (is_valid) {
            detected_styles <- c(detected_styles, style)
        }
    }

    # Since chromosome naming styles are mutually exclusive,
    # we should never detect more than one style
    if (length(detected_styles) > 1) {
        stop(
            "Multiple chromosome styles detected: ",
            paste(detected_styles, collapse = ", "),
            ". This should not happen as naming conventions are mutually exclusive. ",
            "Please report this as a bug."
        )
    }

    if (length(detected_styles) == 1) {
        return(detected_styles[1])
    }

    "unknown"
}

#' Normalize chromosome names to canonical form
#'
#' Converts chromosome names from any recognized style to the canonical UCSC form (chr1, chr2, chrX).
#' Recognized styles include: numeric, UCSC (chr-prefix), RefSeq, and GenBank accessions
#' for mouse and human genomes.
#'
#' @param chr_names Character vector of chromosome names
#' @param from_style Style of input chromosome names. If "auto" (default), the style
#'   will be auto-detected. Can be one of: "auto", "numeric", "ucsc", "refseq_mouse",
#'   "genbank_mouse", "refseq_human", "genbank_human", or "unknown"
#'
#' @return Character vector of chromosome names in canonical UCSC form
#'
#' @examples
#' \dontrun{
#' normalize_chr_names(c("1", "2", "X"))  # Returns c("chr1", "chr2", "chrX")
#' normalize_chr_names(c("NC_000067.6", "NC_000086.7"), from_style = "refseq_mouse")
#' # Returns c("chr1", "chrX")
#' }
#'
#' @export
normalize_chr_names <- function(chr_names, from_style = "auto") {
    if (from_style == "auto") {
        from_style <- detect_chr_style(chr_names)
    }

    if (from_style == "unknown") {
        warning("Chromosome style is unknown. Returning original names.")
        return(chr_names)
    }

    if (from_style == "ucsc") {
        return(chr_names)
    }

    chr_table <- .load_chr_table()

    # Create lookup from from_style -> ucsc (canonical)
    lookup <- stats::setNames(chr_table$ucsc, chr_table[[from_style]])

    # Map chromosomes, keeping NA for unmatched values
    normalized <- lookup[chr_names]

    # For any NA values, keep original
    normalized[is.na(normalized)] <- chr_names[is.na(normalized)]

    as.character(normalized)
}

#' Convert chromosome names between styles
#'
#' Converts chromosome names from one naming convention to another. Supports conversions
#' between numeric, UCSC, RefSeq, and GenBank naming styles for mouse and human genomes.
#'
#' @param chr_names Character vector of chromosome names
#' @param from_style Style of input chromosome names. If "auto" (default), the style
#'   will be auto-detected
#' @param to_style Target style for chromosome names. One of: "numeric", "ucsc",
#'   "refseq_mouse", "genbank_mouse", "refseq_human", "genbank_human"
#'
#' @return Character vector of chromosome names in the target style
#'
#' @examples
#' \dontrun{
#' convert_chr_style(c("1", "2", "X"), to_style = "ucsc")
#' # Returns c("chr1", "chr2", "chrX")
#'
#' convert_chr_style(c("chr1", "chrX"), from_style = "ucsc", to_style = "refseq_mouse")
#' # Returns c("NC_000067.6", "NC_000086.7")
#' }
#'
#' @export
convert_chr_style <- function(chr_names, from_style = "auto", to_style = "ucsc") {
    # First normalize to canonical UCSC form
    canonical <- normalize_chr_names(chr_names, from_style)

    if (to_style == "ucsc") {
        return(canonical)
    }

    chr_table <- .load_chr_table()

    # Create lookup from ucsc -> to_style
    lookup <- stats::setNames(chr_table[[to_style]], chr_table$ucsc)

    # Map chromosomes
    converted <- lookup[canonical]

    # For any NA values (not in table), keep canonical form
    converted[is.na(converted)] <- canonical[is.na(converted)]

    as.character(converted)
}

#' Validate that chromosome style is known
#'
#' Helper function to check if a SNPData object has a known chromosome style.
#' This is required for functions that need to perform chromosome-specific operations.
#'
#' @param x A SNPData object
#' @param operation_name Name of the operation requiring known chr_style (for error message)
#'
#' @return Invisibly returns TRUE if chr_style is known
#' @keywords internal
.validate_chr_style <- function(x, operation_name = "This operation") {
    if (!methods::is(x, "SNPData")) {
        stop("Expected a SNPData object")
    }

    # Use accessor method for backwards compatibility
    style <- chr_style(x)

    if (style == "unknown") {
        stop(
            sprintf(
                "%s requires a known chromosome naming style, but chr_style is 'unknown'. ",
                operation_name
            ),
            "This typically happens with novel genomes or non-standard chromosome names. ",
            "You may need to manually standardize chromosome names before using this function."
        )
    }

    invisible(TRUE)
}

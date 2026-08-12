# Suppress R CMD check notes about global variables used in data.frame operations
utils::globalVariables(c(
    ".data",
    ".env",
    ".snp_row",
    "adj_p_val",
    "allele_on_x1",
    "alt",
    "barcode",
    "cell",
    "cell_id",
    "gene",
    "chrom",
    "chrom_canonical",
    "clonotype",
    "coverage",
    "donor",
    "donor_id",
    "end",
    "gene_name",
    "active_x",
    "pi_g",
    "escape_fraction",
    "xci_informative",
    "xci_informative_donor",
    "xci_allele_on_x1_by_donor",
    "xci_escape_fraction_by_donor",
    "xci_escape_fraction",
    "xci_post_X1_active",
    "xci_fit_unit",
    "unit_id",
    "library_id",
    "library_size",
    "maf",
    "minor_allele_count",
    "n",
    "p_val",
    "pos",
    "ref",
    "assignment",
    "ll_x1_inactive",
    "ll_x2_inactive",
    "llr",
    "lor",
    "pi",
    "post_X1_active",
    "post_X2_active",
    "xi_ref_count",
    "xi_total",
    "seqnames",
    "snp_id",
    "start",
    "tested",
    "total_count",
    "active_count",
    "inactive_count",
    "dominant_allele",
    "same_allele_dominant",
    "phase_contradiction",
    "phase_likely_inverted",
    "snp_coverage",
    "escapes",
    "zygosity",
    "zygosity_source",
    "umi",
    "allele",
    "n_calls",
    "n_files",
    "qname",
    "base",
    "base_quality",
    "read_idx",
    "snp_idx",
    "snp_a",
    "snp_b",
    "same",
    "n_same",
    "relation",
    "consistency",
    "n_genes",
    "allele_on_h1",
    "allele_on_x1_em",
    "allele_on_x1_molecule",
    "implies_h1_is_x1",
    "block",
    "n_anchors",
    "n_x1",
    "n_x2",
    "orientation",
    "own_vote",
    "is_outlier_anchor",
    "phase_conflict",
    "phase_source",
    "phase_block",
    "is_x1",
    "haplotype",
    "molecules",
    "dominant_molecules",
    "n_stranded_molecules",
    "ambiguous",
    "gene_strand",
    "group_key",
    "n_genes_same_strand",
    "n_reads",
    "snps",
    "strand",
    "transcript_strand",
    "xci_median_pi_g",
    "xci_rho",
    "n_cells"
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

groupedRowSums <- function(x, groups) {
    if (length(groups) != ncol(x)) {
        stop(
            "Length of groups must match the number of columns in x. ",
            "Got ",
            length(groups),
            " groups for ",
            ncol(x),
            " columns."
        )
    }

    if (anyNA(groups)) {
        stop("groups must not contain NA values.")
    }

    groups <- as.character(groups)
    group_levels <- sort(unique(groups))
    out <- matrix(
        NA_real_,
        nrow = nrow(x),
        ncol = length(group_levels),
        dimnames = list(rownames(x), group_levels)
    )

    for (group_name in group_levels) {
        out[, group_name] <- Matrix::rowSums(x[, groups == group_name, drop = FALSE], na.rm = TRUE)
    }

    out
}

#' Generate standardised SNP IDs from genomic coordinates
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

#' Vectorised exact binomial test using the beta distribution
#'
#' Computes one-sided exact binomial test p-values matching
#' `stats::binom.test(x, n, p, alternative)$p.value`, but vectorised over
#' `x` and `n` and substantially faster.
#'
#' The equivalence relies on the identity between the binomial CDF and the
#' regularised incomplete beta function:
#'
#' \deqn{P(X \ge x) = I_p(x, n - x + 1) = \mathrm{pbeta}(p, x, n - x + 1)}
#' \deqn{P(X \le x) = 1 - I_p(x + 1, n - x) = \mathrm{pbeta}(p, x + 1, n - x, lower = FALSE)}
#'
#' where \eqn{X \sim \mathrm{Binomial}(n, p)}. So:
#' \itemize{
#'   \item `alternative = "greater"`:
#'     `binom.test` p-value is `P(X >= x)` =
#'     `pbeta(p, x, n - x + 1)`. Edge case: when `x == 0`, the p-value is 1.
#'   \item `alternative = "less"`:
#'     `binom.test` p-value is `P(X <= x)` =
#'     `pbeta(p, x + 1, n - x, lower.tail = FALSE)`. Edge case: when `x == n`,
#'     the p-value is 1.
#' }
#'
#' Two-sided tests are not supported because `binom.test`'s two-sided p-value
#' sums tail probabilities `<=` the observed point probability, which cannot
#' be expressed as a single `pbeta` call.
#'
#' @param x Integer vector of successes.
#' @param n Integer vector of trials. Recycled with `x`.
#' @param p Numeric scalar or vector giving the hypothesised probability of
#'   success under the null. Recycled with `x` and `n`.
#' @param alternative One of `"greater"` or `"less"`.
#'
#' @return Numeric vector of p-values, the same length as the recycled inputs.
#'
#' @keywords internal
binom_test <- function(x, n, p, alternative = c("greater", "less")) {
    alternative <- match.arg(alternative)

    if (any(x < 0, na.rm = TRUE) || any(x > n, na.rm = TRUE)) {
        stop("x must satisfy 0 <= x <= n")
    }
    if (any(p < 0 | p > 1, na.rm = TRUE)) {
        stop("p must be in [0, 1]")
    }

    if (alternative == "greater") {
        # P(X >= x) under Binomial(n, p); pbeta returns 0 at shape1 = 0,
        # but the true p-value at x = 0 is 1.
        p_val <- stats::pbeta(p, x, n - x + 1, lower.tail = TRUE)
        p_val[x == 0] <- 1
    } else {
        # P(X <= x); pbeta returns 0 at shape2 = 0, but the true p-value at
        # x = n is 1.
        p_val <- stats::pbeta(p, x + 1, n - x, lower.tail = FALSE)
        p_val[x == n] <- 1
    }

    p_val
}

#' Vectorised exact beta-binomial test
#'
#' Computes one-sided beta-binomial test p-values, the overdispersed analogue
#' of \code{\link{binom_test}}. Uses \code{VGAM::pbetabinom} directly rather
#' than a closed-form identity (unlike \code{binom_test}, no beta-CDF shortcut
#' is available for the beta-binomial), but is still fully vectorised over
#' `x`, `n`, `p` and `rho`.
#'
#' \itemize{
#'   \item `alternative = "greater"`: \eqn{P(X \ge x) = 1 - \mathrm{pbetabinom}(x - 1, n, p, \rho)}.
#'     Edge case: when `x == 0`, the p-value is 1.
#'   \item `alternative = "less"`: \eqn{P(X \le x) = \mathrm{pbetabinom}(x, n, p, \rho)}.
#'     Edge case: when `x == n`, the p-value is 1.
#' }
#'
#' `rho = 0` reduces to the exact binomial test (VGAM's beta-binomial CDF is
#' well-defined at `rho = 0`).
#'
#' @param x Integer vector of successes.
#' @param n Integer vector of trials. Recycled with `x`.
#' @param p Numeric scalar or vector giving the hypothesised probability of
#'   success under the null. Recycled with `x` and `n`.
#' @param rho Numeric scalar or vector in `[0, 1)` giving the beta-binomial
#'   overdispersion. Recycled with `x`, `n` and `p`.
#' @param alternative One of `"greater"` or `"less"`.
#'
#' @return Numeric vector of p-values, the same length as the recycled inputs.
#'
#' @keywords internal
betabinom_test <- function(x, n, p, rho, alternative = c("greater", "less")) {
    alternative <- match.arg(alternative)

    if (any(x < 0, na.rm = TRUE) || any(x > n, na.rm = TRUE)) {
        stop("x must satisfy 0 <= x <= n")
    }
    if (any(p < 0 | p > 1, na.rm = TRUE)) {
        stop("p must be in [0, 1]")
    }
    if (any(rho < 0 | rho >= 1, na.rm = TRUE)) {
        stop("rho must be in [0, 1)")
    }

    # The boundary quantile is pinned to 1 rather than evaluated: P(X >= 0) and
    # P(X <= n) are 1 by definition, and passing the out-of-support quantile
    # (x - 1 = -1, or x = n for "less") to VGAM makes it warn from its own
    # internal indexing. Recycle first so the mask lines up with the result.
    arg_lengths <- lengths(list(x, n, p, rho))
    # R's own recycling rule: any zero-length argument makes the result empty,
    # rather than the longest argument's length.
    len <- if (any(arg_lengths == 0)) 0L else max(arg_lengths)
    x <- rep_len(x, len)
    n <- rep_len(n, len)
    p <- rep_len(p, len)
    rho <- rep_len(rho, len)

    # An NA count is not on the boundary; it falls through to VGAM and comes
    # back NA, as it did before this branch existed.
    at_boundary <- if (alternative == "greater") x == 0 else x == n
    boundary <- !is.na(at_boundary) & at_boundary

    p_val <- rep(1, len)
    if (any(!boundary)) {
        keep <- !boundary
        quantile <- if (alternative == "greater") x[keep] - 1 else x[keep]
        cdf <- VGAM::pbetabinom(quantile, n[keep], p[keep], rho[keep])
        p_val[keep] <- if (alternative == "greater") 1 - cdf else cdf
    }

    # VGAM::pbetabinom can overshoot 1 by floating-point error (observed up to
    # ~4e-16) for some (x, n, p, rho) combinations, which turns 1 - pbetabinom(...)
    # slightly negative. Clamp to the documented [0, 1] range.
    pmin(pmax(p_val, 0), 1)
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
#' @family chromosome naming functions
#' @keywords internal
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

#' Normalise chromosome names to canonical form
#'
#' Converts chromosome names from any recognised style to the canonical UCSC form (chr1, chr2, chrX).
#' Recognised styles include: numeric, UCSC (chr-prefix), RefSeq, and GenBank accessions
#' for mouse and human genomes.
#'
#' @param chr_names Character vector of chromosome names
#' @param from_style Style of input chromosome names. If "auto" (default), the style
#'   will be auto-detected. Can be one of: "auto", "numeric", "ucsc", "refseq_mouse",
#'   "genbank_mouse", "refseq_human", "genbank_human", or "unknown"
#'
#' @return Character vector of chromosome names in canonical UCSC form
#'
#' @importFrom magrittr set_names
#' @noRd
normalise_chr_names <- function(chr_names, from_style = "auto") {
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
    lookup <- magrittr::set_names(chr_table$ucsc, chr_table[[from_style]])

    # Map chromosomes, keeping NA for unmatched values
    normalised <- lookup[chr_names]

    # For any NA values, keep original
    normalised[is.na(normalised)] <- chr_names[is.na(normalised)]

    as.character(normalised)
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
#' @family chromosome naming functions
#' @keywords internal
convert_chr_style <- function(chr_names, from_style = "auto", to_style = "ucsc") {
    # First normalise to canonical UCSC form
    canonical <- normalise_chr_names(chr_names, from_style)

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

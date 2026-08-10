# ==============================================================================
# Test Suite: SNPData S4 Class
# Description: Tests for SNPData class constructor, accessors, and core methods
# ==============================================================================

library(testthat)
library(Matrix)

# ------------------------------------------------------------------------------
# Test Data Setup
# ------------------------------------------------------------------------------

# Create test matrices
test_alt_count <- Matrix::Matrix(matrix(c(1, 2, 3, 4), nrow = 2, ncol = 2))
dimnames(test_alt_count) <- list(c("snp_1", "snp_2"), c("cell_1", "cell_2"))
test_ref_count <- Matrix::Matrix(matrix(c(5, 6, 7, 8), nrow = 2, ncol = 2))
dimnames(test_ref_count) <- list(c("snp_1", "snp_2"), c("cell_1", "cell_2"))

# Create test metadata
test_snp_info <- data.frame(
    snp_id = c("snp_1", "snp_2"),
    pos = c(100, 200),
    stringsAsFactors = FALSE
)

test_barcode_info <- data.frame(
    cell_id = c("cell_1", "cell_2"),
    donor = c("donor_1", "donor_1"),
    clonotype = c("clonotype_1", "clonotype_2"),
    stringsAsFactors = FALSE
)

# Create test SNPData object using the standard fixtures above
create_test_snpdata <- function() {
    SNPData(
        ref_count = test_ref_count,
        alt_count = test_alt_count,
        snp_info = test_snp_info,
        barcode_info = test_barcode_info
    )
}

# ==============================================================================

test_that("SNPData() creates valid S4 object with correct class", {
    snp_data <- create_test_snpdata()

    # Verify SNPData object is created successfully
    expect_s4_class(snp_data, "SNPData")
})

test_that("SNPData() sets correct dimensions", {
    snp_data <- create_test_snpdata()

    # Check that dim() method returns correct dimensions
    expect_equal(dim(snp_data), c(2, 2))
    # Verify nrow() method works correctly
    expect_equal(nrow(snp_data), 2)
    # Verify ncol() method works correctly
    expect_equal(ncol(snp_data), 2)
})

test_that("SNPData() sets row and column names from metadata", {
    snp_data <- create_test_snpdata()

    # Check that SNP IDs are used as row names
    expect_equal(rownames(snp_data), c("snp_1", "snp_2"))
    # Check that cell IDs are used as column names
    expect_equal(colnames(snp_data), c("cell_1", "cell_2"))
})

test_that("ref_count() returns reference count matrix with correct values and names", {
    snp_data <- create_test_snpdata()

    ref_count_matrix <- ref_count(snp_data)

    # Verify ref_count accessor returns correct matrix values
    expect_equal(as.matrix(ref_count_matrix), as.matrix(test_ref_count))
    # Check that ref_count matrix has correct row names
    expect_equal(rownames(ref_count_matrix), c("snp_1", "snp_2"))
    # Check that ref_count matrix has correct column names
    expect_equal(colnames(ref_count_matrix), c("cell_1", "cell_2"))
})

test_that("alt_count() returns alternate count matrix with correct values and names", {
    snp_data <- create_test_snpdata()

    alt_count_matrix <- alt_count(snp_data)

    # Verify alt_count accessor returns correct matrix values
    expect_equal(as.matrix(alt_count_matrix), as.matrix(test_alt_count))
    # Check that alt_count matrix has correct row names
    expect_equal(rownames(alt_count_matrix), c("snp_1", "snp_2"))
    # Check that alt_count matrix has correct column names
    expect_equal(colnames(alt_count_matrix), c("cell_1", "cell_2"))
})

test_that("snp_info() returns SNP metadata with computed coverage metrics", {
    snp_data <- create_test_snpdata()

    snp_info <- snp_info(snp_data)

    # Verify SNP IDs are preserved from input
    expect_equal(snp_info$snp_id, test_snp_info$snp_id)
    # Verify SNP positions are preserved from input
    expect_equal(snp_info$pos, test_snp_info$pos)
    # Check that coverage column is automatically added
    expect_true("coverage" %in% colnames(snp_info))
    # Check that non_zero_samples column is automatically added
    expect_true("non_zero_samples" %in% colnames(snp_info))
    # Verify coverage calculation: rowSums of alt_count + ref_count
    expect_equal(as.numeric(snp_info$coverage), c(16, 20))
    # Verify non_zero_samples count: all samples have counts
    expect_equal(as.numeric(snp_info$non_zero_samples), c(2, 2))
})

test_that("barcode_info() returns barcode metadata with computed library size metrics", {
    snp_data <- create_test_snpdata()

    barcode_info <- barcode_info(snp_data)

    # Verify cell IDs are preserved from input
    expect_equal(barcode_info$cell_id, test_barcode_info$cell_id)
    # Verify donor information is preserved from input
    expect_equal(barcode_info$donor, test_barcode_info$donor)
    # Verify clonotype information is preserved from input
    expect_equal(barcode_info$clonotype, test_barcode_info$clonotype)
    # Check that library_size column is automatically added
    expect_true("library_size" %in% colnames(barcode_info))
    # Check that non_zero_snps column is automatically added
    expect_true("non_zero_snps" %in% colnames(barcode_info))
    # Verify library_size calculation: colSums of alt_count + ref_count
    expect_equal(as.numeric(barcode_info$library_size), c(14, 22))
    # Verify non_zero_snps count: all SNPs have counts
    expect_equal(as.numeric(barcode_info$non_zero_snps), c(2, 2))
})

test_that("SNPData() adds an all-NA library_id column when none is supplied", {
    snp_data <- create_test_snpdata()

    barcode_info <- barcode_info(snp_data)

    # Check that library_id is always present, so merging code can key on it
    # without first testing whether the object carries the column
    expect_true("library_id" %in% colnames(barcode_info))
    # Verify unlabelled libraries are recorded as NA rather than a placeholder
    expect_equal(barcode_info$library_id, rep(NA_character_, 2))
})

test_that("SNPData() preserves a caller-supplied library_id column", {
    barcode_info_labelled <- test_barcode_info
    barcode_info_labelled$library_id <- c("lib_A", "lib_A")

    snp_data <- SNPData(
        ref_count = test_ref_count,
        alt_count = test_alt_count,
        snp_info = test_snp_info,
        barcode_info = barcode_info_labelled
    )

    # Verify an existing library_id is carried through untouched
    expect_equal(barcode_info(snp_data)$library_id, c("lib_A", "lib_A"))
})

test_that("SNPData() coerces a non-character library_id to character", {
    barcode_info_numeric <- test_barcode_info
    barcode_info_numeric$library_id <- c(1L, 1L)

    snp_data <- SNPData(
        ref_count = test_ref_count,
        alt_count = test_alt_count,
        snp_info = test_snp_info,
        barcode_info = barcode_info_numeric
    )

    # Ensure the column type is stable regardless of how the caller labelled
    # libraries, since it is compared against NA_character_ during merges
    expect_type(barcode_info(snp_data)$library_id, "character")
})

test_that("[() subsetting carries library_id through to the subset", {
    barcode_info_labelled <- test_barcode_info
    barcode_info_labelled$library_id <- c("lib_A", "lib_B")
    snp_data <- SNPData(
        ref_count = test_ref_count,
        alt_count = test_alt_count,
        snp_info = test_snp_info,
        barcode_info = barcode_info_labelled
    )

    subset_data <- snp_data[, 2]

    # Verify the retained cell keeps its own library label
    expect_equal(barcode_info(subset_data)$library_id, "lib_B")
})

test_that("[() subsetting by numeric index returns SNPData object with correct dimensions", {
    snp_data <- create_test_snpdata()

    subset_data <- snp_data[1, 1]

    # Verify subsetting by index returns SNPData object
    expect_s4_class(subset_data, "SNPData")
    # Check that single element subset has correct dimensions
    expect_equal(dim(subset_data), c(1, 1))
})

test_that("[() subsetting by name returns SNPData object with correct dimensions", {
    snp_data <- create_test_snpdata()

    subset_data <- snp_data["snp_1", "cell_1"]

    # Verify subsetting by name returns SNPData object
    expect_s4_class(subset_data, "SNPData")
    # Check that named subset has correct dimensions
    expect_equal(dim(subset_data), c(1, 1))
})

test_that("[() ignores drop parameter and always returns SNPData object", {
    snp_data <- create_test_snpdata()

    subset_data <- snp_data[1, 1, drop = TRUE]

    # Verify drop=TRUE is ignored and still returns SNPData object
    expect_s4_class(subset_data, "SNPData")
    # Check that drop=TRUE doesn't affect dimensions
    expect_equal(dim(subset_data), c(1, 1))
})

test_that("[() subsetting to single row returns SNPData object", {
    snp_data <- create_test_snpdata()

    subset_data <- snp_data[1, ]

    # Verify single row subset returns SNPData object
    expect_s4_class(subset_data, "SNPData")
    # Check that single row subset has correct dimensions
    expect_equal(dim(subset_data), c(1, 2))
})

test_that("[() subsetting to single column returns SNPData object", {
    snp_data <- create_test_snpdata()

    subset_data <- snp_data[, 1]

    # Verify single column subset returns SNPData object
    expect_s4_class(subset_data, "SNPData")
    # Check that single column subset has correct dimensions
    expect_equal(dim(subset_data), c(2, 1))
})

test_that("coverage() returns sum of alt_count and ref_count matrices", {
    snp_data <- create_test_snpdata()
    expected_coverage <- Matrix::Matrix(
        matrix(c(6, 8, 10, 12), nrow = 2, ncol = 2),
        dimnames = list(c("snp_1", "snp_2"), c("cell_1", "cell_2"))
    )

    result <- coverage(snp_data)

    # Verify coverage calculation: alt_count + ref_count
    expect_equal(result, expected_coverage)
})
test_that("SNPData() auto-generates snp_id column when missing from snp_info", {
    test_snp_info_no_id <- data.frame(pos = c(100, 200))
    test_barcode_info_no_id <- data.frame(
        donor = c("donor_1", "donor_1"),
        clonotype = c("clonotype_1", "clonotype_2")
    )

    snp_data <- SNPData(
        alt_count = test_alt_count,
        ref_count = test_ref_count,
        snp_info = test_snp_info_no_id,
        barcode_info = test_barcode_info_no_id
    )

    # Verify snp_id column was auto-generated when missing
    expect_true("snp_id" %in% colnames(snp_info(snp_data)))
    # Check that auto-generated SNP IDs follow expected pattern
    expect_equal(snp_info(snp_data)$snp_id, c("snp_1", "snp_2"))
})

test_that("SNPData() auto-generates cell_id column when missing from barcode_info", {
    test_snp_info_no_id <- data.frame(pos = c(100, 200))
    test_barcode_info_no_id <- data.frame(
        donor = c("donor_1", "donor_1"),
        clonotype = c("clonotype_1", "clonotype_2")
    )

    snp_data <- SNPData(
        alt_count = test_alt_count,
        ref_count = test_ref_count,
        snp_info = test_snp_info_no_id,
        barcode_info = test_barcode_info_no_id
    )

    # Verify cell_id column was auto-generated when missing
    expect_true("cell_id" %in% colnames(barcode_info(snp_data)))
    # Check that auto-generated cell IDs follow expected pattern
    expect_equal(barcode_info(snp_data)$cell_id, c("cell_1", "cell_2"))
})

test_that("SNPData() throws error when alt_count and ref_count have mismatched row dimensions", {
    wrong_dim_alt_count <- Matrix::Matrix(matrix(c(1, 2, 3), nrow = 3, ncol = 1))

    # Verify error when alt_count and ref_count have different row counts
    expect_error(
        SNPData(
            ref_count = test_ref_count,
            alt_count = wrong_dim_alt_count,
            snp_info = test_snp_info,
            barcode_info = test_barcode_info
        ),
        "nrow\\(alt_count\\) == nrow\\(ref_count\\)"
    )
})

test_that("SNPData() throws error when snp_info rows don't match matrix rows", {
    wrong_dim_snp_info <- data.frame(snp_id = c("snp_1", "snp_2", "snp_3"))

    # Verify error when snp_info rows don't match matrix rows
    expect_error(
        SNPData(
            ref_count = test_ref_count,
            alt_count = test_alt_count,
            snp_info = wrong_dim_snp_info,
            barcode_info = test_barcode_info
        ),
        "nrow\\(ref_count\\) == nrow\\(snp_info\\)"
    )
})

test_that("SNPData() throws error when barcode_info rows don't match matrix columns", {
    wrong_dim_barcode_info <- data.frame(cell_id = c("cell_1", "cell_2", "cell_3"))

    # Verify error when barcode_info rows don't match matrix columns
    expect_error(
        SNPData(
            alt_count = test_alt_count,
            ref_count = test_ref_count,
            snp_info = test_snp_info,
            barcode_info = wrong_dim_barcode_info
        ),
        "ncol\\(alt_count\\) == nrow\\(barcode_info\\)"
    )
})

test_that("get_sample_info() returns same result as barcode_info()", {
    snp_data <- create_test_snpdata()

    # Verify get_sample_info is an alias for barcode_info
    expect_equal(get_sample_info(snp_data), barcode_info(snp_data))
})

test_that("barcode_info<- replaces barcode_info with valid data", {
    snp_data <- create_test_snpdata()

    # donor is left unchanged (donor_1, matching test_barcode_info): once
    # donor_info carries data, barcode_info<- rejects changing 'donor' -- see
    # the dedicated tests below for that restriction.
    new_barcode_info <- data.frame(
        cell_id = c("cell_1", "cell_2"),
        donor = c("donor_1", "donor_1"),
        clonotype = c("clonotype_3", "clonotype_4"),
        stringsAsFactors = FALSE
    )

    barcode_info(snp_data) <- new_barcode_info

    # Verify barcode_info was updated
    expect_equal(barcode_info(snp_data)$donor, c("donor_1", "donor_1"))
    # Verify clonotype was updated
    expect_equal(barcode_info(snp_data)$clonotype, c("clonotype_3", "clonotype_4"))
})

test_that("barcode_info<- throws error when dimensions don't match", {
    snp_data <- create_test_snpdata()

    wrong_dim_barcode_info <- data.frame(
        cell_id = c("cell_1", "cell_2", "cell_3"),
        donor = c("donor_1", "donor_1", "donor_1")
    )

    # Verify error when number of rows doesn't match
    expect_error(
        barcode_info(snp_data) <- wrong_dim_barcode_info,
        "Number of rows in barcode_info must match number of columns in count matrices"
    )
})

test_that("barcode_info<- throws error when cell_id column is missing", {
    snp_data <- create_test_snpdata()

    no_cell_id <- data.frame(
        donor = c("donor_1", "donor_1"),
        clonotype = c("clonotype_1", "clonotype_2")
    )

    # Verify error when cell_id column is missing
    expect_error(
        barcode_info(snp_data) <- no_cell_id,
        "barcode_info must contain a 'cell_id' column"
    )
})

test_that("barcode_info<- throws error when cell_id doesn't match matrix column names", {
    snp_data <- create_test_snpdata()

    mismatched_cell_id <- data.frame(
        cell_id = c("cell_3", "cell_4"),
        donor = c("donor_1", "donor_1")
    )

    # Verify error when cell_id doesn't match column names
    expect_error(
        barcode_info(snp_data) <- mismatched_cell_id,
        "barcode_info\\$cell_id must match column names of count matrices"
    )
})

test_that("snp_info<- replaces snp_info with valid data", {
    snp_data <- create_test_snpdata()

    new_snp_info <- data.frame(
        snp_id = c("snp_1", "snp_2"),
        pos = c(150, 250),
        chr = c("chr1", "chr2"),
        stringsAsFactors = FALSE
    )

    snp_info(snp_data) <- new_snp_info

    # Verify snp_info was updated
    expect_equal(snp_info(snp_data)$pos, c(150, 250))
    # Verify chr column was added
    expect_equal(snp_info(snp_data)$chr, c("chr1", "chr2"))
})

test_that("snp_info<- throws error when dimensions don't match", {
    snp_data <- create_test_snpdata()

    wrong_dim_snp_info <- data.frame(
        snp_id = c("snp_1", "snp_2", "snp_3"),
        pos = c(100, 200, 300)
    )

    # Verify error when number of rows doesn't match
    expect_error(
        snp_info(snp_data) <- wrong_dim_snp_info,
        "Number of rows in snp_info must match number of rows in count matrices"
    )
})

test_that("snp_info<- throws error when snp_id column is missing", {
    snp_data <- create_test_snpdata()

    no_snp_id <- data.frame(
        pos = c(100, 200),
        chr = c("chr1", "chr2")
    )

    # Verify error when snp_id column is missing
    expect_error(
        snp_info(snp_data) <- no_snp_id,
        "snp_info must contain a 'snp_id' column"
    )
})

test_that("snp_info<- throws error when snp_id doesn't match matrix row names", {
    snp_data <- create_test_snpdata()

    mismatched_snp_id <- data.frame(
        snp_id = c("snp_3", "snp_4"),
        pos = c(100, 200)
    )

    # Verify error when snp_id doesn't match row names
    expect_error(
        snp_info(snp_data) <- mismatched_snp_id,
        "snp_info\\$snp_id must match row names of count matrices"
    )
})

test_that("donor_info() auto-derives one row per distinct donor in barcode_info", {
    snp_data <- create_test_snpdata()

    # Verify the single donor in test_barcode_info produces one donor_info row
    expect_equal(donor_info(snp_data)$donor, "donor_1")
})

test_that("donor_snp_info() defaults to an empty tibble with the expected columns", {
    snp_data <- create_test_snpdata()

    donor_snp_info <- donor_snp_info(snp_data)
    # Verify no rows are present before any donor-level data has been added
    expect_equal(nrow(donor_snp_info), 0)
    # Verify the schema is already in place for later writers to join against
    expect_true(all(c("snp_id", "donor", "zygosity", "zygosity_source") %in% colnames(donor_snp_info)))
})

test_that("SNPData() throws error when donor_snp_info has duplicate (snp_id, donor) rows", {
    duplicate_donor_snp_info <- data.frame(
        snp_id = c("snp_1", "snp_1"),
        donor = c("donor_1", "donor_1"),
        zygosity = c("het", "hom")
    )

    # Verify validity rejects a duplicate key in donor_snp_info
    expect_error(
        SNPData(
            ref_count = test_ref_count,
            alt_count = test_alt_count,
            snp_info = test_snp_info,
            barcode_info = test_barcode_info,
            donor_snp_info = duplicate_donor_snp_info
        ),
        "duplicate"
    )
})

test_that("SNPData() allows the same (snp_id, donor) pair once per distinct zygosity_source", {
    two_source_donor_snp_info <- data.frame(
        snp_id = c("snp_1", "snp_1"),
        donor = c("donor_1", "donor_1"),
        zygosity = c("het", "hom"),
        zygosity_source = c("vireo_gt", "binomial")
    )

    snp_data <- SNPData(
        ref_count = test_ref_count,
        alt_count = test_alt_count,
        snp_info = test_snp_info,
        barcode_info = test_barcode_info,
        donor_snp_info = two_source_donor_snp_info
    )

    # Verify both source-tagged rows for the same pair are retained
    expect_equal(nrow(donor_snp_info(snp_data, source = "all")), 2)
})

test_that("SNPData() throws error when donor_snp_info has duplicate (snp_id, donor, zygosity_source) rows", {
    duplicate_donor_snp_info <- data.frame(
        snp_id = c("snp_1", "snp_1"),
        donor = c("donor_1", "donor_1"),
        zygosity = c("het", "hom"),
        zygosity_source = c("vireo_gt", "vireo_gt")
    )

    # Verify validity rejects a duplicate key once zygosity_source is part of it
    expect_error(
        SNPData(
            ref_count = test_ref_count,
            alt_count = test_alt_count,
            snp_info = test_snp_info,
            barcode_info = test_barcode_info,
            donor_snp_info = duplicate_donor_snp_info
        ),
        "duplicate"
    )
})

test_that("SNPData() throws error when donor_snp_info references an unknown snp_id", {
    unknown_snp_donor_info <- data.frame(
        snp_id = "snp_does_not_exist",
        donor = "donor_1",
        zygosity = "het"
    )

    # Verify validity rejects a snp_id absent from snp_info
    expect_error(
        SNPData(
            ref_count = test_ref_count,
            alt_count = test_alt_count,
            snp_info = test_snp_info,
            barcode_info = test_barcode_info,
            donor_snp_info = unknown_snp_donor_info
        ),
        "not present in snp_info"
    )
})

test_that("SNPData() throws error when donor_snp_info references an unknown donor", {
    unknown_donor_snp_info <- data.frame(
        snp_id = "snp_1",
        donor = "donor_does_not_exist",
        zygosity = "het"
    )

    # Verify validity rejects a donor absent from donor_info
    expect_error(
        SNPData(
            ref_count = test_ref_count,
            alt_count = test_alt_count,
            snp_info = test_snp_info,
            barcode_info = test_barcode_info,
            donor_snp_info = unknown_donor_snp_info
        ),
        "not present in donor_info"
    )
})

test_that("SNPData() throws error when a zygosity call has no zygosity_source", {
    unsourced_zygosity <- data.frame(
        snp_id = "snp_1",
        donor = "donor_1",
        zygosity = "het",
        zygosity_source = NA_character_
    )

    # Verify validity requires zygosity_source whenever zygosity is called
    expect_error(
        SNPData(
            ref_count = test_ref_count,
            alt_count = test_alt_count,
            snp_info = test_snp_info,
            barcode_info = test_barcode_info,
            donor_snp_info = unsourced_zygosity
        ),
        "zygosity call but no zygosity_source"
    )
})

test_that("zygosity_source() defaults to NA when no donor_snp_info is supplied", {
    snp_data <- create_test_snpdata()

    # Verify a fresh object with no genotype calls has no active source
    expect_true(is.na(zygosity_source(snp_data)))
})

test_that("zygosity_source() is auto-derived at construction from a single-source donor_snp_info", {
    snp_data <- SNPData(
        ref_count = test_ref_count,
        alt_count = test_alt_count,
        snp_info = test_snp_info,
        barcode_info = test_barcode_info,
        donor_snp_info = data.frame(
            snp_id = "snp_1",
            donor = "donor_1",
            zygosity = "het",
            zygosity_source = "vireo_gt"
        )
    )

    # Verify the single distinct source is picked up automatically
    expect_equal(zygosity_source(snp_data), "vireo_gt")
})

test_that("zygosity_source() stays NA at construction when donor_snp_info carries more than one source", {
    snp_data <- SNPData(
        ref_count = test_ref_count,
        alt_count = test_alt_count,
        snp_info = test_snp_info,
        barcode_info = test_barcode_info,
        donor_snp_info = data.frame(
            snp_id = c("snp_1", "snp_1"),
            donor = c("donor_1", "donor_1"),
            zygosity = c("het", "hom"),
            zygosity_source = c("vireo_gt", "binomial")
        )
    )

    # Verify ambiguity between sources is left unresolved rather than guessed
    expect_true(is.na(zygosity_source(snp_data)))
})

test_that("zygosity_source<- switches the active source and validates against what's present", {
    snp_data <- SNPData(
        ref_count = test_ref_count,
        alt_count = test_alt_count,
        snp_info = test_snp_info,
        barcode_info = test_barcode_info,
        donor_snp_info = data.frame(
            snp_id = c("snp_1", "snp_1"),
            donor = c("donor_1", "donor_1"),
            zygosity = c("het", "hom"),
            zygosity_source = c("vireo_gt", "binomial")
        )
    )

    zygosity_source(snp_data) <- "binomial"
    # Verify the active source switched
    expect_equal(zygosity_source(snp_data), "binomial")

    # Verify switching to a source not present in donor_snp_info errors
    expect_error(
        zygosity_source(snp_data) <- "bulk_wgs",
        "not found in donor_snp_info"
    )
})

test_that("donor_snp_info() filters to the active source by default", {
    snp_data <- SNPData(
        ref_count = test_ref_count,
        alt_count = test_alt_count,
        snp_info = test_snp_info,
        barcode_info = test_barcode_info,
        donor_snp_info = data.frame(
            snp_id = c("snp_1", "snp_1"),
            donor = c("donor_1", "donor_1"),
            zygosity = c("het", "hom"),
            zygosity_source = c("vireo_gt", "binomial")
        )
    )
    # Two sources are present, so construction leaves the active source
    # unresolved (see the ambiguity test above); set it explicitly.
    zygosity_source(snp_data) <- "vireo_gt"

    # Verify the default view only surfaces the active source's row
    default_view <- donor_snp_info(snp_data)
    expect_equal(nrow(default_view), 1)
    expect_equal(default_view$zygosity_source, "vireo_gt")

    # Verify source = "all" surfaces every row regardless of the active source
    expect_equal(nrow(donor_snp_info(snp_data, source = "all")), 2)

    # Verify naming a specific source filters to just that one
    binomial_view <- donor_snp_info(snp_data, source = "binomial")
    expect_equal(nrow(binomial_view), 1)
    expect_equal(binomial_view$zygosity_source, "binomial")
})

test_that("updateObject() derives zygosity_source for a legacy object missing the slot", {
    snp_data <- SNPData(
        ref_count = test_ref_count,
        alt_count = test_alt_count,
        snp_info = test_snp_info,
        barcode_info = test_barcode_info,
        donor_snp_info = data.frame(
            snp_id = "snp_1",
            donor = "donor_1",
            zygosity = "het",
            zygosity_source = "vireo_gt"
        )
    )
    # Simulate a legacy object serialised before the zygosity_source slot
    # existed, same technique as the packed *_by_donor migration test below
    legacy_attrs <- attributes(snp_data)
    legacy_attrs$zygosity_source <- NULL
    attributes(snp_data) <- legacy_attrs

    updated <- updateObject(snp_data)

    # Verify the slot was derived from the existing donor_snp_info on update
    expect_equal(zygosity_source(updated), "vireo_gt")
})

test_that("[() drops a donor's donor_info and donor_snp_info rows once its cells are subsetted out", {
    two_donor_barcode_info <- data.frame(
        cell_id = c("cell_1", "cell_2"),
        donor = c("donor_1", "donor_2"),
        stringsAsFactors = FALSE
    )
    donor_snp_info <- data.frame(
        snp_id = c("snp_1", "snp_2"),
        donor = c("donor_1", "donor_2"),
        zygosity = c("het", "het"),
        zygosity_source = c("vireo_gt", "vireo_gt")
    )
    snp_data <- SNPData(
        ref_count = test_ref_count,
        alt_count = test_alt_count,
        snp_info = test_snp_info,
        barcode_info = two_donor_barcode_info,
        donor_snp_info = donor_snp_info
    )

    subset_data <- snp_data[, 1]

    # Verify only the surviving cell's donor remains in donor_info
    expect_equal(donor_info(subset_data)$donor, "donor_1")
    # Verify the dropped donor's donor_snp_info row was removed too
    expect_equal(donor_snp_info(subset_data)$donor, "donor_1")
})

test_that("updateObject() migrates a legacy object's packed *_by_donor snp_info columns", {
    snp_data <- create_test_snpdata()
    legacy_snp_info <- snp_info(snp_data)
    legacy_snp_info$xci_informative <- c(TRUE, FALSE)
    legacy_snp_info$xci_informative_donor <- c("donor_1", NA)
    legacy_snp_info$xci_allele_on_x1_by_donor <- c("REF", NA)
    legacy_snp_info$xci_escape_fraction_by_donor <- c("0.02", NA)
    snp_data@snp_info <- legacy_snp_info

    # Strip the slots this object predates, simulating a deserialized legacy object
    legacy_attrs <- attributes(snp_data)
    legacy_attrs$donor_info <- NULL
    legacy_attrs$donor_snp_info <- NULL
    legacy_attrs$chr_style <- NULL
    attributes(snp_data) <- legacy_attrs

    updated <- updateObject(snp_data)

    # Verify the packed column was unpacked into donor_snp_info
    expect_equal(donor_snp_info(updated)$snp_id, "snp_1")
    # Verify the phase carried over correctly
    expect_equal(donor_snp_info(updated)$allele_on_x1, "REF")
    # Verify the packed columns were removed from snp_info
    expect_false("xci_informative_donor" %in% colnames(snp_info(updated)))
    # Verify xci_informative itself was migrated out of snp_info too, since it
    # now belongs solely to donor_snp_info (informativeness is donor-specific)
    expect_false("xci_informative" %in% colnames(snp_info(updated)))
    # Verify the object passes validity after migration
    expect_true(methods::validObject(updated))
})

test_that("SNPData() donor_map relabels barcode_info$donor at construction time", {
    snp_data <- SNPData(
        ref_count = test_ref_count,
        alt_count = test_alt_count,
        snp_info = test_snp_info,
        barcode_info = test_barcode_info,
        donor_map = c(PatientA = "donor_1")
    )

    # Verify the raw donor label was relabelled before donor_info was derived
    expect_equal(unique(barcode_info(snp_data)$donor), "PatientA")
    expect_equal(donor_info(snp_data)$donor, "PatientA")
})

test_that("SNPData() donor_map also relabels a caller-supplied donor_snp_info", {
    donor_snp_info <- data.frame(
        snp_id = "snp_1",
        donor = "donor_1",
        zygosity = "het",
        zygosity_source = "vireo_gt"
    )

    snp_data <- SNPData(
        ref_count = test_ref_count,
        alt_count = test_alt_count,
        snp_info = test_snp_info,
        barcode_info = test_barcode_info,
        donor_snp_info = donor_snp_info,
        donor_map = c(PatientA = "donor_1")
    )

    # Verify donor_snp_info was relabelled to match barcode_info/donor_info
    expect_equal(donor_snp_info(snp_data)$donor, "PatientA")
})

test_that("rename_donor() cascades a relabel across barcode_info, donor_info, and donor_snp_info", {
    snp_data <- create_test_snpdata()
    snp_data <- add_donor_snp_metadata(
        snp_data,
        data.frame(snp_id = "snp_1", donor = "donor_1", zygosity = "het", zygosity_source = "vireo_gt")
    )

    renamed <- rename_donor(snp_data, c(PatientA = "donor_1"))

    # Verify all three donor-keyed tables agree on the new label
    expect_equal(unique(barcode_info(renamed)$donor), "PatientA")
    expect_equal(donor_info(renamed)$donor, "PatientA")
    expect_equal(donor_snp_info(renamed)$donor, "PatientA")
})

test_that("rename_donor() swaps two donor labels without corrupting either", {
    barcode_info <- data.frame(cell_id = c("cell_1", "cell_2"), donor = c("A", "B"))
    snp_data <- SNPData(
        ref_count = test_ref_count,
        alt_count = test_alt_count,
        snp_info = test_snp_info,
        barcode_info = barcode_info
    )

    # new = old for each entry: donor A's new label is B, donor B's new label is A
    swapped <- rename_donor(snp_data, c(B = "A", A = "B"))

    # Verify the swap applied simultaneously rather than cascading sequentially
    # (a naive pair-by-pair substitution would collapse both cells onto "A")
    expect_equal(barcode_info(swapped)$donor, c("B", "A"))
})

test_that("rename_donor() errors when donor_map references a donor that doesn't exist", {
    snp_data <- create_test_snpdata()

    # Verify the error names the unknown (old) donor
    expect_error(
        rename_donor(snp_data, c(X = "no_such_donor")),
        "no_such_donor"
    )
})

test_that("rename_donor() errors when donor_map would collide two donors onto one label", {
    barcode_info <- data.frame(cell_id = c("cell_1", "cell_2"), donor = c("A", "B"))
    snp_data <- SNPData(
        ref_count = test_ref_count,
        alt_count = test_alt_count,
        snp_info = test_snp_info,
        barcode_info = barcode_info
    )

    # Verify renaming A to the already-existing label B is rejected, not merged
    expect_error(
        rename_donor(snp_data, c(B = "A")),
        "duplicate donor label"
    )
})

test_that("barcode_info<- rejects a donor column change once donor_info carries data", {
    snp_data <- create_test_snpdata()
    new_barcode_info <- barcode_info(snp_data)
    new_barcode_info$donor <- "donor_2"

    # Verify the error points at rename_donor() instead
    expect_error(
        barcode_info(snp_data) <- new_barcode_info,
        "rename_donor"
    )
})

test_that("barcode_info<- still allows non-donor columns to change once donor_info carries data", {
    snp_data <- create_test_snpdata()
    new_barcode_info <- barcode_info(snp_data)
    new_barcode_info$clonotype <- "clonotype_new"

    barcode_info(snp_data) <- new_barcode_info

    # Verify the unrelated column change went through
    expect_equal(barcode_info(snp_data)$clonotype, c("clonotype_new", "clonotype_new"))
})

test_that("barcode_info<- allows a donor column to be set freely when donor_info is still empty", {
    no_donor_barcode_info <- data.frame(cell_id = c("cell_1", "cell_2"))
    snp_data <- SNPData(
        ref_count = test_ref_count,
        alt_count = test_alt_count,
        snp_info = test_snp_info,
        barcode_info = no_donor_barcode_info
    )
    new_barcode_info <- barcode_info(snp_data)
    new_barcode_info$donor <- c("donor_a", "donor_b")

    barcode_info(snp_data) <- new_barcode_info

    # Verify the donor column was accepted since donor_info had no prior data
    expect_equal(barcode_info(snp_data)$donor, c("donor_a", "donor_b"))
})

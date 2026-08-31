# snplet 0.6.5

* Added a `test_escape()` method taking a SNPData object directly, which draws the counts,
  the null escape fraction and the overdispersion from the donor's own fit rather than
  requiring them to be supplied
* Added gene-level pooling of read-backed phase blocks, so a molecule spanning several of a
  gene's heterozygous SNPs is counted once rather than once per SNP
* Added a `snp_gene_map` slot to SNPData, populated by `import_cellsnp()` when the gene
  annotation carries a `strand` column, giving each SNP-to-gene candidate its own row for
  molecule attribution at multi-gene SNPs
* Changed `haplotype_expression()` to report one row per gene by default, electing a
  representative SNP, with the per-SNP output available via `by_snp = TRUE`
* Changed `import_cellsnp()` to check BAM files when given, and moved `library_id` earlier in
  the argument list
* Fixed the beta-binomial CDF being evaluated outside its support during escape testing
* Removed unused imports flagged by roxygen2 8.0

# snplet 0.6.4

* Added `library_id` to SNPData, so an object can hold cells from more than one sequencing
  library without conflating barcodes that are only unique within a library
* Added a `library_info` slot recording the BAM file(s) each library's reads came from, along
  with `add_library_bams()`, so molecule-level counting can find them without being told again
* Changed `merge_snpdata()` behaviour when barcodes overlap between objects
* Documented that phasing is derived from RNA-seq expression only, and the limits that places
  on interpreting active/inactive haplotypes

# snplet 0.6.3

* Added a hard cap of 0.5 to the escape fraction reported for non-XIST genes
* Changed metadata updates to modify slots in place rather than rebuilding the whole object

# snplet 0.6.2

* Added `test_escape()`, an exact beta-binomial test of the phase-corrected inactive-X count
  against a null escape fraction
* Added `infer_zygosity()` for binomial-derived heterozygosity calls, and refactored
  `zygosity_source` so calls from different sources coexist on one object
* Added faster SNPData merging, with reworked sparse matrix subsetting
* Added logging of XCI model fitting
* Changed escape testing from a fixed to a donor-estimated beta-binomial overdispersion
* Changed the canonical accessors to drop the `get_` prefix, retaining the old names as
  aliases
* Fixed R CMD check issues

# snplet 0.6.1

* Added strand handling to molecule-level counting, so reads reoriented by demultiplexing are
  attributed to the correct gene
* Added an assertion that every EM observation carries a cell posterior
* Changed `aggregate_count_df()` to accept more than one grouping column
* Changed inactive/active haplotype terminology for consistency
* Refactored the XCI beta-binomial EM into a shared engine with its own diagnostics, used by
  both the cell- and clonotype-level entry points

# snplet 0.6.0

* Added read-backed molecule-level phasing: `extract_snp_calls()` for per-(molecule, SNP)
  allele calls from an indexed BAM, `molecule_snp_alleles()`, `phase_snps()` for phasing
  heterozygous SNPs from the molecules that span them, and `add_molecule_phase()` to orient
  the resulting blocks against the EM phase and store them on the object
* Added `haplotype_expression_by_molecule()`, which counts each molecule once per gene instead
  of once per SNP it covers
* Added exclusion of `"doublet"` and `"unassigned"` donors from `assign_xci()`, neither of
  which has a genotype to fit against

# snplet 0.5.3

* Added `donor_info` and `donor_snp_info` slots to SNPData, with `add_donor_metadata()` and
  `add_donor_snp_metadata()`, and support for both in `merge_snpdata()`
* Added import of Vireo genotypes: `import_cellsnp()` now takes `vireo_folder` in place of
  `vireo_file` and reads `GT_donors.vireo.vcf.gz` when present to populate per-(SNP, donor)
  zygosity, with `donor_map` for relabelling Vireo's arbitrary donor names at import
* Added `zygosity_source` to key zygosity calls by where they came from
* Added `haplotype_expression()` for per-gene active and inactive haplotype counts
* Added `updateObject()` to backfill `chr_style` on objects made by earlier versions
* Changed inactive-X terminology to XCI throughout (`assign_xci()`,
  `assign_xci_by_clonotype()`, `plot_xci_heatmap()`)
* Changed spelling to British English throughout prose, messages and identifiers
* Removed `refit_after_filter()`, and made `test_maf()`, `detect_chr_style()` and
  `convert_chr_style()` internal

# snplet 0.5.2

* Added a test data construction helper and substantially expanded the test suite
* Changed the default `min_cov` to 1
* Fixed refitting when no escapee genes are present
* Refactored the cell- and clonotype-level assignment paths onto shared code

# snplet 0.5.1

* Fixed inactive X assignment to use only X chromosome SNPs
* Added `@family` tags across the documentation

# snplet 0.5.0

* Optimised inactive X assignment

# snplet 0.4.5

* Added `xci_haplotypes()` for the model's inferred per-SNP phase and escape fraction
* Refactored the binomial test behind a single wrapper

# snplet 0.4.3

* Changed the XCI heatmap to label rows with gene names

# snplet 0.4.2

* Added NA removal to expression matrix construction
* Added checks to grouped row sums
* Changed random seed handling to local scope, so fitting no longer disturbs the global RNG
  state
* Moved ComplexHeatmap from Suggests to Imports
* Refactored count aggregation to remove duplication between the `*_count_df()` functions

# snplet 0.4.1

* Fixed a possible division by zero in the MAD calculation
* Fixed gene filtering during assignment

# snplet 0.4.0

* Fixed gene filtering when refitting the assignment model

# snplet 0.3.35

* Reimplemented the inactive X assignment algorithm
* Added parallelism via furrr/future
* Fixed missing gene filtering, and updated the log-likelihood calculation
* Removed unused heatmap plotting functions

# snplet 0.3.30

* Fixed a filter expression check that wrongly accepted environment variables in place of
  column names
* Improved robustness of inactive X assignment and its clustering step
* Reduced the volume of logging emitted during assignment

# snplet 0.3.29

* Improved test coverage across import, filtering, SNPData methods, `test_maf()`,
  `to_expr_matrix()` and the chromosome name utilities
* Changed inactive X assignment to use the raw expression matrix value rather than its sign

# snplet 0.3.28

* Changed aggregated SNP counts to carry gene name and chromosome position

# snplet 0.3.27

* Added `get_donor_het_snpdata()` to pull the heterozygous SNPs of a single donor
* Added clonotype-level inactive X assignment
* Removed the requirement that a barcode metadata update table cover every `cell_id`

# snplet 0.3.26

* Reimplemented `assign_inactive_x()` on a simpler basis
* Changed the expression matrix calculation

# snplet 0.3.25

* Changed the cell cap from random sampling to taking the highest-library-size cells
* Reverted to the previous inactive X assignment after a regression

# snplet 0.3.24

* Changed SNP filtering during assignment from a quantile to a fixed threshold

# snplet 0.3.23

* Added clonotype-level inactive X assignment
* Added removal of constant rows and columns before fitting

# snplet 0.3.22

* Fixed `assign_inactive_x()` to work in the presence of multiple donors

# snplet 0.3.20

* Fixed `assign_inactive_x()` to use only heterozygous sites

# snplet 0.3.19

* Added upper limits on the number of cells and SNPs used for X inactivation assignment
* Fixed the canonical X chromosome name, and renamed `chr_canonical` to `chrom_canonical`

# snplet 0.3.18

* Added `assign_inactive_x()`
* Added handling for differing chromosome naming conventions, recorded as `chr_style` on the
  object, with backwards compatibility for objects lacking it
* Fixed a bad logger threshold call

# snplet 0.3.16

* Renamed `assign_snp_genes()` to `add_snp_gene_names()`

# snplet 0.3.15

* Added `assign_snp_genes()` for mapping SNPs onto overlapping genes
* Refactored `test_maf()` onto a vectorised `pbeta()` call instead of mapping `binom.test()`
* Refactored the SNP zygosity calculation

# snplet 0.3.14

* Added `donor_het_status_df()`, assigning a zygosity to each SNP for each donor
* Improved SNPData merging efficiency
* Fixed validation when reassigning metadata

# snplet 0.3.13

* Changed merging to key on barcodes rather than cell IDs
* Improved handling of duplicate SNPs

# snplet 0.3.12

* Fixed handling of duplicate SNPs in `import_cellsnp()`

# snplet 0.3.10

* Added `merge_snpdata()` for combining SNPData objects
* Changed how `snp_id` is generated
* Fixed R CMD check issues

# snplet 0.3.9

* Added metadata update utilities (`add_barcode_metadata()`, `add_snp_metadata()`) and setters
  for the `barcode_info` and `snp_info` slots
* Added handling for missing VDJ information when exporting a SNPData object, and clearer
  errors for missing clonotype data
* Changed `vdj_file` to be optional when constructing a SNPData object

# snplet 0.3.7

* Fixed the cellSNP import logic for column name correction
* Fixed filter variable checking so function arguments are resolved correctly

# snplet 0.3.6

* Added a `.lintr` configuration and applied Air formatting across the package
* Fixed missing exports

# snplet 0.3.5

* Added VCF handling: the `VCFData` class and `read_vcf()`

# snplet 0.3.1

* Removed the `ref_fraction()` and `major_allele_fraction()` methods

# snplet 0.3.0

* Addressed BiocCheck notes

# snplet 0.2.0

* Added `aggregate_count_df()` for aggregating counts by any `barcode_info` column
* Added GitHub Actions for checks and coverage
* Changed terminology from sample to barcode throughout the codebase, renaming
  `get_sample_info()` to `get_barcode_info()` with a backwards-compatible alias
* Fixed importing without Vireo files
* Fixed barcode export and namespacing

# snplet 0.1.0

* Added OTH count import, giving the object its third count matrix
* Added bundled example data and `get_example_snpdata()`
* Added `export_cellsnp()`, and changed `vireo_file` to be optional
* Refactored the SNPData class and its methods
* Removed the redundant `cell_count_df()` function

# snplet 0.0.7

* Added `to_expr_matrix()` for converting allele counts to an expression-like matrix
* Added a barcode-level count function
* Added input validation for gene annotation, a variable checker for SNP filtering, and an
  error when filtering removes all data
* Added the MAF distribution plot
* Changed the MAF test to be more conservative

# snplet 0.0.6

* Added `import_cellsnp()` for reading cellSNP-lite output
* Added the filtering functions, including doublet removal
* Parallelised `test_maf()`

# snplet 0.0.5

* Added the first tests
* Changed matrix index operations to never drop dimensions
* Changed the `"unassigned"` and `"doublet"` columns to be dropped only when present
* Fixed an incorrect `is.na()` call

# snplet 0.0.1

* Initial version: the `SNPData` S4 class holding REF/ALT count matrices with SNP and cell
  metadata, its accessors, and coverage calculation

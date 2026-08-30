# snplet

An R package for single-cell allele-specific expression analysis. `snplet` integrates SNP
count matrices from cellSNP-lite with donor assignments from Vireo, optional clonotype
annotations from cellranger VDJ (TCR/BCR-seq), and optional BAM files for molecule-level
counting, to quantify allelic imbalance and X-chromosome inactivation (XCI) across individual
cells.

## Inputs

| Input | Source | Required? |
| --- | --- | --- |
| REF/ALT/OTH count matrices, SNP and barcode lists | [cellSNP-lite](https://github.com/single-cell-genetics/cellSNP-lite) | Required |
| Gene annotation (`chrom`, `start`, `end`, `gene_name`, optionally `strand`) | Any annotation source | Required |
| Donor assignments, and per-donor genotypes if `GT_donors.vireo.vcf.gz` is present | [Vireo](https://github.com/single-cell-genetics/vireo) | Optional |
| Clonotype annotations (`filtered_contig_annotations.csv`) | cellranger VDJ (TCR/BCR-seq) | Optional |
| Aligned reads for molecule-level counting | Indexed BAM file(s) | Optional |

Only cellSNP-lite output and a gene annotation are needed to build a `SNPData` object.
Vireo output adds donor identity and genotype-derived zygosity; **TCR/BCR-seq (VDJ) data is
optional** and is used purely to group cells by clonotype for aggregation and XCI assignment —
the package does not perform clonal lineage or clonal evolution analysis. **BAM files are also
optional**, and unlock read-backed phasing and molecule-level counts (see
[Molecule-level counting](#molecule-level-counting-optional)). Each analysis layer only requires
the inputs it uses, so an object built from cellSNP-lite alone still supports filtering,
aggregation, MAF testing, and expression-matrix export.

## Features

### Core data structures
- **`SNPData` S4 class**: REF/ALT/OTH sparse count matrices (SNPs x cells) with SNP metadata
  (`snp_info`), cell metadata (`barcode_info`), per-donor metadata (`donor_info`),
  per-(SNP, donor) zygosity and XCI diagnostics (`donor_snp_info`), per-library metadata
  including recorded BAM paths (`library_info`), and a SNP-to-gene map (`snp_gene_map`).
- **`VCFData` S4 class**: lightweight VCF container (`read_vcf()`), with `header()`,
  `samples()`, and `variants()` accessors.
- **Automatic metrics**: coverage, library size, and minor allele frequency are recomputed as
  the object is filtered and merged.

### Import and export
- `import_cellsnp()` builds a `SNPData` object from a cellSNP-lite directory, taking optional
  `vdj_file`, `vireo_folder`, and `bam_files` arguments.
- `merge_snpdata()` combines runs, using each object's `library_id` to keep cells that share a
  10x barcode by chance apart.
- `export_cellsnp()` writes an object back out in cellSNP-lite format; `to_expr_matrix()`
  converts allele counts into an expression-like matrix.

### Processing and aggregation
- **Filtering**: `filter_snps()`, `filter_barcodes()` (`filter_samples()` is a
  backwards-compatible alias), `remove_doublets()`, `remove_na_clonotypes()`,
  `remove_na_genes()`.
- **Aggregation**: `barcode_count_df()`, `donor_count_df()`, `clonotype_count_df()`, and
  `aggregate_count_df()` for any `barcode_info` column, each with an optional exact binomial
  test of allele usage against a null minor allele frequency.
- **Metadata**: `add_barcode_metadata()`, `add_snp_metadata()`, `add_donor_metadata()`,
  `add_snp_gene_names()`, `assign_snp_genes()`, `add_library_bams()`, `rename_donor()`.

### Zygosity
- Vireo genotypes populate per-(SNP, donor) zygosity at import when `GT_donors.vireo.vcf.gz` is
  present.
- `infer_zygosity()` derives binomial zygosity calls from the counts themselves when genotypes
  are unavailable; `donor_het_status_df()` and `get_donor_het_snpdata()` inspect and subset on
  them. Calls are keyed by source, and `zygosity_source()` selects which is active.

### X-chromosome inactivation
- `assign_xci()` fits a beta-binomial EM model per donor to call the active X per cell;
  `assign_xci_by_clonotype()` runs the same model at clonotype level (requires VDJ data) and
  projects assignments back to cells.
- `xci_assignments()` and `xci_haplotypes()` extract the stored per-cell calls and per-SNP
  inferred phase.
- `haplotype_expression()` reports active- and inactive-haplotype counts per gene, and
  `test_escape()` tests each gene for escape with an exact beta-binomial test against the
  donor's own fitted null and overdispersion.

### Molecule-level counting (optional)
Supplying indexed BAM files enables read-backed phasing and per-molecule counting, which avoids
counting one read once per SNP it covers:
- `extract_snp_calls()` reads per-(molecule, SNP) allele calls from a BAM, using threaded
  `samtools` when available and falling back to Rsamtools/GenomicAlignments otherwise.
- `phase_snps()` phases heterozygous SNPs from the molecules that span them, and
  `add_molecule_phase()` orients the resulting blocks against the EM phase and writes them into
  the object.
- `haplotype_expression_by_molecule()` then counts each molecule once per gene, which
  `test_escape()` prefers automatically when present.

## Installation

```r
# Install from GitHub
devtools::install_github("Shians/snplet")
```

Upstream tools (cellSNP-lite, Vireo, cellranger VDJ) must be run separately to produce the input
files `snplet` imports. [samtools](http://www.htslib.org/) on `PATH` is optional and accelerates
molecule-level read extraction.

## Quick start

```r
library(snplet)

# Load a small bundled example SNPData object
snp_data <- get_example_snpdata()

# Filter low-coverage SNPs and cells
snp_data <- filter_snps(snp_data, coverage > 5)
snp_data <- filter_barcodes(snp_data, library_size > 100)

# Aggregate allele counts by donor; includes a binomial test of allele usage
# against a null minor allele frequency (p_val / adj_p_val columns)
donor_counts <- donor_count_df(snp_data)
```

## Importing your own data

```r
# cellSNP-lite output alone is enough
snp_data <- import_cellsnp(
  cellsnp_dir = "path/to/cellsnp_output",
  gene_annotation = gene_anno_df,
  library_id = "run1"
)

# Add donor assignments, TCR/BCR clonotypes and BAM paths, all optional
snp_data <- import_cellsnp(
  cellsnp_dir = "path/to/cellsnp_output",
  gene_annotation = gene_anno_df,
  library_id = "run1",
  vireo_folder = "path/to/vireo_output",
  vdj_file = "path/to/filtered_contig_annotations.csv",
  bam_files = "path/to/possorted_genome_bam.bam"
)

# Multiple libraries merge on distinct library_id labels
combined <- merge_snpdata(run1, run2)
```

## XCI and escape workflow

```r
# Fit the active-X model per donor
snp_data <- assign_xci(snp_data)

# Inspect assignments
xci_assignments(snp_data)
plot_xci_heatmap(snp_data, donor = "donor1")

# Per-gene active/inactive haplotype counts, and a test for escape
hap <- haplotype_expression(snp_data)
escape <- test_escape(snp_data)
```

If BAM files were recorded at import (or added later with `add_library_bams()`), read-backed
phase and molecule-level counts can be added before testing:

```r
snp_data <- add_molecule_phase(snp_data)
hap_mol <- haplotype_expression_by_molecule(snp_data)

# test_escape() uses the molecule counts automatically when they are present
escape <- test_escape(snp_data)
```

See the function documentation (`?import_cellsnp`, `?assign_xci`, `?haplotype_expression`,
`?add_molecule_phase`, `?test_escape`) for full details.

## License

This project is licensed under the Apache License (>= 2).

# snplet

An R package for single-cell SNP data analysis, focusing on allele-specific expression. Integrates SNP count matrices (REF/ALT), donor assignments from Vireo, and clonotype metadata from cellranger VDJ to support comprehensive allele-specific analysis workflows.

## Features

### Core Data Structure
- **SNPData S4 Class**: Central object containing REF/ALT count matrices, SNP metadata, and cell metadata
- **Automatic metrics**: Coverage, library size, and minor allele frequency calculations

### Data Processing
- **Quality filtering**: Filter SNPs and cells by coverage, library size, and metadata
- **Data aggregation**: Summarize counts at cell, donor, and clonotype levels

## Installation

```r
# Install from GitHub
devtools::install_github("Shians/snplet")
```

`snplet` reads output from [cellSNP-lite](https://github.com/single-cell-genetics/cellSNP-lite)
(SNP genotype calls), [Vireo](https://github.com/single-cell-genetics/vireo) (donor
assignments), and cellranger VDJ (clonotype annotations). These upstream tools must be run
separately to produce the input files `snplet` imports. `extract_snp_calls()` can optionally use
[samtools](http://www.htslib.org/) on `PATH` to accelerate molecule-level read extraction; it
falls back to pure R (Rsamtools/GenomicAlignments) if samtools is not available.

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

For a full workflow — importing raw cellSNP-lite/Vireo/cellranger output with
`import_cellsnp()`, filtering, aggregation, and X-chromosome inactivation or haplotype
expression analysis — see the function documentation (`?import_cellsnp`, `?assign_xci`,
`?haplotype_expression`).

## License

This project is licensed under the Apache License (>= 2).

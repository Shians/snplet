# snplet

**This package is currently work-in-progress. The functions are expected to change dramatically throughout development and nothing should be considered stable.**

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
devtools::install_github("shians/snplet")
```

## License

This project is licensed under the Apache License (>= 2).

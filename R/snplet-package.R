#' snplet: Allele-Specific Expression Analysis for Single-Cell RNA-Seq
#'
#' \code{snplet} enables downstream analysis of allele-specific expression in single-cell RNA
#' sequencing data, using output from cellSNP-lite, Vireo, and cellranger VDJ. It integrates SNP
#' genotype calls, donor assignments, and clonotype annotations to assess allelic imbalance
#' across individual cells.
#'
#' @section Typical workflow:
#' \enumerate{
#'   \item Import data with \code{\link{import_cellsnp}} (optionally merge multiple runs with
#'     \code{\link{merge_snpdata}}).
#'   \item Filter SNPs and cells with \code{\link{filter_snps}} and \code{\link{filter_barcodes}}.
#'   \item Aggregate counts with \code{\link{barcode_count_df}}, \code{\link{donor_count_df}}, or
#'     \code{\link{clonotype_count_df}}, and visualise with the \code{plot_*} functions.
#'   \item Analyse allele usage (\code{\link{test_maf}}), X-chromosome inactivation
#'     (\code{\link{assign_xci}}), or haplotype expression (\code{\link{haplotype_expression}}).
#' }
#'
#' \code{\link{get_example_snpdata}()} loads a small bundled dataset useful for trying out the
#' package without any external files.
#'
#' @keywords internal
"_PACKAGE"

#' @importFrom Matrix Matrix rowSums colSums readMM writeMM sparseMatrix
#' @importFrom dplyr mutate inner_join filter left_join select rename distinct
#' @importFrom dplyr first summarise if_else count arrange desc slice any_of pull
#' @importFrom dplyr slice_head slice_sample bind_rows case_when group_by
#' @importFrom dplyr relocate all_of
#' @importFrom tidyr pivot_longer contains everything
#' @importFrom tibble as_tibble tibble
#' @importFrom rlang enquos get_expr sym
#' @importFrom scales percent label_comma label_number cut_short_scale
#' @importFrom logger log_info log_success log_warn
#' @importFrom glue glue
#' @importFrom stringr str_remove str_to_title str_detect str_split_fixed
#' @importFrom methods as is new setClass setGeneric setMethod show
#' @importFrom BiocGenerics nrow ncol rownames colnames updateObject start end
#' @importFrom magrittr %>%
#' @importFrom ggplot2 ggplot aes geom_density scale_y_continuous scale_x_log10
#' @importFrom ggplot2 geom_segment scale_x_continuous theme_void theme
#' @importFrom ggplot2 element_blank element_line element_text margin labs
#' @importFrom ggplot2 geom_hline geom_point facet_grid vars expansion
#' @importFrom ggplot2 coord_cartesian theme_classic theme_bw
#' @importFrom ggrepel geom_text_repel
#' @importFrom readr read_tsv read_csv write_tsv write_csv cols col_character col_integer
#' @importFrom plyranges as_granges join_overlap_left
#' @importFrom fs path
#' @importFrom furrr future_map future_map2
#' @importFrom purrr map map2 map2_dbl map_dbl
#' @importFrom R.utils gzip
#' @importFrom grid unit
#' @importFrom VGAM dbetabinom pbetabinom
#' @importFrom stats median quantile p.adjust pbeta setNames cor as.dist hclust cutree dist mad na.omit
#' @importFrom utils head combn write.table
#' @importFrom Rsamtools ScanBamParam scanBamFlag
#' @importFrom GenomicAlignments readGAlignments cigar cigarRangesAlongReferenceSpace
#' @importFrom GenomicAlignments cigarRangesAlongQuerySpace
#' @importFrom GenomicRanges reduce
#' @importFrom Seqinfo seqnames seqlevels "seqlevels<-"
#' @importFrom IRanges IRanges findOverlaps PartitioningByWidth togroup
#' @importFrom S4Vectors mcols queryHits subjectHits
#' @importFrom Biostrings subseq
NULL

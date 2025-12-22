#' Assign overlapping gene names to SNPs
#'
#' Takes SNP coordinates and gene annotations, then adds a gene_name column
#' to each SNP with all overlapping genes collapsed by ", ".
#'
#' @param snp_df Data frame with columns \code{chrom} and \code{pos}.
#' @param gene_anno Data frame with columns \code{chrom}, \code{start},
#'   \code{end}, and \code{gene_name}.
#'
#' @return A data frame matching \code{snp_df} with an added \code{gene_name} column.
#'
#' @examples
#' snp_df <- data.frame(chrom = c("chr1", "chr1"), pos = c(100, 500))
#' gene_anno <- data.frame(
#'     chrom = c("chr1", "chr1"),
#'     start = c(50, 450),
#'     end = c(150, 550),
#'     gene_name = c("GENE1", "GENE2")
#' )
#' assign_snp_genes(snp_df, gene_anno)
#'
#' @export
assign_snp_genes <- function(snp_df, gene_anno) {
    required_snp_cols <- c("chrom", "pos")
    missing_snp_cols <- setdiff(required_snp_cols, colnames(snp_df))
    if (length(missing_snp_cols) > 0) {
        stop(sprintf(
            "snp_df is missing required columns: %s",
            paste(missing_snp_cols, collapse = ", ")
        ))
    }

    required_gene_cols <- c("chrom", "start", "end", "gene_name")
    missing_gene_cols <- setdiff(required_gene_cols, colnames(gene_anno))
    if (length(missing_gene_cols) > 0) {
        stop(sprintf(
            "gene_anno is missing required columns: %s",
            paste(missing_gene_cols, collapse = ", ")
        ))
    }

    if (nrow(snp_df) == 0) {
        snp_df$gene_name <- character(0)
        return(snp_df)
    }

    snp_tbl <- snp_df
    snp_tbl$gene_name <- NULL
    snp_tbl$.snp_row <- seq_len(nrow(snp_tbl))

    snps_gr <- plyranges::as_granges(
        snp_tbl,
        seqnames = chrom,
        start = pos,
        end = pos
    )

    gene_anno_gr <- plyranges::as_granges(
        gene_anno,
        seqnames = chrom,
        start = start,
        end = end
    )

    gene_summary <- plyranges::join_overlap_left(snps_gr, gene_anno_gr) %>%
        tibble::as_tibble() %>%
        dplyr::summarise(
            gene_name = {
                unique_genes <- unique(gene_name)
                unique_genes <- unique_genes[!is.na(unique_genes)]
                if (length(unique_genes) == 0) {
                    NA_character_
                } else {
                    paste(unique_genes, collapse = ", ")
                }
            },
            .by = .snp_row
        )

    out <- dplyr::left_join(snp_tbl, gene_summary, by = ".snp_row")
    out$.snp_row <- NULL
    out
}

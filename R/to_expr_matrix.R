#' Convert SNPData to expression-like matrix
#'
#' @param x SNPData object
#' @param level Aggregation level: "barcode", "clonotype", or "donor"
#' @return Matrix with snp_id as rows and cell/clonotype/donor as columns
#' @export
setGeneric("to_expr_matrix", function(x, level = c("barcode", "clonotype", "donor")) standardGeneric("to_expr_matrix"))

#' @rdname to_expr_matrix
#' @include SNPData-class.R
setMethod("to_expr_matrix", signature(x = "SNPData"), function(x, level = c("barcode", "clonotype", "donor")) {
    level <- match.arg(level)
    barcode_info <- get_barcode_info(x)

    if (level != "barcode" && !level %in% colnames(barcode_info)) {
        if (level == "clonotype") {
            stop("Clonotype information not available. Add clonotype data using add_barcode_metadata() or import_cellsnp() with vdj_file parameter.")
        } else {
            stop("No ", level, " column in barcode_info.")
        }
    }

    if (level == "barcode") {
        ref <- ref_count(x)
        alt <- alt_count(x)
        mat <- sign(ref - alt) * log1p(abs(ref - alt))
        rownames(mat) <- rownames(ref)
        colnames(mat) <- colnames(ref)
        return(mat)
    } else if (level == "clonotype") {
        clono <- barcode_info$clonotype
        # Check if all clonotype values are NA
        if (all(is.na(clono))) {
            stop("All clonotype values are NA. Cannot perform clonotype-level aggregation. Add clonotype data using add_barcode_metadata() or import_cellsnp() with vdj_file parameter.")
        }
        ref <- groupedRowSums(ref_count(x), clono)
        alt <- groupedRowSums(alt_count(x), clono)
        mat <- sign(ref - alt) * log1p(abs(ref - alt))
        rownames(mat) <- rownames(ref)
        colnames(mat) <- colnames(ref)
        return(mat)
    } else if (level == "donor") {
        donor <- barcode_info$donor
        ref <- groupedRowSums(ref_count(x), donor)
        alt <- groupedRowSums(alt_count(x), donor)
        mat <- sign(ref - alt) * log1p(abs(ref - alt))
        rownames(mat) <- rownames(ref)
        colnames(mat) <- colnames(ref)
        return(mat)
    }
})

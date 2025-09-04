#' Convert SNPData to expression-like matrix
#'
#' @param x SNPData object
#' @param level Aggregation level: "barcode", "clonotype", or "donor"
#' @return Matrix with snp_id as rows and cell/clonotype/donor as columns
#' @export
setGeneric("to_expr_matrix", function(x, level = c("barcode", "clonotype", "donor")) standardGeneric("to_expr_matrix"))

#' @rdname to_expr_matrix
setMethod("to_expr_matrix", signature(x = "SNPData"), function(x, level = c("barcode", "clonotype", "donor")) {
    level <- match.arg(level)
    barcode_info <- get_barcode_info(x)

    if (level != "barcode" && !level %in% colnames(barcode_info)) {
        stop("No ", level, " column in barcode_info.")
    }

    if (level == "barcode") {
        ref <- ref_count(x)
        alt <- alt_count(x)
        mat <- sign(ref - alt) * log1p(abs(ref - alt))
        rownames(mat) <- rownames(ref)
        colnames(mat) <- colnames(ref)
        return(mat)
    } else if (level == "clonotype") {
        clono <- get_barcode_info(x)$clonotype
        ref <- groupedRowSums(ref_count(x), clono)
        alt <- groupedRowSums(alt_count(x), clono)
        mat <- sign(ref - alt) * log1p(abs(ref - alt))
        rownames(mat) <- rownames(ref)
        colnames(mat) <- colnames(ref)
        return(mat)
    } else if (level == "donor") {
        donor <- get_barcode_info(x)$donor
        ref <- groupedRowSums(ref_count(x), donor)
        alt <- groupedRowSums(alt_count(x), donor)
        mat <- sign(ref - alt) * log1p(abs(ref - alt))
        rownames(mat) <- rownames(ref)
        colnames(mat) <- colnames(ref)
        return(mat)
    }
})

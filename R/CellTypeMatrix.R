#' Cell Type Identification Matrix
#' 
#' This function transforms the cell type information for each cell into an 
#' identification matrix that can be used in our downstream function `FRETCOM()`.
#'
#' @param stmeta A `data.frame` contains cell type information for each cell, where
#'              each row represents a cell. Please make sure that the order of rows 
#'              corresponds to the columns of the gene expression matrix.
#' @param celltype A variable name that represents the cell type information in `stmata`.
#'
#' @return A cell type identification matrix.
#' @export
#' @importFrom fastDummies dummy_cols
#' 
CellTypeMatrix <- function(stmeta, celltype) {
    
    if (!is.character(celltype)) {
        
        stop("Please provide the variable (column) name of the cell type in dataset!")
        
    }
    
    numcol <- ncol(stmeta)
    M <- as.matrix(dummy_cols(stmeta, select_columns=celltype)[, -c(1:numcol)])
    colnames(M) <- gsub(".*_", "", colnames(M))
    
    return(M)
    
}
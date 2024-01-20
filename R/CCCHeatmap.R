#' Cell-cell Communiction Heatmap
#'
#' @param restab A `data.frame` returned by `FRETCOM()`.
#' @param ligand A specific ligand name in `restab`, should be paired with `receptor`.
#' @param receptor A specific receptor name in `restab`, should be paired with `ligand`.
#' @param pathway A specific pathway name in `restab` (if any). Once `pathway` is 
#'              defined, `ligand` and `receptor` will not be considered.
#' @param celltype A vector including all cell type names that the user wants to 
#'                  include in the heatmap. If `celltype=NULL`, the cell types in 
#'                  `restab` will be used. The default is `NULL`.
#' @param bound A numeric value indicates the lower and upper bounds of the 
#'              coefficient estimates displayed in the heatmap. All values 
#'              greater than `abs(bound)` will be truncated to `abs(bound)`, 
#'              and all values smaller than `-abs(bound)` will be truncated 
#'              to `-abs(bound)`. The default is `1`.
#' @param cutoff A numeric value within the range of `(0, 1)` indicates the cutoff 
#'              for the adjusted p-values. If the adjusted p-value is less than 
#'              the specified cutoff, we treat it as significant, and it will be 
#'              displayed in the heatmap. The default is `0.2`.
#' @param ... Other arguments passed on to `ComplexHeatmap::draw()`.
#'
#' @return A heatmap produced by `ComplexHeatmap`.
#' 
#' @importFrom circlize colorRamp2
#' @importFrom grid unit gpar
#' @import ComplexHeatmap
#'
#' @export
CCCHeatmap <- function(restab, ligand=NULL, receptor=NULL, pathway=NULL, celltype=NULL, bound=1,
                       cutoff=0.2, ...) {
    
    if ((is.null(ligand) | is.null(receptor)) & is.null(pathway)) {
        
        stop("Please provide the ligand and receptor names or the signaling pathway name.")
        
    } else if (!is.null(pathway)) {
        
        if (pathway %in% restab$Pathway) {
            
            tempres <- restab[restab$Pathway == pathway, ]
            temptitle <- pathway
            
        } else {
            
            stop("The signaling pathway does not exist in restab! Please check the pathway name.")
            
        }
        
    } else {
        
        if (ligand %in% restab$Ligand & receptor %in% restab$Receptor) {
            
            tempres <- restab[restab$Ligand == ligand & restab$Receptor == receptor, ]
            temptitle <- paste0(ligand, " to ", receptor)
            
        } else {
            
            stop("The LR pair does not exist in restab! Please check the ligand and receptor names.")
            
        }
        
    }
    
    if (max(abs(tempres$Estimate)) > abs(bound)) {
        
        tempres$Estimate[tempres$Estimate < -abs(bound)] <- -abs(bound)
        tempres$Estimate[tempres$Estimate > abs(bound)] <- abs(bound)
        
    } else {
        
        bound <- max(abs(tempres$Estimate))
        
    }
    
    tempres$Qvalue[is.nan(tempres$Qvalue)] <- 1
    
    if (is.null(celltype)) {
        
        celltype <- unique(c(tempres$CellType_Sender, tempres$CellType_Receiver))
        
    }
    
    heatmat <- matrix(NA, nrow=length(celltype), ncol=length(celltype))
    rownames(heatmat) <- celltype
    colnames(heatmat) <- celltype
    
    for (i in 1:nrow(tempres)) {
        if (tempres$Qvalue[i] < cutoff) {
            heatmat[rownames(heatmat) == tempres$CellType_Sender[i],
                    colnames(heatmat) == tempres$CellType_Receiver[i]] <- tempres$Estimate[i]
        } else {
            heatmat[rownames(heatmat) == tempres$CellType_Sender[i],
                    colnames(heatmat) == tempres$CellType_Receiver[i]] <- 0
        }
    }
    
    selfcolor <- colorRamp2(c(-abs(bound), 0, abs(bound)), 
                            c("#86aaeb", "#FFFFFF", "#fc8383"))
    htmp <- Heatmap(heatmat, cluster_rows=FALSE, cluster_columns=FALSE, col=selfcolor,
                    row_title="Sender Cell Type", column_title="Receiver Cell Type",
                    heatmap_legend_param=list(title="Communication",
                                              title_position = "leftcenter-rot",
                                              labels_gp = gpar(fontsize = 5),
                                              title_gp = gpar(fontsize = 8, fontface = "bold"),
                                              legend_heigth = unit(1, "cm")),
                    rect_gp = gpar(col = "gray30", lwd = 1),
                    column_names_gp = gpar(fontsize = 8),
                    row_names_gp = gpar(fontsize = 8),
                    column_title_gp = gpar(fontsize = 10, fontface = "bold"),
                    row_title_gp = gpar(fontsize = 10, fontface = "bold"))
    draw(htmp, heatmap_legend_side="right", annotation_legend_side="right",
         column_title = temptitle,
         column_title_gp = gpar(fontsize = 12, fontface = "bold"),
         legend_grouping = "original", ...)
    
}

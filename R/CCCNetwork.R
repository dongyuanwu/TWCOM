#' Cell-cell Communiction Network
#'
#' @param restab A `data.frame` returned by `FRETCOM()`.
#' @param ligand A specific ligand name in `restab`, should be paired with `receptor`.
#' @param receptor A specific receptor name in `restab`, should be paired with `ligand`.
#' @param pathway A specific pathway name in `restab` (if any). Once `pathway` is 
#'              defined, `ligand` and `receptor` will not be considered.
#' @param cutoff A numeric value within the range of `(0, 1)` indicates the cutoff 
#'              for the adjusted p-values. If the adjusted p-value is less than 
#'              the specified cutoff, we treat it as significant, and it will be 
#'              displayed in the heatmap. The default is `0.2`.
#' @param weight A numeric value indicates the weight that determines the edge width. 
#'              Edge width is calculated as the product of the estimate and the weight.
#'              The default is `10`.
#' @param ... Other arguments passed on to `igraph::plot()`.
#'
#' @return A network produced by `igraph`.
#' 
#' @import igraph
#' @importFrom graphics title
#'
#' @export
CCCNetwork <- function(restab, ligand=NULL, receptor=NULL, pathway=NULL, cutoff=0.2, weight=10, ...) {
    
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
    
    tempres$Qvalue[is.nan(tempres$Qvalue)] <- 1

    tempres <- tempres[tempres$Qvalue < cutoff, ]
    
    netdf <- subset(tempres, select=c("CellType_Sender", "CellType_Receiver", "Estimate"))
    netdf$Estimate[netdf$Estimate < -1] <- -1
    netdf$Estimate[netdf$Estimate > 1] <- 1
    
    net <- graph_from_data_frame(netdf, directed=TRUE)
    E(net)$width <- abs(netdf$Estimate) * weight
    
    E(net)$color[netdf$Estimate < 0] <- '#a0bbeb'
    E(net)$color[netdf$Estimate > 0] <- '#FF9797'
    
    test.layout <- layout.reingold.tilford(net, circular=TRUE)
    
    plot(net, vertex.color="#C5E0B4", edge.color = E(net)$color,
         vertex.label.color="black", vertex.size=15,
         vertex.label.font=1.8, vertex.label.family="sans",
         vertex.frame.color="gray50", edge.color="gray", 
         vertex.label.cex=0.8, 
         edge.arrow.size=0.5, edge.curved=0.3,
         layout=test.layout, ...)
    
    title(temptitle, cex.main = 1.2)
    
}

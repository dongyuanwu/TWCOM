#' Frequentist Tweedie Modeling of Cell-cell Communication
#'
#' @param stdat A normalized gene expression matrix where each row represents a 
#'              gene and each column represents a cell/spot. Please make sure row
#'              names are gene names.
#' @param M A cell type proportion/identification matrix where each row represents
#'          a cell/spot and each column represents a cell type. Please ensure that 
#'          the row names (and order) match the column names (and order) of `stdat`.
#' @param ligands A list or vector of gene names for ligands corresponding to LR 
#'                pairs. Please make sure the order is the same as `receptors`.
#' @param receptors A list or vector of gene names for receptors corresponding to LR 
#'                pairs. Please make sure the order is the same as `ligands`.
#' @param ngridx An integer to indicate the number of grids that will be create 
#'              for the X coordinate if `re=TRUE`. Default is `5`.
#' @param ngridy An integer to indicate the number of grids that will be create 
#'              for the Y coordinate if `re=TRUE`. Default is `5`.
#' @param coordx A vector of X coordinates for cells/spots. Please make sure the 
#'              order can match the order of columns in `stdat`.
#' @param coordy A vector of Y coordinates for cells/spots. Please make sure the 
#'              order can match the order of columns in `stdat`.
#' @param maxdist The maximum distance allowed for communication. This number 
#'                should be based on the units of `coordx` and `coordy`. Usually, 
#'                we would recommend 200-300 μm.
#' @param sampleid A vector of sample IDs for cells/spots. Please make sure the 
#'              order can match the order of columns in `stdat`. Default is `NULL`,
#'              i.e., only considering one sample. 
#' @param rholist A vector that includes all possible rho values (communication 
#'              constraint parameter) that need to be tried in model fitting. 
#'              The best one will be automatically chosen based on AIC. Default
#'              value is `c(0.2, 0.5, 0.8)`.
#' @param parallel A logical value to indicate whether or not use parallel computing.
#'                  Default is `FALSE`.
#' @param ncores An integer to indicate the number of cores to be used for parallel 
#'              computing. This is only considered if `parallel=TRUE`. Default is `4`.
#'
#' @return A `data.frame` contains the inference results of cell-cell communication.
#' 
#' @references{
#'   \insertRef{wu2023inferring}{TWCOM}
#' }
#' 
#' @importFrom Rdpack reprompt
#' @import mgcv
#' @importFrom stats as.formula p.adjust quantile
#' @import parallel
#' @import doParallel
#' @import foreach
#' @import Rcpp
#' 
#' @export
#' 
FRETCOM <- function(stdat, M, ligands, receptors, ngridx=5, ngridy=5,
                    coordx, coordy, maxdist, sampleid=NULL, 
                    rholist=c(0.2, 0.5, 0.8),
                    parallel=FALSE, ncores=4) {
    
    celltype <- colnames(M)
    celltype_interaction <- paste0(rep(celltype, each=length(celltype)), 
                                   "_to_", rep(celltype, times=length(celltype)))
    
    ###### check ligand and receptor ######
    
    checkl <- sapply(ligands, function(x) any(x %in% rownames(stdat)))
    checkr <- sapply(receptors, function(x) any(x %in% rownames(stdat)))
    checklr <- checkl & checkr
    
    if (sum(!checklr) > 0) {
        
        ligands <- ligands[checklr]
        receptors <- receptors[checklr]
        message(sprintf("%s out of %s LR pairs have been removed due to their non-existence in SRT data.", 
                        sum(!checklr), length(checklr)))
        
    }
    
    ##
    
    if (is.null(sampleid)) {
        
        check_neighbor <- if_neighbors_samesample(coordx, coordy, dist=maxdist)
        
    } else {
        
        if (length(sampleid) != ncol(stdat)) {
            
            stop("The length of sampleid is different from the number of columns of stdat.")
            
        } else {
            
            check_neighbor <- if_neighbors(coordx, coordy, sampleid, dist=maxdist)
            
        }
        
    }
    
    ###### construct grids ######
    
    gridx <- cut(coordx, quantile(coordx, prob=(0:ngridx)/ngridx, names=FALSE, na.rm=TRUE), 
                 labels=FALSE, include.lowest=TRUE)
    gridy <- cut(coordy, quantile(coordy, prob=(0:ngridy)/ngridy, names=FALSE, na.rm=TRUE), 
                 labels=FALSE, include.lowest=TRUE)
    gridid <- as.numeric(as.factor(paste0(gridx, "_", gridy)))
    
    ###### find cell pairs that can communicate with each others ######
    
    spot <- data.frame(spoti=rep(1:ncol(stdat), each=ncol(stdat)), spotj=rep(1:ncol(stdat), times=ncol(stdat)))
    spot1 <- data.frame(spoti=rep(gridid, each=ncol(stdat)), spotj=rep(gridid, times=ncol(stdat)))
    
    spot <- as.matrix(spot[check_neighbor, ])
    spot1 <- as.matrix(spot1[check_neighbor, ])
    
    ##### calculate cell-to-cell communication scores Cij #####
    
    expr1 <- as.matrix(stdat)
    genenames <- rownames(stdat)
    
    Cspots <- vector(mode="list", length(ligands))
    
    for(k in 1:length(ligands)) {
        
        expr <- expr1[genenames %in% c(ligands[[k]], receptors[[k]]), ]
        
        if(class(expr)[1] != "matrix") {
            
            Cspots[[k]] <- get_Ck_v(expr, spot)
            
        } else {
            
            ligand <- which(rownames(expr) %in% ligands[[k]])
            receptor <- which(rownames(expr) %in% receptors[[k]])
            
            Cspots[[k]] <- get_Ck(expr, spot, ligand, receptor)
            
        }
        
    }
    
    D <- get_dist(coordx, coordy, spot)
    posi <- sapply(Cspots, function(x) sum(x) != 0)
    ligands <- ligands[posi]
    receptors <- receptors[posi]
    Cspots <- Cspots[posi]
    
    nonzero_prop <- sapply(Cspots, function(x) sum(x != 0) / nrow(spot))
    checkprop <- nonzero_prop >= 0.02
    
    if (sum(!checkprop) > 0) {
        
        ligands <- ligands[checkprop]
        receptors <- receptors[checkprop]
        Cspots <- Cspots[checkprop]
        message(sprintf("%s out of %s LR pairs have been removed due to a lack of gene expression information in SRT data.", 
                        sum(!checkprop), length(checkprop)))
        
    }
    
    ##### M * M #####
    
    if (!is.matrix(M)) {
        
        Ms <- get_Ms(as.matrix(M), spot)
        
    } else {
        
        Ms <- get_Ms(M, spot)
        
    }
    rmposi <- which(apply(Ms, 2, function(x) sum(x > 0)) < (nrow(Ms)*0.0001))
    if (length(rmposi) > 0) {
        
        if (length(rmposi) == ncol(Ms)) {
            
            stop("Not enough information for cell type interactions can be detectable.")
            
        } else {
            
            celltype_interaction <- celltype_interaction[-rmposi]
            Ms <- Ms[, -rmposi]
            
        }
        
    }
    
    ##### Modeling #####
    
    celltype_interaction_split <- strsplit(celltype_interaction, "_to_")
    
    if (length(rholist) == 1) {
        
        if (!parallel) {
            
            restab <- NULL
            
            for(k in 1:length(Cspots)) {
                
                message(sprintf("LR pairs: %s", k))
                
                Y <- Cspots[[k]]
                tempres <- gamfun_single(Ms=Ms, dist=D, Y=Y, spot=spot1, rho=rholist)
                
                temp_ligand <- paste0(ligands[[k]], collapse="_")
                temp_receptor <- paste0(receptors[[k]], collapse="_")
                
                temptab <- data.frame(Ligand=temp_ligand,
                                      Receptor=temp_receptor,
                                      CellType_Sender=sapply(celltype_interaction_split, function(x) x[1]),
                                      CellType_Receiver=sapply(celltype_interaction_split, function(x) x[2]),
                                      Estimate=tempres$res$coefest[-1],
                                      Pvalue=tempres$res$pval[-1],
                                      Qvalue=tempres$res$qval[-1])
                
                restab <- rbind(restab, temptab)
                
            }
            
        } else {
            
            n_cores <- min(ncores, detectCores())
            registerDoParallel(n_cores)
            
            restab <- foreach(k=1:length(Cspots), .combine=rbind,
                              .packages="mgcv", .export="gamfun_single") %dopar% {
                                  
                                  Y <- Cspots[[k]]
                                  tempres <- gamfun_single(Ms=Ms, dist=D, Y=Y, spot=spot1, rho=rholist)
                                  
                                  temp_ligand <- paste0(ligands[[k]], collapse="_")
                                  temp_receptor <- paste0(receptors[[k]], collapse="_")
                                  
                                  temptab <- data.frame(Ligand=temp_ligand,
                                                        Receptor=temp_receptor,
                                                        CellType_Sender=sapply(celltype_interaction_split, function(x) x[1]),
                                                        CellType_Receiver=sapply(celltype_interaction_split, function(x) x[2]),
                                                        Estimate=tempres$res$coefest[-1],
                                                        Pvalue=tempres$res$pval[-1],
                                                        Qvalue=tempres$res$qval[-1])
                                  
                                  temptab
                                  
                              }
            
        }
        
    } else {
        
        if (!parallel) {
            
            restab <- NULL
            
            for(k in 1:length(Cspots)) {
                
                message(sprintf("LR pairs: %s", k))
                
                Y <- Cspots[[k]]
                tempres <- gamfun(Ms=Ms, dist=D, Y=Y, spot=spot1, rholist=rholist)
                
                temp_ligand <- paste0(ligands[[k]], collapse="_")
                temp_receptor <- paste0(receptors[[k]], collapse="_")
                
                temptab <- data.frame(Ligand=temp_ligand,
                                      Receptor=temp_receptor,
                                      CellType_Sender=sapply(celltype_interaction_split, function(x) x[1]),
                                      CellType_Receiver=sapply(celltype_interaction_split, function(x) x[2]),
                                      Estimate=tempres$res$coefest[-1],
                                      Pvalue=tempres$res$pval[-1],
                                      Qvalue=tempres$res$qval[-1])
                
                restab <- rbind(restab, temptab)
                
            }
            
        } else {
            
            n_cores <- min(ncores, detectCores())
            registerDoParallel(n_cores)
            
            restab <- foreach(k=1:length(Cspots), .combine=rbind,
                              .packages="mgcv", .export="gamfun") %dopar% {
                                  
                                  message(sprintf("LR pairs: %s", k))
                                  
                                  Y <- Cspots[[k]]
                                  tempres <- gamfun(Ms=Ms, dist=D, Y=Y, spot=spot1, rholist=rholist)
                                  
                                  temp_ligand <- paste0(ligands[[k]], collapse="_")
                                  temp_receptor <- paste0(receptors[[k]], collapse="_")
                                  
                                  temptab <- data.frame(Ligand=temp_ligand,
                                                        Receptor=temp_receptor,
                                                        CellType_Sender=sapply(celltype_interaction_split, function(x) x[1]),
                                                        CellType_Receiver=sapply(celltype_interaction_split, function(x) x[2]),
                                                        Estimate=tempres$res$coefest[-1],
                                                        Pvalue=tempres$res$pval[-1],
                                                        Qvalue=tempres$res$qval[-1])
                                  
                                  temptab
                                  
                              }
            
        }
        
    }
    
    return(restab)
    
}



#' Frequentist Tweedie Modeling of Cell-cell Communication for Signaling Pathways
#'
#' @param stdat A normalized gene expression matrix where each row represents a 
#'              gene and each column represents a cell/spot. Please make sure row
#'              names are gene names.
#' @param M A cell type proportion/identification matrix where each row represents
#'          a cell/spot and each column represents a cell type. Please ensure that 
#'          the row names (and order) match the column names (and order) of `stdat`.
#' @param ligands A list or vector of gene names for ligands corresponding to LR 
#'                pairs. Please make sure the order is the same as `receptors`.
#' @param receptors A list or vector of gene names for receptors corresponding to LR 
#'                pairs. Please make sure the order is the same as `ligands`.
#' @param pathways A vector of signaling pathway names for corresponding LR 
#'                pairs. Please make sure the order is the same as `ligands` and `receptors`.
#' @param ngridx An integer to indicate the number of grids that will be create 
#'              for the X coordinate if `re=TRUE`. Default is `5`.
#' @param ngridy An integer to indicate the number of grids that will be create 
#'              for the Y coordinate if `re=TRUE`. Default is `5`.
#' @param coordx A vector of X coordinates for cells/spots. Please make sure the 
#'              order can match the order of columns in `stdat`.
#' @param coordy A vector of Y coordinates for cells/spots. Please make sure the 
#'              order can match the order of columns in `stdat`.
#' @param maxdist The maximum distance allowed for communication. This number 
#'                should be based on the units of `coordx` and `coordy`. Usually, 
#'                we would recommend 200-300 μm.
#' @param sampleid A vector of sample IDs for cells/spots. Please make sure the 
#'              order can match the order of columns in `stdat`. Default is `NULL`,
#'              i.e., only considering one sample. 
#' @param rholist A vector that includes all possible rho values (communication 
#'              constraint parameter) that need to be tried in model fitting. 
#'              The best one will be automatically chosen based on AIC. Default
#'              value is `c(0.2, 0.5, 0.8)`.
#' @param parallel A logical value to indicate whether or not use parallel computing.
#'                  Default is `FALSE`.
#' @param ncores An integer to indicate the number of cores to be used for parallel 
#'              computing. This is only considered if `parallel=TRUE`. Default is `4`.
#'
#' @return A `data.frame` contains the inference results of cell-cell communication.
#' 
#' @references{
#'   \insertRef{wu2023inferring}{TWCOM}
#' }
#' 
#' @importFrom Rdpack reprompt
#' @import mgcv
#' @importFrom stats as.formula p.adjust quantile
#' @import parallel
#' @import doParallel
#' @import foreach
#' @import Rcpp
#' 
#' @export
#' 
FRETCOMPathway <- function(stdat, M, ligands, receptors, pathways, 
                           ngridx=5, ngridy=5,
                           coordx, coordy, maxdist, sampleid=NULL, 
                           rholist=c(0.2, 0.5, 0.8),
                           parallel=FALSE, ncores=4) {
    
    celltype <- colnames(M)
    celltype_interaction <- paste0(rep(celltype, each=length(celltype)), 
                                   "_to_", rep(celltype, times=length(celltype)))
    
    ###### check ligand and receptor ######
    
    checkl <- sapply(ligands, function(x) any(x %in% rownames(stdat)))
    checkr <- sapply(receptors, function(x) any(x %in% rownames(stdat)))
    checklr <- checkl & checkr
    
    if (sum(!checklr) > 0) {
        
        ligands <- ligands[checklr]
        receptors <- receptors[checklr]
        pathways <- pathways[checklr]
        message(sprintf("%s out of %s LR pairs have been removed due to their non-existence in SRT data.", 
                        sum(!checklr), length(checklr)))
        
    }
    
    ##
    
    if (is.null(sampleid)) {
        
        check_neighbor <- if_neighbors_samesample(coordx, coordy, dist=maxdist)
        
    } else {
        
        if (length(sampleid) != ncol(stdat)) {
            
            stop("The length of sampleid is different from the number of columns of stdat.")
            
        } else {
            
            check_neighbor <- if_neighbors(coordx, coordy, sampleid, dist=maxdist)
            
        }
        
    }
    
    ###### construct grids ######
    
    gridx <- cut(coordx, quantile(coordx, prob=(0:ngridx)/ngridx, names=FALSE, na.rm=TRUE), 
                 labels=FALSE, include.lowest=TRUE)
    gridy <- cut(coordy, quantile(coordy, prob=(0:ngridy)/ngridy, names=FALSE, na.rm=TRUE), 
                 labels=FALSE, include.lowest=TRUE)
    gridid <- as.numeric(as.factor(paste0(gridx, "_", gridy)))
    
    ###### find cell pairs that can communicate with each others ######
    
    spot <- data.frame(spoti=rep(1:ncol(stdat), each=ncol(stdat)), spotj=rep(1:ncol(stdat), times=ncol(stdat)))
    spot1 <- data.frame(spoti=rep(gridid, each=ncol(stdat)), spotj=rep(gridid, times=ncol(stdat)))
    
    spot <- as.matrix(spot[check_neighbor, ])
    spot1 <- as.matrix(spot1[check_neighbor, ])
    
    ##### calculate cell-to-cell communication scores Cij #####
    
    expr1 <- as.matrix(stdat)
    genenames <- rownames(stdat)
    
    Cspots <- vector(mode="list", length(ligands))
    
    for(k in 1:length(ligands)) {
        
        expr <- expr1[genenames %in% c(ligands[[k]], receptors[[k]]), ]
        
        if(class(expr)[1] != "matrix") {
            
            Cspots[[k]] <- get_Ck_v(expr, spot)
            
        } else {
            
            ligand <- which(rownames(expr) %in% ligands[[k]])
            receptor <- which(rownames(expr) %in% receptors[[k]])
            
            Cspots[[k]] <- get_Ck(expr, spot, ligand, receptor)
            
        }
        
    }
    
    D <- get_dist(coordx, coordy, spot)
    posi <- sapply(Cspots, function(x) sum(x) != 0)
    # ligands <- ligands[posi]
    # receptors <- receptors[posi]
    Cspots <- Cspots[posi]
    pathways <- pathways[posi]
    
    pathways_uniq <- unique(pathways)
    Cspots_aggregate <- vector(mode="list", length=length(pathways_uniq))
    
    for(k in 1:length(pathways_uniq)) {
        
        Cspots_temp <- matrix(unlist(Cspots[pathways == pathways_uniq[k]]), ncol=length(posi))
        Cspots_aggregate[[k]] <- rowSums(Cspots_temp)
        
    }
    
    nonzero_prop <- sapply(Cspots_aggregate, function(x) sum(x != 0) / nrow(spot))
    checkprop <- nonzero_prop >= 0.02
    
    if (sum(!checkprop) > 0) {
        
        # ligands <- ligands[checkprop]
        # receptors <- receptors[checkprop]
        Cspots_aggregate <- Cspots_aggregate[checkprop]
        pathways <- pathways[checkprop]
        message(sprintf("%s out of %s signaling pathways have been removed due to a lack of gene expression information in SRT data.", 
                        sum(!checkprop), length(checkprop)))
        
    }
    
    if (length(Cspots_aggregate) > 1) {
        
        Cspots_temp <- matrix(unlist(Cspots_aggregate), ncol=length(posi))
        Cspots_aggregate <- c(Cspots_aggregate, rowSums(Cspots_temp))
        
    }
    
    ##### M * M #####
    
    if (!is.matrix(M)) {
        
        Ms <- get_Ms(as.matrix(M), spot)
        
    } else {
        
        Ms <- get_Ms(M, spot)
        
    }
    rmposi <- which(apply(Ms, 2, function(x) sum(x > 0)) < (nrow(Ms)*0.0001))
    if (length(rmposi) > 0) {
        
        if (length(rmposi) == ncol(Ms)) {
            
            stop("Not enough information for cell type interactions can be detectable.")
            
        } else {
            
            message(sprintf("%s out of %s cell type interactions have been removed due to a lack of information.", 
                            length(rmposi), length(celltype_interaction)))
            celltype_interaction <- celltype_interaction[-rmposi]
            Ms <- Ms[, -rmposi]
            
        }
        
    }
    
    ##### Modeling #####
    
    celltype_interaction_split <- strsplit(celltype_interaction, "_to_")
    
    if (length(rholist) == 1) {
        
        if (length(Cspots_aggregate) == 1) {
            
            message("Pathway: 1")
            
            Y <- Cspots_aggregate[[1]]
            tempres <- gamfun_single(Ms=Ms, dist=D, Y=Y, spot=spot1, rho=rholist)
            
            temp_pathway <- pathways_uniq
            
            restab <- data.frame(Pathway=temp_pathway,
                                  CellType_Sender=sapply(celltype_interaction_split, function(x) x[1]),
                                  CellType_Receiver=sapply(celltype_interaction_split, function(x) x[2]),
                                  Estimate=tempres$res$coefest[-1],
                                  Pvalue=tempres$res$pval[-1],
                                  Qvalue=tempres$res$qval[-1])
            
        } else {
            
            
            if (!parallel) {
                
                restab <- NULL
                
                for(k in 1:length(Cspots_aggregate)) {
                    
                    message(sprintf("Pathway: %s", k))
                    
                    Y <- Cspots_aggregate[[k]]
                    tempres <- gamfun_single(Ms=Ms, dist=D, Y=Y, spot=spot1, rho=rholist)
                    
                    temp_pathway <- ifelse(k == length(Cspots_aggregate), "Overall", pathways_uniq[k])
                    
                    temptab <- data.frame(Pathway=temp_pathway,
                                          CellType_Sender=sapply(celltype_interaction_split, function(x) x[1]),
                                          CellType_Receiver=sapply(celltype_interaction_split, function(x) x[2]),
                                          Estimate=tempres$res$coefest[-1],
                                          Pvalue=tempres$res$pval[-1],
                                          Qvalue=tempres$res$qval[-1])
                    
                    restab <- rbind(restab, temptab)
                    
                }
                
            } else {
                
                n_cores <- min(ncores, detectCores())
                registerDoParallel(n_cores)
                
                restab <- foreach(k=1:length(Cspots_aggregate), .combine=rbind,
                                  .packages="mgcv", .export="gamfun_single") %dopar% {
                                      
                                      message(sprintf("Pathway: %s", k))
                                      
                                      Y <- Cspots_aggregate[[k]]
                                      tempres <- gamfun_single(Ms=Ms, dist=D, Y=Y, spot=spot1, rho=rholist)
                                      
                                      temp_pathway <- ifelse(k == length(Cspots_aggregate), "Overall", pathways_uniq[k])
                                      
                                      temptab <- data.frame(Pathway=temp_pathway,
                                                            CellType_Sender=sapply(celltype_interaction_split, function(x) x[1]),
                                                            CellType_Receiver=sapply(celltype_interaction_split, function(x) x[2]),
                                                            Estimate=tempres$res$coefest[-1],
                                                            Pvalue=tempres$res$pval[-1],
                                                            Qvalue=tempres$res$qval[-1])
                                      
                                      temptab
                                      
                                  }
                
            }
            
        }
        
    } else {
        
        if (length(Cspots_aggregate) == 1) {
            
            message("Pathway: 1")
            
            Y <- Cspots_aggregate[[1]]
            tempres <- gamfun(Ms=Ms, dist=D, Y=Y, spot=spot1, rholist=rholist)
            
            temp_pathway <- pathways_uniq
            
            restab <- data.frame(Pathway=temp_pathway,
                                 CellType_Sender=sapply(celltype_interaction_split, function(x) x[1]),
                                 CellType_Receiver=sapply(celltype_interaction_split, function(x) x[2]),
                                 Estimate=tempres$res$coefest[-1],
                                 Pvalue=tempres$res$pval[-1],
                                 Qvalue=tempres$res$qval[-1])
            
        } else {
            
            if (!parallel) {
                
                restab <- NULL
                
                for(k in 1:length(Cspots_aggregate)) {
                    
                    message(sprintf("Pathway: %s", k))
                    
                    Y <- Cspots_aggregate[[k]]
                    tempres <- gamfun(Ms=Ms, dist=D, Y=Y, spot=spot1, rholist=rholist)
                    
                    temp_pathway <- ifelse(k == length(Cspots_aggregate), "Overall", pathways_uniq[k])
                    
                    temptab <- data.frame(Pathway=temp_pathway,
                                          CellType_Sender=sapply(celltype_interaction_split, function(x) x[1]),
                                          CellType_Receiver=sapply(celltype_interaction_split, function(x) x[2]),
                                          Estimate=tempres$res$coefest[-1],
                                          Pvalue=tempres$res$pval[-1],
                                          Qvalue=tempres$res$qval[-1])
                    
                    restab <- rbind(restab, temptab)
                    
                }
                
            } else {
                
                n_cores <- min(ncores, detectCores())
                registerDoParallel(n_cores)
                
                restab <- foreach(k=1:length(Cspots_aggregate), .combine=rbind,
                                  .packages="mgcv", .export="gamfun") %dopar% {
                                      
                                      message(sprintf("Pathway: %s", k))
                                      
                                      Y <- Cspots_aggregate[[k]]
                                      tempres <- gamfun_single(Ms=Ms, dist=D, Y=Y, spot=spot1, rho=rholist)
                                      
                                      temp_pathway <- ifelse(k == length(Cspots_aggregate), "Overall", pathways_uniq[k])
                                      
                                      temptab <- data.frame(Pathway=temp_pathway,
                                                            CellType_Sender=sapply(celltype_interaction_split, function(x) x[1]),
                                                            CellType_Receiver=sapply(celltype_interaction_split, function(x) x[2]),
                                                            Estimate=tempres$res$coefest[-1],
                                                            Pvalue=tempres$res$pval[-1],
                                                            Qvalue=tempres$res$qval[-1])
                                      
                                      temptab
                                      
                                  }
                
            }
            
        }
            
    }
    
    return(restab)
    
}



gamfun <- function(Ms, dist, Y, spot, rholist) {
    
    aic_best <- Inf
    rho_best <- 0
    
    tempdat <- as.data.frame(cbind(Y, Ms, spot))
    tempdat[, ncol(tempdat)-1] <- as.factor(tempdat[, ncol(tempdat)-1])
    tempdat[, ncol(tempdat)] <- as.factor(tempdat[, ncol(tempdat)])
    
    formu <- paste0(colnames(tempdat)[1], "~", 
                    paste0(colnames(tempdat)[-c(1, c(ncol(tempdat)-1, ncol(tempdat)))], 
                           collapse="+"))
    formu <- paste0(formu, "+s(", colnames(tempdat)[ncol(tempdat)-1],
                    ", bs='re') + s(", colnames(tempdat)[ncol(tempdat)], ", bs='re')")
    
    for (rho in rholist) {
        
        X <- scale(Ms * exp(-rho*dist), center=TRUE, scale=TRUE)
        tempdat[, 2:(ncol(tempdat)-2)] <- X
        
        tempaic <- NA
        tryCatch({
            tempaic <- bam(as.formula(formu), family=tw(), data=tempdat, method="ML")$aic
        }, error=function(e){})
        if(is.na(tempaic)) {
            tempaic <- gam(as.formula(formu), family=tw(), data=tempdat, method="ML")$aic
        }
        message(sprintf("Rho = %s: AIC = %s", rho, tempaic))
        if (tempaic < aic_best) {
            
            aic_best <- tempaic
            rho_best <- rho
            
        }
        
    }
    
    X <- scale(Ms * exp(-rho_best*dist), center=TRUE, scale=TRUE)
    tempdat[, 2:(ncol(tempdat)-2)] <- X
    
    gamfit <- tryCatch({summary(bam(as.formula(formu), family=tw(), data=tempdat))},
                       error=function(e){summary(gam(as.formula(formu), family=tw(), data=tempdat))})
    
    coefest <- as.vector(gamfit$p.coeff)
    pval <- as.vector(gamfit$p.pv)
    qval <- p.adjust(pval, method="fdr")
    
    return(list(res=data.frame(coefest=coefest, pval=pval, qval=qval),
                rho_best=rho_best))
    
}


gamfun_single <- function(Ms, dist, Y, spot, rho) {
    
    X <- scale(Ms * exp(-rho*dist), center=TRUE, scale=TRUE)
    tempdat <- as.data.frame(cbind(Y, X, spot))
    tempdat[, ncol(tempdat)-1] <- as.factor(tempdat[, ncol(tempdat)-1])
    tempdat[, ncol(tempdat)] <- as.factor(tempdat[, ncol(tempdat)])
    
    formu <- paste0(colnames(tempdat)[1], "~", 
                    paste0(colnames(tempdat)[-c(1, c(ncol(tempdat)-1, ncol(tempdat)))], 
                           collapse="+"))
    formu <- paste0(formu, "+s(", colnames(tempdat)[ncol(tempdat)-1],
                    ", bs='re') + s(", colnames(tempdat)[ncol(tempdat)], ", bs='re')")
    
    gamfit <- tryCatch({summary(bam(as.formula(formu), family=tw(), data=tempdat))},
                       error=function(e){summary(gam(as.formula(formu), family=tw(), data=tempdat))})
    
    coefest <- as.vector(gamfit$p.coeff)
    pval <- as.vector(gamfit$p.pv)
    qval <- p.adjust(pval, method="fdr")
    
    return(list(res=data.frame(coefest=coefest, pval=pval, qval=qval),
                rho_best=rho))
    
}


#' @importFrom mgcv ldTweedie
#' @export
mgcv::ldTweedie
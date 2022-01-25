SAClustering <- function(S.obj,NN_range=c(3,30), assay="RNA", slot="scale.data", 
                         clust_alg=1, type.fun="mean.silhouette",reduction=FALSE,
                         control=NULL, verbose=TRUE){
  
  ######
  require(GenSA)
  require(Seurat)
  
  cat("Function currently works with Seurat objects consider\n",
  "extending to more generic object types and generating the\n",
  "corresponding Seurat object from that.\n")
  
  cat("Comment and describe all of the parameters and allowed values.\n")
  # type.fun <- "mean.silhouette" # "mean.silhouette" "median.silhouette" "group.mean.silhouette" "group.median.silhouette"
  
  cat("Add other distance types rather than the sole correlation distance.\n",
      "Add possibility to personalize also the number of features that can be used within an assay.\n")
  
  # Process inputs to function 
  
  if (reduction==FALSE){# use original features 
    
    X <- switch(slot,
                "counts"={
                  X <- S.obj[[assay]]@counts
                },
                "data"={
                  X <- S.obj[[assay]]@data
                },
                "scale.data"={
                  X <- S.obj[[assay]]@data
                })
      
    
    cell.dims <- 2 # cells are along columns
    d <- 1 - cor(X)
    
    } else if (reduction==TRUE){
  
    # use principal components
    
    X <- S.obj@reductions$pca@cell.embeddings 
    num_PCs <- ncol(X)
    
    cell.dims <- 1 # cells are along rows 
    d <- 1 - cor(t(X))
  
    } else {
      
      stop("reduction must be logical.")
    
    }
  
  
  
  
  
  
  
  # Call to GenSA to solve the optimization problem
  
  
  fn.call <<-0
  
  par.history <<- matrix(c(0,0,0,0),nrow=1)
  colnames(par.history) <- c("res", "NN", "num.clusters","obj fun")
  
  
  out.SA <- GenSA(fn=obj.features,
                     par=NULL,
                     lower=lower,
                     upper=upper,
                     control = control,
                     d, S.obj,NN_range, assay.name, clust_alg, type.fun) # other parameters
  
  par.history <- par.history[-1,]
  
  out.SA$par[2] <- as.integer(floor(out.SA$par[2]*NN_range[2]))
  
  cat("Make also case for PCA-based obj.fn.\n Should you also add a require for Seurat and cluster and factoextra within obj?\n",
      "Clust alg for the moment is just Louvain, either remove it or consider adding many more.\n",
      "Also consider passing the parameter random.seed=some number in FindClusters.\n",
      "Also set control options, especially temperature and stopping conditions")
  
  cat("Optimization completed.\nCall functions", fn.calls, "times.\n")
  
  return(out.GenSA)
  
}




obj.features <- function(x,d,S.obj,NN_range, assay.name, clust_alg, type.fun){
  
  # first argument: x are the parameters SA is optimizing over
  # first element in x is resolution value, second element is num NN
  # additional aguments: additional parameters needed for computation
  # d = distance matrix
  # S.obj = Seurat object
  # min and max number of NNs
  # assay.name = perhaps needed for PCA-based clustering
  # clust_alg=Louvain, Leiden etc
  # type.fn=for the switch case statement
  # add random.seed=seed in the FindClusters
  
  fn.call <<- fn.call + 1
  
  NN <- as.integer(floor(x[2]*NN_range[2]))
  
  
  suppressWarnings(
    S.obj@graphs <- Seurat::FindNeighbors(d,
                                          distance.matrix = TRUE,
                                          verbose = FALSE,
                                          k.param = NN,
                                          annoy.metric = "euclidean",
                                          dims=NULL,
                                          reduction=NULL,
                                          #assay=assay.name,
                                          compute.SNN = TRUE)
  )
  
  names(S.obj@graphs) <- c("nn","snn")
  
  suppressWarnings(
    S.obj <- Seurat::FindClusters(S.obj,
                                  graph.name="snn",
                                  resolution=x[1],
                                  verbose=FALSE,
                                  modularity.fxn=1,
                                  algorithm=clust_alg)
  )
  
  num_clusts <- nlevels(S.obj$seurat_clusters)
  
  if (num_clusts == 1) {return(1)}
  
  s <- cluster::silhouette( as.integer(S.obj$seurat_clusters) , d)
  
  obj.fn <- obj.functions(sil=s,type.fun=type.fun)
  
  if (nlevels(S.obj$seurat_clusters) == 1){ obj.fn <- -1 }
  
  
  obj.fn <- -obj.fn # - (obj.fn) for optimization
  
  par.history <<- rbind(par.history, c(x[1],NN,num_clusts,obj.fn))
  
  cat(c(x[1],NN,num_clusts,obj.fn),"\n")
  
  
  return(obj.fn)
  
}


obj.reduction <- function(x,d,S.obj,NN_range, numPCs, assay.name, clust_alg, type.fun,optim.pc=FALSE){
  
  
  # describe inputs to all function
  fn.call <<- fn.call + 1
  
  NN <- as.integer(floor(x[2]*NN_range[2]))
  
  
  
  S.obj@graphs <- Seurat::FindNeighbors(S.obj,
                                        reduction="pca",
                                        verbose = TRUE,
                                        k.param = NN,
                                        annoy.metric = "euclidean",
                                        dims=numPCs,
                                        #assay=assay.name,
                                        compute.SNN = TRUE)
  
  names(S.obj@graphs) <- c("nn","snn")
  
  S.obj <- Seurat::FindClusters(S.obj,
                                graph.name="snn",
                                resolution=x[1],
                                verbose=TRUE,
                                modularity.fxn=1,
                                algorithm=clust_alg)
  
  
  
  
  
  
}








obj.functions <- function(sil,type.fun="mean.silhouette"){
  
  obj.fn <- switch(type.fun,
                   "mean.silhouette"={
                     obj.fn <- mean(sil[,"sil_width"])
                   },
                   "median.silhouette"={
                     obj.fn <- median(sil[,"sil_width"])
                   },
                   "group.mean.silhouette"={
                     obj.fn <- sapply( unique(sil[,"cluster"]),
                                       function(i) mean( sil[ sil[,1]==i ,"sil_width"] ) )
                     obj.fn <- mean(obj.fn)
                   },
                   "group.median.silhouette"={
                     obj.fn <- sapply( unique(s[,"cluster"]),
                                       function(i) median( sil[ sil[,1]==i ,"sil_width"] ) )
                     obj.fn <- mean(obj.fn)
                   })
  
  return(obj.fn)
  
  
}



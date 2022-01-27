#' Single Cell cluster determination by simulated annealing
#'
#' @description
#' 'SAClustering()' identifies clusters of cells by a shared nearest neighbor (SNN)
#' modularity optimization based clustering algorithm. `SAClustering()` makes use of
#' `FindNeighbors()` and `FindClusters()` from the Seurat package to construct a SNN
#' graph and identify clusters by optimizing a modularity function, respectively.
#' To help the user in the choice of a (close-to) optimal choice of the parameters,
#' number of nearest neighbors and resolution, `SAClustering()` optimizes a user-defined
#' objective function (see \code{type.fun}) using Simulated Annealing.
#'
#' @details
#'
#' @param S.obj a Seurat object
#' @param NN_range numeric vector containing the minimum and maximum number of nearest neighbors to
#' on which optimization should be performed.
#' @param assay assay
SAClustering <- function(S.obj,NN_range=c(3,30), assay="RNA", slot="scale.data",
                         clust_alg=1, type.fun="mean.silhouette",reduction=FALSE,
                         control=NULL, verbose=TRUE, final=TRUE,plot=TRUE)
  {

  ######SA
  #require(GenSA)
  #require(Seurat)

  cat("Function currently works with Seurat objects consider\n",
  "extending to more generic object types and generating the\n",
  "corresponding Seurat object from that.\n")

  cat("Comment and describe all of the parameters and allowed values.\n")
  # type.fun <- "mean.silhouette" # "mean.silhouette" "median.silhouette" "group.mean.silhouette" "group.median.silhouette"

  cat("Add other distance types rather than the sole correlation distance.\n",
      "Add possibility to personalize also the number of features that can be used within an assay.\n",
      "Add also description of value.\n")

  # Process inputs to function

  if (reduction==FALSE){# use original features

    X <- switch(slot,
                "counts"={
                  X <- as.matrix(S.obj[[assay]]@counts)
                },
                "data"={
                  X <- as.matrix(S.obj[[assay]]@data)
                },
                "scale.data"={
                  X <- S.obj[[assay]]@scale.data
                })


    cell.dims <- 2 # cells are along columns
    d <- 1 - stats::cor(X)


    } else if (reduction==TRUE){

    # use principal components

    X <- S.obj@reductions$pca@cell.embeddings
    num_PCs <- ncol(X)

    cell.dims <- 1 # cells are along rows
    d <- 1 - stats::cor(t(X))


    } else {

      stop("reduction must be logical.")

    }

    rm(X)




  # Call to GenSA to solve the optimization problem


  fn.call <<-0

  par.history <<- matrix(c(0,0,0,0),nrow=1)
  colnames(par.history) <- c("res", "NN", "num.clusters","obj fun")


  cat("Distinguish cases also for features-based and PCA-based when calling the function.\n")

  out.SA <- GenSA::GenSA(fn=obj.features,
                     par=NULL,
                     lower=lower,
                     upper=upper,
                     control = control,
                     d, S.obj,NN_range, assay.name, clust_alg, type.fun) # other parameters


  par.history <- par.history[-1,]
  out.SA$par[2] <- as.integer(floor(out.SA$par[2]*NN_range[2]))


  cat("Optimization completed.\nCall functions", fn.call, "times.\n")


  clustering.optimization <- list("out.SA",out.SA,
                                  "par.history", par.history)


  if (final == TRUE){ # compute clustering using parameters from the optimization routine

    x <- out.SA$par

    NN <- as.integer(floor(x[2]*NN_range[2]))


    S.obj@graphs <- Seurat::FindNeighbors(object=d,
                                          distance.matrix = TRUE,
                                          verbose = FALSE,
                                          k.param = NN,
                                          annoy.metric = "euclidean",
                                          #dims=NULL,
                                          #reduction=NULL,
                                          #assay=assay.name,
                                          compute.SNN = TRUE)


    names(S.obj@graphs) <- c("SA_nn","SA_snn")

    suppressWarnings(
      S.obj <- Seurat::FindClusters(object=S.obj,
                                graph.name="SA_snn",
                                resolution=x[1],
                                verbose=FALSE,
                                modularity.fxn=1,
                                algorithm=clust_alg)
    )


    clustering.optimization$Seurat_object <- S.obj

  }












  cat("Make also case for PCA-based obj.fn.\n Should you also add a require for Seurat and cluster and factoextra within obj?\n",
      "Clust alg for the moment is just Louvain, either remove it or consider adding many more.\n",
      "Also consider passing the parameter random.seed=some number in FindClusters.\n",
      "Also set control options, especially t`emperature and stopping conditions")






 return(clustering.optimization)



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


    S.obj@graphs <- Seurat::FindNeighbors(object=d,
                                          distance.matrix = TRUE,
                                          verbose = FALSE,
                                          k.param = NN,
                                          annoy.metric = "euclidean",
                                          #dims=NULL,
                                          #reduction=NULL,
                                          #assay=assay.name,
                                          compute.SNN = TRUE)

  names(S.obj@graphs) <- c("SA_nn","SA_snn")

  suppressWarnings(S.obj <- Seurat::FindClusters(object=S.obj,
                                  graph.name="SA_snn",
                                  resolution=x[1],
                                  verbose=FALSE,
                                  modularity.fxn=1,
                                  algorithm=clust_alg)
  )

  num_clusts <- nlevels(S.obj$seurat_clusters)

  if (num_clusts == 1) {return(1)}

  s <- cluster::silhouette( as.integer(S.obj$seurat_clusters), d)

  obj.fn <- obj.functions(sil=s,type.fun=type.fun)

  if (nlevels(S.obj$seurat_clusters) == 1){
    obj.fn <- -1
    }


  obj.fn <- -obj.fn # -(obj.fn) for optimization

  par.history <<- rbind(par.history, c(x[1],NN,num_clusts,obj.fn))

  cat(c(x[1],NN,num_clusts,obj.fn),"\n")


  return(obj.fn)

}


obj.reduction <- function(x,d,S.obj,NN_range, numPCs, assay.name, clust_alg, type.fun,optim.pc=FALSE){


  # describe inputs to all function
  fn.call <<- fn.call + 1

  NN <- as.integer(floor(x[2]*NN_range[2]))


  S.obj <- Seurat::FindNeighbors(object=S.obj,
                                        reduction="pca",
                                        verbose = TRUE,
                                        k.param = NN,
                                        annoy.metric = "euclidean",
                                        dims=numPCs,
                                        compute.SNN = TRUE)

  names(S.obj@graphs) <- c("SA_nn","SA_snn")

  S.obj <- Seurat::FindClusters(object=S.obj,
                                graph.name="snn",
                                resolution=x[1],
                                verbose=TRUE,
                                modularity.fxn=1,
                                algorithm=clust_alg)



  num_clusts <- nlevels(S.obj$seurat_clusters)

  if (num_clusts == 1) {return(1)}

  s <- cluster::silhouette( as.integer(S.obj$seurat_clusters), d)

  obj.fn <- obj.functions(sil=s,type.fun=type.fun)

  if (nlevels(S.obj$seurat_clusters) == 1){ obj.fn <- -1 }


  obj.fn <- -obj.fn # -(obj.fn) for optimization

  par.history <<- rbind(par.history, c(x[1],NN,num_clusts,obj.fn))

  cat(c(x[1],NN,num_clusts,obj.fn),"\n")


  return(obj.fn)


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



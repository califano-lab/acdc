#' Single Cell cluster determination by generalized simulated annealing
#'
#' @description
#' `SAClustering()` identifies clusters of cells by a shared nearest neighbor (SNN)
#' modularity optimization based clustering algorithm. `SAClustering()` makes use of
#' `FindNeighbors()` and `FindClusters()` from the Seurat package to construct a SNN
#' graph and identify clusters by optimizing a modularity function, respectively.
#' To help the user in the choice of a (close-to) optimal choice of the parameters,
#' number of nearest neighbors and resolution, `SAClustering()` optimizes a user-defined
#' objective function (see \code{type.fun}) using Simulated Annealing.
#'
#' @details
#' ASasssqqw=wqwq Continue with description of details
#'
#' @param S.obj A Seurat object
#' @param NN_range A numeric vector containing the minimum and maximum number of nearest neighbors to
#' on which optimization should be performed.
#' @param assay Assay to use in construction of (S)NN. Default is `"RNA"`; used only when `reduction` is `FALSE`
#' @param slot Slot to use in construction of (S)NN. Default is `scale.data`; used only when `reduction` is `FALSE`
#' @param clust_alg Algorithm for modularity optimization (input to `Seurat::FindClusters`). `1` = Louvain (default); `2` = Louvain
#' with multilevel refinement; `3` = SLM; `4` = Leiden (requires the leidenalg python). See `Seurat::FindClusters()`
#' for further details.
#' @param type.fun Objective function to be optimized by Simulated Annealing. Options include: `"mean.silhouette"` = mean
#' silhouette computed over all cells in the dataset (default); `"median.silhouette"` = median silhouette computed over all cells in
#' the dataset. `"group.mean.silhouette"` = mean of the per-group average silhouettes. `"group.median.silhouette"` = mean of the
#' per-group median silhouettes.
#' @param reduction Logical. Whether to perform silhouette-based optimized clustering using principal components (`TRUE`) or
#' original variables in the provided `slot` of the Seurat object `assay` (`FALSE`). Setting `reduction` = `TRUE` requires
#' a `DimReduc` object of name `"pca"` to be present in `S.obj`.
#' @param control List of options for Simulated Annealing. See Gubian el al. (2018) https://cran.r-project.org/web/packages/GenSA/GenSA.pdf for the
#' complete list of settings accepted. Default is `NULL`, i.e. default settings are employed in the
#' optimization of `type.fun`.
#' @param verbose Whether to print output. Default is `TRUE`.
#' @param final Whether `SAClustering()` should include a Seurat object with optimal clustering
#' results stored under `seurat_clusters` (thus overwritting pre-existent ones).
#' @param plot Whether to plot outcomes from clustering.
#'
#' @return Returns a list with the following fields:
#' \itemize{
#' \item Continue with decription of items
#' \item  weqew2
#' \item woqwepq3
#' }
#'
#' @examples Askosapa


SAClustering <- function(S.obj,NN_range=c(3,30), assay="RNA", slot="scale.data",
                         clust_alg=1, type.fun="mean.silhouette",reduction=FALSE,
                         control=NULL, verbose=TRUE, final=TRUE,plot=FALSE)
  {

  ######SA
  #require(GenSA)
  #require(Seurat)

  cat("Function currently works with Seurat objects consider\n",
  "extending to more generic object types and generating the\n",
  "corresponding Seurat object from that.\n")

  cat("Should I add a ... arguments passed to other methods?\n")

  cat("Comment and describe all of the parameters and allowed values.\n")

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

    if (IsMatrixEmpty(X)==TRUE) {

    }


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

  par.env <- new.env() # to store par values when parent and children functions communicate

  par.env$fn.call <-0
  assign("par.history", matrix(c(0,0,0,0),nrow=1), envir=par.env)


  cat("Distinguish cases also for features-based and PCA-based when calling the function.\n")

  cat("Optimizing ",type.fun," using generalized simulated annealing.\n" )

  out.SA <- GenSA::GenSA(fn=obj.features,
                     par=NULL,
                     lower=lower,
                     upper=upper,
                     control = control,
                     d, S.obj,NN_range, assay.name, clust_alg, type.fun, verbose, par.env) # other parameters


  par.env$par.history <- par.env$par.history[-1,]
  out.SA$par[2] <- as.integer(floor(out.SA$par[2]*NN_range[2]))
  colnames(par.env$par.history) <- c("res", "NN", "num.clusters","obj fun")


  cat("Optimization completed.\nCall functions", par.env$fn.call, "times.\n")


  clustering.optimization <- vector(mode="list")

  clustering.optimization$optim.par <- out.SA$par
  clustering.optimization$optim.value <- out.SA$value
  clustering.optimization$trace.mat <- out.SA$trace.mat
  clustering.optimization$counts <- out.SA$counts
  clustering.optimization$par.history <- par.env$par.history


  if (final == TRUE){ # compute clustering using parameters from the optimization routine

    x <- out.SA$par


    S.obj@graphs <- Seurat::FindNeighbors(object=d,
                                          distance.matrix = TRUE,
                                          verbose = verbose,
                                          k.param = x[2],
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
                                verbose=verbose,
                                modularity.fxn=1,
                                algorithm=clust_alg)
    )


    Seurat::Idents(S.obj) <- "seurat_clusters"



    # Displays silhouette plot

      s <- cluster::silhouette( as.integer(S.obj$seurat_clusters), d)

      #require(factoextra)
      plt.sil <- factoextra::fviz_silhouette(s)

      print(plt.sil)



    clustering.optimization$Seurat_object <- S.obj

  }













  cat("Make also case for PCA-based obj.fn.\n Should you also add a require for Seurat and cluster and factoextra within obj?\n",
      "Clust alg for the moment is just Louvain, either remove it or consider adding many more.\n",
      "Also consider passing the parameter random.seed=some number in FindClusters.\n",
      "Also set control options, especially t`emperature and stopping conditions")



 return(clustering.optimization)



}




obj.features <- function(x,d,S.obj,NN_range, assay.name, clust_alg, type.fun,verbose, par.env){

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

  par.env$fn.call <- par.env$fn.call + 1

  NN <- as.integer(floor(x[2]*NN_range[2]))


    S.obj@graphs <- Seurat::FindNeighbors(object=d,
                                          distance.matrix = TRUE,
                                          verbose = verbose,
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
                                  verbose=verbose,
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

  par.env$par.history <- rbind(par.env$par.history, c(x[1],NN,num_clusts,obj.fn))

  cat(c(x[1],NN,num_clusts,obj.fn),"\n")


  return(obj.fn)

}


obj.reduction <- function(x,d,S.obj,NN_range, numPCs, assay.name, clust_alg, type.fun, verbose,optim.pc=FALSE, par.env){


  # describe inputs to all function
  par.env$fn.call <- par.env$fn.call + 1

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

  par.env$par.history <- rbind(par.env$par.history, c(x[1],NN,num_clusts,obj.fn))

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



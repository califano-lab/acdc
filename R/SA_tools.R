#' Single Cell cluster determination by generalized simulated annealing
#'
#' @description
#' `SAClustering()` identifies clusters of cells by a shared nearest neighbor (SNN)
#' modularity optimization based clustering algorithm. `SAClustering()` makes use of
#' `FindNeighbors()` and `FindClusters()` from the Seurat package to construct a SNN
#' graph and identify clusters by optimizing a modularity function, respectively.
#' To help the user in the choice of a (close-to) optimal choice of the parameters,
#' number of nearest neighbors and resolution, `SAClustering()` optimizes a user-selected
#' objective function (see \code{type.fun}) using Simulated Annealing.
#'
#' @details
#' Continue with description of details
#'
#' @param S.obj A Seurat object
#' @param res.range A numeric vector containing the minimum and maximum values for resolution. Default is `c(0.01,2)`.
#' @param NN.range A numeric vector containing the minimum and maximum number of nearest neighbors to
#' on which optimization should be performed. Default is `c(3,30)`.
#' @param par.init A vector containing the initial values for the optimization parameters. First and second element
#' refer the the initial value of the resolution and number of nearest-neighbors, respectively. Values must be
#' within in the range specified by `res.range` and `NN.range`.
#' @param assay Assay to use in construction of (S)NN. Default is `"RNA"`; used only when `reduction` is `FALSE`
#' @param slot Slot to use in construction of (S)NN. Default is `scale.data`; used only when `reduction` is `FALSE`
#' @param reduction.slot reduction slot to use if `reduction` is set to `TRUE`, ignored otherwise. Default is `"pca"`
#' @param clust.alg Algorithm for modularity optimization (input to `Seurat::FindClusters`). `1` = Louvain (default); `2` = Louvain
#' with multilevel refinement; `3` = SLM; `4` = Leiden (requires the leidenalg python). See `Seurat::FindClusters()`
#' for further details.
#' @param type.fun Objective function to be optimized by Simulated Annealing. Options include: `"mean.silhouette"` = mean
#' silhouette computed over all cells in the dataset (default); `"median.silhouette"` = median silhouette computed over all cells in
#' the dataset. `"group.mean.silhouette"` = mean of the per-group average silhouettes. `"group.median.silhouette"` = mean of the
#' per-group median silhouettes.
#' @param reduction Logical. Whether to perform silhouette-based optimized clustering using principal components (`TRUE`) or
#' original variables in the provided `slot` of the Seurat object `assay` (`FALSE`). Setting `reduction` = `TRUE` requires
#' a `DimReduc` object of name `"pca"` to be present in `S.obj`.
#' @param control List of options for Simulated Annealing. Available options are: `maxit`, `threshold.stop`, `nb.stop.improvement`,
#' `smooth`, `max.call`, `max.time`, `temperature`, `visiting.param`, `acceptance.param`, `verbose`,
#' `simple.function`, `trace.mat`, `seed`. See Gubian el al. (2018) https://cran.r-project.org/web/packages/GenSA/GenSA.pdf for the
#' complete description of the settings. Default is `NULL`, i.e. default settings are employed in the
#' optimization of `type.fun`.
#' @param verbose Whether to print output of each function call. Default is `TRUE`.
#' @param diagnostics whether to print the outcomes of `FindNeighbors` and `FindClusters` at each function call. Default is `FALSE`.
#' @param final Whether `SAClustering()` should include a Seurat object with optimal clustering
#' results stored under `seurat_clusters` (thus overwritting pre-existent ones).
#' @param plot Whether to plot outcomes from clustering.
#' @param rng.seeds Seeds of the random number generators. The first element is used in `GenSA`, the second element is `FindClusters`.
#'
#' @return Returns a list with the following fields:
#' \itemize{
#' \item `optim.par` parameters corresponding to the optimal clustering solution obtained by generalized simulated annealing
#' \item  `optim.value` optimal value of the objective function
#' \item `trace.mat` matrix collecting the history of the algorithm, as produced by the GenSA package.
#' \item `num.evaluations` number of times the objective function is evaluated
#' \item `par.history` matrix collecting resolution, number of nearest neighbors, number of clusters and objective
#  function at each function call
#' \item `Seurat_object` copy of the input argument S.obj with optimal clustering solution stored in `seurat_clusters`.
#' Returned only if `final = TRUE`.
#' }
#'
#' @note Add notes
#'
#' @seealso \code{\link[GenSA]{GenSA}}, \code{\link[Seurat]{FindNeighbors}}, \code{\link[Seurat]{FindClusters}}
#'
#' @author Luca Zanella
#'
#' @examples
#' # Do not run
#' # Just to retrieve example data
#' # (devtools::install_github('satijalab/seurat-data') # if package SeuratData is needed, just for e.g.
#'
#' library(SeuratData) # just to retrieve some example data
#' AvailableData() # to see some example data
#' InstallData("pbmc3k")
#' pbmc3k.final <- LoadData("pbmc3k",type="pbmc3k.final")
#'
#' # Actual example
#' # Run SAClustering on gene expression data with default optimization parameters
#'
#' clust.optimization <- SAClustering(S.obj=pbmc3k.final,
#' res.range=c(0.1,1),
#' NN.range=c(3,15),
#' verbose=FALSE,
#' final=TRUE,
#' plot=TRUE)
#'
#'
#' # Run SAClustering using principal components as features
#' clust.optimization <- SAClustering(S.obj=pbmc3.final,
#'
#'@export



SAClustering <- function(S.obj,res.range=c(0.01,2),NN.range=c(3,30), par.init=NULL, assay="RNA", slot="scale.data", reduction=FALSE,
                         reduction.slot="pca",clust.alg=1, type.fun="mean.silhouette",
                         control=NULL, verbose=TRUE, final=TRUE, plot=FALSE, diagnostics=FALSE, rng.seeds=c(1234,0))
  {

  ######SA
  #require(GenSA)
  #require(Seurat)

  cat("Function currently works with Seurat objects consider\n",
  "extending to more generic object types.")

  cat("Should I add a ... arguments passed to other methods?\n")

  cat("Add other distance types rather than the sole correlation distance.\n",
      "Add possibility to personalize also the number of features that can be used within an assay.\n")

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

      stop("The provided slot is empty. Check combination of assay")

    }


    cell.dims <- 2 # cells are along columns
    d <- 1 - stats::cor(X)


    } else if (reduction==TRUE){

    # use principal components

    X <- S.obj@reductions[[reduction.slot]]@cell.embeddings
    numPCs <- ncol(X)

    cell.dims <- 1 # cells are along rows
    d <- 1 - stats::cor(t(X))


    } else {

      stop("reduction must be logical.")

    }

    rm(X)




  # Call to GenSA to solve the optimization problem

  lower <- c(res.range[1], NN.range[1]/NN.range[2])
  upper <- c(res.range[2], NN.range[2]/NN.range[2])

  if (!is.null(par.init)){
    par.init[2] <- par.init[2]/NN.range[2]
  }



  par.env <- new.env() # to store par values when parent and children functions communicate

  par.env$fn.call <-0
  assign("par.history", matrix(c(0,0,0,0),nrow=1), envir=par.env)


  set.seed(rng.seeds[1])
  cat("Optimizing ",type.fun," using generalized simulated annealing. Reduction set to ", as.character(reduction), ".\n" )

  if (reduction==FALSE) { # original features


    out.SA <- GenSA::GenSA(fn=obj.features,
                           par=par.init,
                           lower=lower,
                           upper=upper,
                           control = control,
                           d, S.obj,NN.range, assay.name, clust.alg, type.fun, verbose, diagnostics, rng.seeds, par.env) # other parameters

  } else if (reduction==TRUE) {

    out.SA <- GenSA::GenSA(fn=obj.reduction,
                           par=par.init,
                           lower=lower,
                           upper=upper,
                           control = control,
                           d,S.obj,NN.range, numPCs, assay.name, clust.alg, type.fun, verbose, diagnostics, optim.pc=FALSE, rng.seeds, par.env) # other parameters

  }




  par.env$par.history <- par.env$par.history[-1,]
  out.SA$par[2] <- as.integer(floor(out.SA$par[2]*NN.range[2]))
  colnames(par.env$par.history) <- c("res", "NN", "num.clusters","obj fun")


  cat("Optimization completed.\nCall functions", par.env$fn.call, "times.\n")


  clustering.optimization <- vector(mode="list")

  clustering.optimization$optim.par <- out.SA$par
  clustering.optimization$optim.value <- out.SA$value
  clustering.optimization$trace.mat <- out.SA$trace.mat
  clustering.optimization$num.evaluations <- out.SA$counts
  clustering.optimization$par.history <- par.env$par.history


  if (final == TRUE){ # compute clustering using parameters from the optimization routine

    x <- out.SA$par


    if (reduction==FALSE) { # original features

    S.obj@graphs <- Seurat::FindNeighbors(object=d,
                                          distance.matrix = TRUE,
                                          verbose = diagnostics,
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
                                verbose=diagnostics,
                                modularity.fxn=1,
                                random.seed=rng.seeds[2],
                                algorithm=clust.alg)
    )

    } else if (reduction==TRUE) { # principal components

      S.obj <- Seurat::FindNeighbors(object=S.obj,
                                     reduction="pca",
                                     verbose = diagnostics,
                                     k.param = x[2],
                                     annoy.metric = "euclidean",
                                     dims=1:numPCs,
                                     compute.SNN = TRUE)

      names(S.obj@graphs) <- c("SA_nn","SA_snn")

      S.obj <- Seurat::FindClusters(object=S.obj,
                                    graph.name="SA_snn",
                                    resolution=x[1],
                                    verbose=diagnostics,
                                    modularity.fxn=1,
                                    random.seed=rng.seeds[2],
                                    algorithm=clust.alg)

    }






    Seurat::Idents(S.obj) <- "seurat_clusters"

    # Displays silhouette plot

      s <- cluster::silhouette( as.integer(S.obj$seurat_clusters), d)

      #require(factoextra)
      plt.sil <- factoextra::fviz_silhouette(s)

      print(plt.sil)



    clustering.optimization$Seurat_object <- S.obj

  }












 return(clustering.optimization)



}




obj.features <- function(x,d,S.obj,NN.range, assay.name, clust.alg, type.fun,verbose, diagnostics, rng.seeds, par.env){

  # x are the parameters SA is optimizing over
  # first element in x is resolution value, second element is num NN
  # additional arguments: additional parameters needed for computation
  # d = distance matrix
  # S.obj = Seurat object
  # NN.range = min and max number of NNs
  # assay.name = perhaps needed for PCA-based clustering?
  # clust.alg=Louvain, Leiden etc
  # type.fn=for the switch case statement
  # add random.seed=seed in the FindClusters

  par.env$fn.call <- par.env$fn.call + 1

  NN <- as.integer(floor(x[2]*NN.range[2]))


    S.obj@graphs <- Seurat::FindNeighbors(object=d,
                                          distance.matrix = TRUE,
                                          verbose = diagnostics,
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
                                  verbose=diagnostics,
                                  modularity.fxn=1,
                                  random.seed=rng.seeds[2],
                                  algorithm=clust.alg)
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

  switch(verbose,"TRUE"={cat(c(x[1],NN,num_clusts,obj.fn),"\n")})


  return(obj.fn)

}


obj.reduction <- function(x,d,S.obj,NN.range, numPCs, assay.name, clust.alg, type.fun, verbose, diagnostics, optim.pc=FALSE, rng.seeds, par.env){


  # describe inputs to all function
  par.env$fn.call <- par.env$fn.call + 1

  NN <- as.integer(floor(x[2]*NN.range[2]))


  S.obj <- Seurat::FindNeighbors(object=S.obj,
                                        reduction="pca",
                                        verbose = diagnostics,
                                        k.param = NN,
                                        annoy.metric = "euclidean",
                                        dims=1:numPCs,
                                        compute.SNN = TRUE)

  names(S.obj@graphs) <- c("SA_nn","SA_snn")

  S.obj <- Seurat::FindClusters(object=S.obj,
                                graph.name="SA_snn",
                                resolution=x[1],
                                verbose=diagnostics,
                                modularity.fxn=1,
                                random.seed=rng.seeds[2],
                                algorithm=clust.alg)



  num_clusts <- nlevels(S.obj$seurat_clusters)

  if (num_clusts == 1) {return(1)}

  s <- cluster::silhouette( as.integer(S.obj$seurat_clusters), d)

  obj.fn <- obj.functions(sil=s,type.fun=type.fun)

  if (nlevels(S.obj$seurat_clusters) == 1){ obj.fn <- -1 }


  obj.fn <- -obj.fn # -(obj.fn) for optimization

  par.env$par.history <- rbind(par.env$par.history, c(x[1],NN,num_clusts,obj.fn))

  switch(verbose, "TRUE"={cat(c(x[1],NN,num_clusts,obj.fn),"\n")})



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



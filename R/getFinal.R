#' Single Cell cluster determination using (S)NN modularity optimization algorithm. Useful to retrieve the exactly the same
#' solution calculated using SAClustering, without running again the whole optimization routine.
#'
#'
#' @description
#' `getFinal()` returns clustering solution and the optimal score corresponding to the given input parameters,
#'  using a modularity optimization based clustering algorithm. `getFinal()` makes use of
#' `FindNeighbors()` and `FindClusters()` from the Seurat package to construct a SNN
#' graph and identify clusters by optimizing a modularity function. To evaluate the quality of the returned solution
#' `getFinal()` returns the value associated to a user-defined metric (see \code{type.fun}).
#'
#' @details
#' Continue with description of details
#'
#' @param S.obj A Seurat object
#' @param res value for resolution (numeric). Default is `0.5`.
#' @param NN number of nearest neighbors (numeric). Default is `15`.
#' @param assay assay to use in construction of (S)NN. Default is `"RNA"`; used only when `reduction` is `FALSE`
#' @param slot slot to use in construction of (S)NN. Default is `scale.data`; used only when `reduction` is `FALSE`
#' @param reduction.slot reduction slot to use if `reduction` is set to `TRUE`, ignored otherwise. Default is `"pca"`
#' @param num.pcs number of principal components to use for the construction of the (S)NN (numeric). Used only if `reduction = TRUE`. Default is `NULL`,
#' meaning that all principal components in the dimensionality reduction slot are employed in (S)NN construction.
#' @param clust.alg Algorithm for modularity optimization (input to `Seurat::FindClusters`). `1` = Louvain (default); `2` = Louvain
#' with multilevel refinement; `3` = SLM; `4` = Leiden (requires the leidenalg python). See `Seurat::FindClusters()`
#' for further details.
#' @param type.fun Metric to evaluate the quality of the clustering solution. Options include: `"mean.silhouette"` = mean
#' silhouette computed over all cells in the dataset (default); `"median.silhouette"` = median silhouette computed over all cells in
#' the dataset. `"group.mean.silhouette"` = mean of the per-group average silhouettes. `"group.median.silhouette"` = mean of the
#' per-group median silhouettes. `"generalized.logistic"` = mean of the transformed-silhouette computed using a generalized logistic
#' transformation. See vignette for further details.
#' @param weights weights assigned to negative silhouette scores in the calculation of the objective function `type.fun`. Possible values are
#' either `"unitary"`, i.e. negative silhouette scores are used in the calculation of the objective function as they are, or `"exp"`, i.e.
#' negative silhouette scores are used in the calculation of the objective function after exponentiation. Default is `"unitary"`.
#' @param reduction Logical. Whether to perform clustering using principal components (`TRUE`) or
#' original variables in the provided `slot` of the Seurat object `assay` (`FALSE`). Setting `reduction` = `TRUE` requires
#' a `DimReduc` object of name `"pca"` to be present in `S.obj`.
#' @param verbose whether to print the outcomes of `FindNeighbors` and `FindClusters` at each function call. Default is `FALSE`.
#' @param rng.seed Seed of the random number generator used in `FindClusters`.
#'
#' @return Returns an object of class Seurat with the with optimal clustering solution stored in the metadata `seurat_clusters`, the
#' corresponding `silhouette` object stored in `Seurat_object[[assay]]@misc$sil` and the value of the metric `type.fun` for the assessment of the cluster quality
#' `Seurat_object[[assay]]@misc$metric`.
#'
#' @note Add notes
#'
#' @seealso \code{\link[Seurat]{FindNeighbors}}, \code{\link[Seurat]{FindClusters}}
#'
#' @author Luca Zanella
#'
#' @examples
#'
#'
#'
#' \dontrun{
#' # Just to retrieve example data
#' # devtools::install_github('satijalab/seurat-data') # if package SeuratData is needed, just for e.g.
#'
#' library(SeuratData) # just to retrieve some example data
#' AvailableData() # to see some example data
#' InstallData("pbmc3k")
#' pbmc3k.final <- LoadData("pbmc3k",type="pbmc3k.final")
#'
#' Actual example
#' Get clustering solution on gene expression data with input parameters
#' Get clustering solution using principal components as features and add the
#' clustering solution to the metadata under the voice `seurat_clusters` and the
#' median silhouette computed across all cells as the output metric in `S.obj[["RNA"]]@misc$metric`.
#'
#' S.obj <- getFinal(S.obj=pbmc3k.final,
#' res=1,
#' NN=30,
#' reduction=TRUE,
#' type.fun="median.silhouette")
#' }
#'
#' \dontrun{
#' The following example uses reduction = FALSE.
#' clustering.output <- getFinal(S.obj=pbmc3k.final,
#' res=0.5,
#' NN=15,
#' reduction=FALSE,
#' verbose=TRUE)
#'}
#'
#'
#'
#'
#'@export


getFinal <- function(
    S.obj,
    res = 0.5,
    NN = 15,
    assay = "RNA",
    slot = "scale.data",
    reduction = TRUE,
    reduction.slot = "pca",
    num.pcs = NULL,
    verbose = FALSE,
    clust.alg = 1,
    type.fun = "mean.silhouette",
    weights = "unitary",
    exp_base = 2.7182,
    rng.seed = 0,
    object.dist = NULL
)

{

  #require(Seurat)
  require(dplyr)
  # Process inputs to function

  if(is.null(object.dist)){
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
      d <- sqrt(1 - stats::cor(X))


    } else if (reduction==TRUE){

      # use principal components

      X <- S.obj@reductions[[reduction.slot]]@cell.embeddings
      numPCs <- ncol(X)
      cell.dims <- 1 # cells are along rows
      d <- sqrt(1 - stats::cor(t(X)))

    } else {

      stop("reduction must be logical.")

    }

    rm(X)
  } else {
    d <- object.dist
    if(reduction==FALSE){
      cell.dims <- 2 # cells are along columns
    } else if(reduction==TRUE) {
      numPCs <- ncol(S.obj@reductions[[reduction.slot]]@cell.embeddings)
      cell.dims <- 1 # cells are along rows
    }
  }




  ######
  # Compute solution with optimal clustering parameters and return Seurat object

  if (reduction==FALSE) { # original features

    S.obj@graphs <- Seurat::FindNeighbors(object=d,
                                          distance.matrix = TRUE,
                                          verbose = verbose,
                                          k.param = NN,
                                          annoy.metric = "euclidean",
                                          #dims=NULL,
                                          #reduction=NULL,
                                          #assay=assay.name,
                                          compute.SNN = TRUE)


    names(S.obj@graphs) <- c("opt_nn","opt_snn")

    suppressWarnings(
      S.obj <- Seurat::FindClusters(object=S.obj,
                                    graph.name="opt_snn",
                                    resolution=res,
                                    verbose=verbose,
                                    modularity.fxn=1,
                                    random.seed=rng.seed,
                                    algorithm=clust.alg)
    )

  } else if (reduction==TRUE) { # principal components

    if (is.null(num.pcs)) {
      S.obj <- Seurat::FindNeighbors(object=S.obj,
                                     reduction=reduction.slot,
                                     verbose = verbose,
                                     k.param = NN,
                                     annoy.metric = "euclidean",
                                     dims=1:numPCs,
                                     compute.SNN = TRUE)

    } else if (is.null(num.pcs) == FALSE) {

      if (num.pcs>numPCs){
        stop("The provided number of principal components (num.pcs) is greater than the number of principal component in the reduction slot.\n")
      }

      S.obj <- Seurat::FindNeighbors(object=S.obj,
                                     reduction=reduction.slot,
                                     verbose = verbose,
                                     k.param = NN,
                                     annoy.metric = "euclidean",
                                     dims=1:num.pcs,
                                     compute.SNN = TRUE)

    }

    names(S.obj@graphs) <- c("opt_nn","opt_snn")

    S.obj <- Seurat::FindClusters(object=S.obj,
                                  graph.name="opt_snn",
                                  resolution=res,
                                  verbose=verbose,
                                  modularity.fxn=1,
                                  random.seed=rng.seed,
                                  algorithm=clust.alg)


  }




  Seurat::Idents(S.obj) <- "seurat_clusters"

  # Displays silhouette plot

  if ( nlevels(S.obj$seurat_clusters) > 1 ) {

    s <- cluster::silhouette( as.integer(S.obj$seurat_clusters), d)


    # sil_neg <- sapply( unique(s[,"cluster"]),
    #                    function(i) { sum( s[s[,1]==i, "sil_width"] < lq ) / nrow( s[s[,1]==i,] ) } )
    #


    S.obj[[assay]]@misc$sil <- s






    #require(factoextra)
    #require(dplyr)

    # plt.sil <- factoextra::fviz_silhouette(s,print.summary = FALSE)

    # switch(verbose, "TRUE"={print(plt.sil)})

    # plt.sil <- plt.sil$data %>%
    #   group_by(cluster) %>%
    #   summarise(size = n(),
    #             ave.sil.width=round(mean(sil_width), 2)) %>%
    #   as.data.frame()








    # Return metric for the given run
    metric <- obj.functions(S.obj = S.obj,
                            d = d,
                            assay.name = assay,
                            slot = slot,
                            type.fun=type.fun,
                            weights=weights,
                            exp_base=exp_base)
    names(metric) <- type.fun

    S.obj[[assay]]@misc$metric <- metric

  } else {
    S.obj[[assay]]@misc$sil <- NA
    S.obj[[assay]]@misc$metric <- NA
  }

  return(S.obj)



}

#' GridSearch2
#' Computes the optimal combination of k and resolution for
#' clustering, by grid-search optimization of the Silhouette score across multiple
#' resamplings
#'
#' @param object Seurat object (genes/proteins along rows and cells on columns)
#' @param assay.name Assay to use
#' @param .resolutions vector of the resolution parameter. Default=seq(0.01,2,by = 0.01).
#' @param .bootstraps vector of bootstraps. Default=1 (i.e. no bootstraps).
#' @param .knns sequence of values for the number of nearest neighbors. Default=seq(3, 31, by=2).
#' @param .pct_cells percentage of cells sample at each bootstrap iteration. Default is 100% (all cells are used).
#' @param .replace: (logical) whether to sample cells with (TRUE) or without replacement (FALSE). Default=FALSE.
#' @param .clust_alg clustering algorithm. Choose among: "Louvain" (default); "Louvain-mult-ref"; "SLM"; "Leiden".
#' @param free_cores number of cores that are not used for the calculation. Default=2.
#' @param my_seed random seed for FindClusters. Default=0.
#' @return A 11-columns tibble containing the outcomes of the calculation for choosing the optimal number of clusters. \cr
#'
#' Columns consist of the following vectors.
#' \code{index}: integer in the interval 1-\emph{n_idx}, where \emph{n_idx}=\code{length(.bootstrap)*}\code{length(.knns)}
#' uniquely assigned to each combination of \code{.bootstraps} and \code{.knns}. \cr
#' \code{boostrap}: integer assigned to the given resampling. Possible values are those in \code{.bootstraps}. \cr
#' \code{knn}: number of nearest neighbors used as input parameters to FindNeighbors. Possible values are those in \code{.knns}. \cr
#' \code{resolution}: resolution for FindClusters. Possible values are those in \code{.resolutions}. \cr
#' \code{tot_sil_neg}: \cr
#' \code{lowest_sil_clust}: \cr
#' \code{sil_avg}: \cr
#' \code{sil_mean_median}: \cr
#' \code{n_clust}: number of clusters resulting for the given combination of \code{knn}, \code{resolution} with cells sampled according to \code{bootstrap}. \cr
#' \code{random.seed}: random seed used to initialize the rng for cell subsamplings. Values are the same as \code{index}.
#'
#' @export
GridSearch2 <- function(object,
                       assay.name,
                       .resolutions=seq(0.01,2,by = 0.01),
                       .bootstraps=1,
                       .knns=seq(3,31,by = 2),
                       .pct_cells=100,
                       .replace=FALSE,
                       .clust_alg="Louvain",
                       free.cores=2,
                       my_seed=0)
{

  ## Grid Search ----
  print("Selection of parameters for optimal clustering solution")

  {


    if (Sys.info()['sysname'] == "Windows"){
      clust.type <- "PSOCK"
    } else if ( (Sys.info()['sysname'] == "Linux") | (Sys.info()['sysname'] == "Darwin") ) {
      clust.type <- "FORK"
    }

    n_cores <- parallel::detectCores()-free.cores
    myCluster <- parallel::makeCluster(n_cores,type = clust.type)
    doParallel::registerDoParallel(myCluster)
    n_cells_to_subsample <- round(ncol(object)*.pct_cells/100)
    n_cells_to_subsample


    .clust_alg = switch(.clust_alg,
                        "Louvain"= 1,
                        "Louvain-mult-ref"=2,
                        "SLM"=3,
                        "Leiden"=4)


     print(system.time({
    result <- foreach::foreach( a_knn = .knns, .combine = 'rbind' ) %dopar% {
        resolution = rep( .resolutions, length(.bootstraps) )
        bootstrap = rep( .bootstraps, length(.resolutions) )
        index <- 1:length(resolution)
        my_sil.df <- dplyr::tibble( index = index,
                                    bootstrap = bootstrap, knn = 0, resolution = resolution,
                                    tot_sil_neg = 0, lowest_sil_clust = 0, max_sil_clust = 0,
                                    sil_avg = 0, sil_mean_median = 0, n_clust = 0, random.seed = 0 )
        nrow(my_sil.df)
        selected_samples <- sample(colnames(object),size=n_cells_to_subsample,replace=.replace)
        seurat_obj <- object[ , colnames(object) %in% selected_samples ]
        # x.dist <- as.matrix(vpdist)[selected_samples,selected_samples]
        # x.dist <- as.matrix(as.dist( 1-cor( x@assays$VIPER@scale.data ,
        #                                     method = “spe” )))
        # x.dist <- as.dist(viperSimilarity(x@assays$VIPER@scale.data)) %>% as.matrix()
        x.dist <- as.matrix(as.dist( 1-cor( seurat_obj[[assay.name]]@scale.data,method = "pea" )))
        my_ssn_graph <- Seurat::FindNeighbors( x.dist,
                                               # dims = pcs_to_use,
                                               # assay = “VIPER”,
                                               distance.matrix = TRUE,
                                               verbose = TRUE,
                                               # k.param = knn_n_neighbor,
                                               k.param = a_knn,
                                               # annoy.metric = “cosine”,
                                               annoy.metric = "euclidean",
                                               compute.SNN = TRUE )
        seurat_obj@graphs <- my_ssn_graph
        for ( a_index in my_sil.df$index )
        {
          set.seed(a_index)
          my_sil.df$random.seed[a_index] <- a_index
          print(paste0("--- Silhouette score computation knn/bootstrap : ", a_knn , "|" , my_sil.df$bootstrap[a_index] ))
          seurat_obj <- Seurat::FindClusters( seurat_obj,
                                     graph.name = "snn",
                                     resolution = my_sil.df$resolution[a_index],
                                     verbose = FALSE,
                                     modularity.fxn = 1,
                                     algorithm = .clust_alg,
                                     random.seed = my_seed
          )
          table(seurat_obj$seurat_clusters)
          if ( nlevels(seurat_obj$seurat_clusters) == 1 ) next ;
          s <- cluster::silhouette( as.integer(seurat_obj$seurat_clusters) , x.dist )
          # pdf( file.path(reports.dir,paste0(“sil-res-“,a_res,“.pdf”)) )
          x <- factoextra::fviz_silhouette(s,print.summary = FALSE)
          y <- sapply( levels(x$data$cluster) , function(i) mean( x$data$sil_width[ x$data$cluster == i ] ) )
          z <- sapply( levels(x$data$cluster) , function(i) median( x$data$sil_width[ x$data$cluster == i ] ) )
          tot_sil_neg <- sapply( levels(x$data$cluster) , function(i) sum( x$data$sil_width[ x$data$cluster == i ] < 0.25 ) )
          my_sil.df$sil_avg[ a_index ] = mean(y)
          my_sil.df$sil_mean_median[ a_index ] = mean(z)
          my_sil.df$tot_sil_neg[ a_index ] = sum(tot_sil_neg)
          my_sil.df$lowest_sil_clust[ a_index ] = min(y)
          my_sil.df$max_sil_clust[a_index] = max(y)
          my_sil.df$n_clust[ a_index ] = nlevels(x$data$cluster)
          # dev.off()
          # View(my_sil.df)
        }
        my_sil.df$knn = a_knn
        return(my_sil.df)
      } # End of dopar
    })) # End of print
    stopCluster(myCluster)
    my_sil.df <- result

  }


}

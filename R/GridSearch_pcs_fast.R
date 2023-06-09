#' Determination of the optimal combination of k and resolution for
#' single-cell clustering.
#'
#' @description
#' `GridSearch_pcs_fast()` optimizes clustering solution of scRNA-seq databy
#' grid-search optimization of the Silhouette score across multiple
#' resamplings. Runs faster by computing the distance matrix once if not bootstrapping
#' and placing NN in the outer loop and resolution in the inner loop, reducing the number
#' of times that neighbors are computed.
#'
#' @details
#' Input arguments
#'
#' @param object Seurat object (genes/proteins along rows and cells on columns)
#' @param assay.name Assay to use
#' @param .resolutions vector of the resolution parameter. Default=seq(0.01,2,by = 0.01).
#' @param .bootstraps vector of bootstraps. Default=1 (i.e. no bootstraps).
#' @param .knns sequence of values for the number of nearest neighbors. Default=seq(3, 31, by=2).
#' @param .pct_cells percentage of cells sample at each bootstrap iteration. All cells are used by default.
#' @param .replace: (logical) whether to sample cells with (TRUE) or without replacement (FALSE). Default=FALSE.
#' @param .clust_alg clustering algorithm. Choose among: "Louvain" (default); "Louvain-mult-ref"; "SLM"; "Leiden".
#' @param type.fun Objective function to be optimized when computing optimal solution. Options include: `"mean.silhouette"` = mean
#' silhouette computed over all cells in the dataset (default); `"median.silhouette"` = median silhouette computed over all cells in
#' the dataset. `"group.mean.silhouette"` = mean of the per-group average silhouettes. `"group.median.silhouette"` = mean of the
#' per-group median silhouettes. `"generalized.logistic"` = mean of the transformed-silhouette computed using a generalized logistic
#' @param free_cores number of cores that are not used for the calculation. Default=2.
#' @param my_seed random seed for FindClusters. Default=0.
#' @param weights unitary (unitary) or exponential (exp) way of weighing the silhouette scores. Default=unitary.
#'
#' @return Returns an object of class Seurat with the with optimal clustering solution stored in the metadata `seurat_clusters`, the
#' corresponding `silhouette` object stored in `Seurat_object[[assay]]@misc$sil` and a tibble containing  the history of the optimization
#' algorithm stored in `Seurat_object[[assay]]@misc$GS.history`. The 11-columns tibble containing the outcomes of the calculation for choosing the optimal number of clusters.
#' Columns consist of the following vectors.
#' \itemize{
#' \item `index`: integer in the interval 1-`n_idx`, where `n_idx`=`length(.bootstrap)*length(.knns)`
#' uniquely assigned to each combination of `.bootstraps` and `.knns`.
#' \item `boostrap`: integer assigned to the given resampling. Possible values are those in `.bootstraps`.
#' \item `knn`: number of nearest neighbors used as input parameters to `FindNeighbors`. Possible values are those in `.knns`.
#' \item `resolution`: resolution for `FindClusters`. Possible values are those in `.resolutions`.
#' \item `tot_sil_neg`: negative silhouette score
#' \item `lowest_sil_clust`: lowest silhouette
#' \item `sil_avg`: average silhouette score
#' \item `sil_mean_median`: mean of the median silhouette score
#' \item `n_clust`: number of clusters resulting for the given combination of `knn`, `resolution` with cells sampled according to `bootstrap`.
#' \item `random.seed`: random seed used to initialize the rng for cell subsamplings. Values are the same as `index`.
#' }
#'
#'
#' @export



GridSearch_pcs_fast <- function(object,
                                assay.name,
                                .resolutions = seq(0.1, 1.9, by = 0.2),
                                .bootstraps = 1,
                                .knns = seq(11, 101, by = 10),
                                .pct_cells = 100,
                                .replace = FALSE,
                                .clust_alg = "Louvain",
                                type.fun = "group.mean.silhouette",
                                free.cores = 2,
                                .type = "genes",
                                .dims = NULL,
                                my_seed = 0,
                                weights = "unitary",
                                exp_base = 2.718282,
                                show_progress_bar = TRUE)

{

  cat("This is a beta version. acdc is currently under development.\n")

  suppressMessages(require(foreach))
  ## Silhouette Analysis ----
  print("Selection of parameters for optimal clustering solution")

  {


    if (Sys.info()['sysname'] == "Windows"){
      clust.type <- "PSOCK"
    } else if ( (Sys.info()['sysname'] == "Linux") | (Sys.info()['sysname'] == "Darwin") ) {
      clust.type <- "FORK"
    }

    n_cores <- parallel::detectCores()-free.cores
    n_cells_to_subsample <- round(ncol(object)*.pct_cells/100)
    n_cells_to_subsample


    .clust_alg = switch(.clust_alg,
                        "Louvain"= 1,
                        "Louvain-mult-ref"=2,
                        "SLM"=3,
                        "Leiden"=4)


    if (.type=="genes"){
      object.dist <- as.matrix(as.dist( 1-cor( object[[assay.name]]@scale.data,method = "pea" )))
    } else if (.type=="PCA"){
      object.dist <- as.matrix(as.dist( 1-cor( t(object@reductions$pca@cell.embeddings),method = "pea" )))
      numPCs <- ncol(object@reductions$pca@cell.embeddings)
      .dims = 1:max(.dims, numPCs)
    }

    if(n_cores==1){
      result <- matrix(NA, nrow = 0, ncol = 11)
      result <- as_tibble(data.frame(result))
      colnames(result) <- c(
        "index",
        "bootstrap",
        "knn",
        "resolution",
        "tot_sil_neg",
        "lowest_sil_clust",
        "max_sil_clust",
        "sil_avg",
        "sil_mean_median",
        "n_clust",
        "random.seed"
      )
      # print(system.time({
      total_iters <- length(.knns)*length(.resolutions)*.bootstraps
      if(show_progress_bar) pb = txtProgressBar(min = 0, max = total_iters, initial = 0, style = 3)
      iter_cur <- 0
      for(a_knn in .knns) {
        res = rep( .resolutions , length(.bootstraps) )
        bootstrap <- c()
        for(i in 1:length(.bootstraps)){
          bootstrap <- c(bootstrap, rep(i, length(.resolutions)))
        }

        index <- 1:length(res)

        my_sil.df <- dplyr::tibble( index = index,
                                    bootstrap = bootstrap,
                                    knn = 0,
                                    resolution = res,
                                    tot_sil_neg = 0, lowest_sil_clust = 0, max_sil_clust = 0,
                                    sil_avg = 0, sil_mean_median = 0, n_clust = 0, random.seed = 0 )
        set.seed(a_knn)
        selected_samples <- sample(colnames(object),size=n_cells_to_subsample,replace=.replace)
        x <- object[ , colnames(object) %in% selected_samples ]

        if (.type=="genes"){
          if(.pct_cells < 100){
            x.dist <- as.matrix(as.dist( 1-cor( x[[assay.name]]@scale.data,method = "pea" )))
          } else {
            x.dist <- object.dist
          }
          x@graphs <- Seurat::FindNeighbors(x.dist,
                                            distance.matrix = TRUE,
                                            verbose = FALSE,
                                            # k.param = my_sil.df$knn[a_index], # iteration over KNNs
                                            k.param = a_knn,
                                            annoy.metric = "euclidean",
                                            dims=.dims,
                                            compute.SNN = TRUE)
        } else if (.type=="PCA"){
          if(.pct_cells < 100){
            x.dist <- as.matrix(as.dist( 1-cor( t(x@reductions$pca@cell.embeddings),method = "pea" )))
          } else {
            x.dist <- object.dist
          }
          x <- Seurat::FindNeighbors(x,
                                     reduction = "pca",
                                     verbose = FALSE,
                                     # k.param = my_sil.df$knn[a_index], # iteration over KNNs
                                     k.param = a_knn,
                                     annoy.metric = "euclidean",
                                     dims=.dims,
                                     compute.SNN = TRUE)
        }

        names(x@graphs) <- c("nn", "snn")

        for ( a_index in my_sil.df$index )
        {
          set.seed(a_index)
          my_sil.df$random.seed[a_index] <- a_index
          # print(paste0("--- Silhouette score computation resolution/bootstrap : ",
          #              res[a_index] ,
          #              "|" ,
          #              my_sil.df$bootstrap[a_index] ))
          x <- Seurat::FindClusters( x,
                                     graph.name = "snn",
                                     resolution = res[a_index],
                                     verbose = FALSE,
                                     modularity.fxn = 1,
                                     algorithm = .clust_alg,
                                     random.seed = my_seed
          )
          # table(x$seurat_clusters)
          iter_cur <- iter_cur + 1
          if(show_progress_bar) setTxtProgressBar(pb, iter_cur)

          if ( nlevels(x$seurat_clusters) == 1 ) next ;
          s <- cluster::silhouette( as.integer(x$seurat_clusters) , x.dist )
          if (weights=="exp"){
            neg.sil <- (s[,"sil_width"] < 0)
            # exp_base <- exp(1)
            s[neg.sil,"sil_width"] <- -1*(exp_base^abs(s[neg.sil,"sil_width"]))
          }
          # pdf( file.path(reports.dir,paste0(“sil-res-“,a_res,“.pdf”)) )
          p <- factoextra::fviz_silhouette(s,print.summary = FALSE)
          y <- sapply( levels(p$data$cluster) , function(i) mean( p$data$sil_width[ p$data$cluster == i ] ) )
          z <- sapply( levels(p$data$cluster) , function(i) median( p$data$sil_width[ p$data$cluster == i ] ) )
          tot_sil_neg <- sapply( levels(p$data$cluster) , function(i) sum( p$data$sil_width[ p$data$cluster == i ] < 0.25 ) )
          my_sil.df$sil_avg[ a_index ] = mean(y)
          my_sil.df$sil_mean_median[ a_index ] = mean(z)
          my_sil.df$tot_sil_neg[ a_index ] = sum(tot_sil_neg)
          my_sil.df$lowest_sil_clust[ a_index ] = min(y)
          my_sil.df$max_sil_clust[a_index] = max(y)
          my_sil.df$n_clust[ a_index ] = nlevels(p$data$cluster)
          # dev.off()
          # View(my_sil.df)
        }
        my_sil.df$knn = a_knn
        result <- rbind(result, my_sil.df)
      }
      if(show_progress_bar) close(pb)
      #})) # End of print
    } else {
      myCluster <- parallel::makeCluster(n_cores,type = clust.type)
      doParallel::registerDoParallel(myCluster)
      print(system.time({
        result <- foreach::foreach( a_knn = .knns, .combine = 'rbind' ) %dopar% {
          # result <- foreach::foreach( a_res = .resolutions, .combine = 'rbind' ) %dopar% { # Iterates using a_res
          # knn = rep( .knns , length(.bootstraps) )
          # bootstrap = rep( .bootstraps, length(.knns) )
          # index <- 1:length(knn)

          res = rep( .resolutions , length(.bootstraps) )
          # bootstrap = rep( .bootstraps, length(.resolutions) )
          bootstrap <- c() #rep(NA, length(.resolutions)*length(.bootstraps))
          for(i in 1:length(.bootstraps)){
            bootstrap <- c(bootstrap, rep(i, length(.resolutions)))
          }

          index <- 1:length(res)


          # my_sil.df <- dplyr::tibble( index = index,
          #                             bootstrap = bootstrap,
          #                             knn = knn, # KNN set to values of .knn repeated X bootstraps
          #                             resolution = 0, # Resolutions set to 0
          #                             tot_sil_neg = 0, lowest_sil_clust = 0, max_sil_clust = 0,
          #                             sil_avg = 0, sil_mean_median = 0, n_clust = 0, random.seed = 0 )


          my_sil.df <- dplyr::tibble( index = index,
                                      bootstrap = bootstrap,
                                      knn = 0,
                                      resolution = res,
                                      tot_sil_neg = 0, lowest_sil_clust = 0, max_sil_clust = 0,
                                      sil_avg = 0, sil_mean_median = 0, n_clust = 0, random.seed = 0 )

          nrow(my_sil.df)

          set.seed(a_knn)
          selected_samples <- sample(colnames(object),size=n_cells_to_subsample,replace=.replace)
          x <- object[ , colnames(object) %in% selected_samples ]

          if (.type=="genes"){
            if(.pct_cells < 100){
              x.dist <- as.matrix(as.dist( 1-cor( x[[assay.name]]@scale.data,method = "pea" )))
            } else {
              x.dist <- object.dist
            }
            x@graphs <- Seurat::FindNeighbors(x.dist,
                                              distance.matrix = TRUE,
                                              verbose = FALSE,
                                              # k.param = my_sil.df$knn[a_index], # iteration over KNNs
                                              k.param = a_knn,
                                              annoy.metric = "euclidean",
                                              dims=.dims,
                                              compute.SNN = TRUE)
          } else if (.type=="PCA"){
            if(.pct_cells < 100){
              x.dist <- as.matrix(as.dist( 1-cor( t(x@reductions$pca@cell.embeddings),method = "pea" )))
            } else {
              x.dist <- object.dist
            }
            x <- Seurat::FindNeighbors(x,
                                       reduction = "pca",
                                       verbose = FALSE,
                                       # k.param = my_sil.df$knn[a_index], # iteration over KNNs
                                       k.param = a_knn,
                                       annoy.metric = "euclidean",
                                       dims=.dims,
                                       compute.SNN = TRUE)
          }

          names(x@graphs) <- c("nn", "snn")

          for ( a_index in my_sil.df$index )
          {
            set.seed(a_index)
            my_sil.df$random.seed[a_index] <- a_index
            print(paste0("--- Silhouette score computation resolution/bootstrap : ",
                         res[a_index] ,
                         "|" ,
                         my_sil.df$bootstrap[a_index] ))
            x <- Seurat::FindClusters( x,
                                       graph.name = "snn",
                                       resolution = res[a_index],
                                       verbose = FALSE,
                                       modularity.fxn = 1,
                                       algorithm = .clust_alg,
                                       random.seed = my_seed
            )
            table(x$seurat_clusters)
            if ( nlevels(x$seurat_clusters) == 1 ) next ;
            s <- cluster::silhouette( as.integer(x$seurat_clusters) , x.dist )
            if (weights=="exp"){
              neg.sil <- (s[,"sil_width"] < 0)
              # exp_base <- exp(1)
              s[neg.sil,"sil_width"] <- -1*(exp_base^abs(s[neg.sil,"sil_width"]))
            }
            # pdf( file.path(reports.dir,paste0(“sil-res-“,a_res,“.pdf”)) )
            p <- factoextra::fviz_silhouette(s,print.summary = FALSE)
            y <- sapply( levels(p$data$cluster) , function(i) mean( p$data$sil_width[ p$data$cluster == i ] ) )
            z <- sapply( levels(p$data$cluster) , function(i) median( p$data$sil_width[ p$data$cluster == i ] ) )
            tot_sil_neg <- sapply( levels(p$data$cluster) , function(i) sum( p$data$sil_width[ p$data$cluster == i ] < 0.25 ) )
            my_sil.df$sil_avg[ a_index ] = mean(y)
            my_sil.df$sil_mean_median[ a_index ] = mean(z)
            my_sil.df$tot_sil_neg[ a_index ] = sum(tot_sil_neg)
            my_sil.df$lowest_sil_clust[ a_index ] = min(y)
            my_sil.df$max_sil_clust[a_index] = max(y)
            my_sil.df$n_clust[ a_index ] = nlevels(p$data$cluster)
            # dev.off()
            # View(my_sil.df)
          }
          my_sil.df$knn = a_knn
          return(my_sil.df)
        } # End of dopar
      })) # End of print
      parallel::stopCluster(myCluster)
    }



    my_sil.df <- result
    # return(my_sil.df)
  }
  opt_params <- my_sil.df %>%
    filter(sil_avg == max(sil_avg, na.rm = TRUE)) %>%
    select(knn, resolution) %>%
    slice_head()
  opt_params$knn
  opt_params$resolution

  object <- getFinal(
    object,
    res = opt_params$resolution,
    NN = opt_params$knn,
    assay = assay.name,
    slot = "scale.data",
    reduction = .type=="PCA",
    reduction.slot = "pca",
    num.pcs = max(.dims),
    verbose = FALSE,
    clust.alg = .clust_alg,
    type.fun = "group.mean.silhouette",
    weights = weights,
    exp_base = exp_base,
    rng.seed = my_seed,
    object.dist = object.dist
  )

  Seurat::Idents(object) <- "seurat_clusters"
  # if(verbose){
  #   s <- cluster::silhouette( as.integer(object$seurat_clusters), object.dist)
  #   if (weights=="exp"){
  #     neg.sil <- (s[,"sil_width"] < 0)
  #     s[neg.sil,"sil_width"] <- -1*(exp_base^abs(s[neg.sil,"sil_width"]))
  #   }
  #   plt.sil <- factoextra::fviz_silhouette(s)
  #   print(plt.sil)
  # }
  object[[assay.name]]@misc$GS.history <- my_sil.df
  return(object)
}

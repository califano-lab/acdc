GridSearch_iter <- function(object,
                            assay.name,
                            .resolutions = seq(0.01, 2, by = 0.01),
                            .knns = seq(3, 31, by = 2),
                            .clust_alg = "Louvain",
                            .dims = NULL,
                            my_seed = 0,
                            slot = "scale.data",
                            reduction = TRUE,
                            reduction.slot = "pca",
                            type.fun = "mean.silhouette",
                            weights = "unitary",
                            exp_base = 2.718282) {
  GS_mat <- getGSMat(.resolutions, .knns, my_seed)
  .clust_alg = switch(.clust_alg, Louvain = 1, `Louvain-mult-ref` = 2,
                      SLM = 3, Leiden = 4)
  a_index <- 0
  # you first need to create the progressbar object outside the loop
  pb = txtProgressBar(min = 0, max = nrow(GS_mat), initial = 0, style = 3)
  for(k in 1:length(.knns)){
    nn <- .knns[k]
    for(j in 1:length(.resolutions)){
      res <- .resolutions[j]
      a_index <- a_index + 1
      capture.output(
        sObj_i_j <- getFinal(
          object,
          res = res,
          NN = nn,
          assay = assay.name,
          slot = slot,
          reduction = reduction,
          reduction.slot = reduction.slot,
          num.pcs = ifelse(is.null(dims), NULL, max(dims)),
          verbose = FALSE,
          clust.alg = .clust_alg,
          type.fun = type.fun,
          weights = weights,
          rng.seed = my_seed
        )
      )
      GS_mat <- addIterationDataToGSMat(GS_mat,
                                        sObj_i_j,
                                        assay.name,
                                        a_index,
                                        res,
                                        nn,
                                        weights,
                                        exp_base)
      # then inside you need to update with every iteration
      setTxtProgressBar(pb, a_index)
    }
  }
  # Remember to close the progress bar to output the newline character. From the documentation:
  # The progress bar should be closed when finished with: this outputs the final newline character.
  # Simply add at the end of your loop:
  close(pb)
  return(GS_mat)
}
getGSMat <- function(.resolutions, .knns, my_seed){
  GS_mat <- as_tibble(data.frame(matrix(0,
                                        nrow = length(.resolutions) * length(.knns),
                                        ncol = 11)))
  colnames(GS_mat) <- getGSColnames()
  GS_mat$index <- seq(1, nrow(GS_mat))
  GS_mat$bootstrap <- rep(1, nrow(GS_mat))
  GS_mat$random.seed <- rep(my_seed, nrow(GS_mat))
  return(GS_mat)
}
getGSColnames <- function(){
  c(
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
}
addIterationDataToGSMat <- function(GS_mat,
                                    sObj_i_j,
                                    assay.name,
                                    a_index,
                                    res,
                                    nn,
                                    weights,
                                    exp_base){
  GS_mat$resolution[a_index] <- res
  GS_mat$knn[a_index] <- nn
  if (nlevels(sObj_i_j$seurat_clusters) > 1){
    s <- sObj_i_j[[assay.name]]@misc$sil

    if (weights=="exp"){
      neg.sil <- (s[,"sil_width"] < 0)
      # exp_base <- exp(1)
      s[neg.sil,"sil_width"] <- -1*(exp_base^abs(s[neg.sil,"sil_width"]))
    }
    # pdf( file.path(reports.dir,paste0(“sil-res-“,a_res,“.pdf”)) )
    x <- factoextra::fviz_silhouette(s, print.summary = FALSE)
    y <- sapply(levels(x$data$cluster) , function(i)
      mean(x$data$sil_width[x$data$cluster == i]))
    z <- sapply(levels(x$data$cluster) , function(i)
      median(x$data$sil_width[x$data$cluster == i]))
    tot_sil_neg <- sapply(levels(x$data$cluster) , function(i)
      sum(x$data$sil_width[x$data$cluster == i] < 0.25))
    GS_mat$sil_avg[a_index] = mean(y)
    GS_mat$sil_mean_median[a_index] = mean(z)
    GS_mat$tot_sil_neg[a_index] = sum(tot_sil_neg)
    GS_mat$lowest_sil_clust[a_index] = min(y)
    GS_mat$max_sil_clust[a_index] = max(y)
    GS_mat$n_clust[a_index] = nlevels(x$data$cluster)
  } else {
    GS_mat$n_clust[a_index] <- 1
  }
  return(GS_mat)
}

# ACDC 🤘 - beta version 

**A**utomated **C**ommunity **D**etection of **C**ell Populations

This repo contains the current beta version of ```acdc```, an optimization-based framework to automatize clustering of cell populations from scRNA-seq data using community detection algorithms. 
```acdc``` is currently **under development** and new functionalities will be released, following completion and benchmarking. 
```acdc``` is deployed as an R package and fully compatible with ```Seurat```, the main scRNA-seq analysis pipeline in R.

<div align="center">
  <img width="240" alt="image" src="https://github.com/califano-lab/acdc-beta/assets/92543296/09feabaf-d868-48d7-b830-933210db6005">
  <img width="240" alt="image" src="https://github.com/califano-lab/acdc-beta/assets/92543296/28952fc8-841e-4d3a-80bd-d1a3a92c5a07"> 
  <img width="240" alt="image" src="https://github.com/califano-lab/acdc-beta/assets/92543296/41678fd3-c583-4b7b-939e-dbd443d44c97">
</div>

- Several graph-based clustering algorithms are available within ```acdc```, including Leiden and Louvain. 
- 2 optimization routines for parameter tuning are available, Grid Search and(generalized) Simulated Annealing.
- Optimization variables include the number of nearest neighbors, *k*, resolution, *res*, and the number of principal components, *PCs*.
- Several objective functions are available, including the Silhouette Score (default).


New releases will expand functionalities to new features, including the possibility to iteratively sub-cluster cell populations to find fine grain and biologically meaningful clustering solutions.

``` 
STAY TUNED FOR UPDATES AND NOVEL DEVELOPMENTS!🤘🏾
```

**Please, be aware that while this project is "work in progress" and outcomes are continuously benchmarked, cross-platform compability might not be guaranteed. The current beta version of `acdc` has been installed and tested on several systems running MacOS and Windows, but not on Linux-based systems. As such, currently we are not able to guarantee that all functionalities will be available on Linux-based platforms"


# Installation 
1. Start R
2. Run the following commands
```
install.packages("devtools")
devtools::install_github("califano-lab/acdc")
```


# Usage and Core functions
1. Load the acdc-beta version into R as:
```
library("acdc")
```

2. ... Start playing around! 🎸


The core functions of the package are:
1. `SAClustering`: optimizes the clustering solution to find the best set of parameters (k, resolution, PCs) using a Simulated Annealing-based optimization 
2. `GridSearch_pcs_fast`: optimizes the clustering solution to find best set of parameters (k, resolution, PCs) using a Grid Search 
3. `getFinal`: returns the optimal clustering solution with a user-defined set of parameters (k, resolution, PCs). Useful when a set of optimal parameters has been identified by `SAClustering` or `GridSearch`, and one is willing to store the optimal parameters into the Seurat object without re-running an optimization routine.


# Examples

Given a Seurat Object named `S.obj`, automatically identify an optimal clustering solution by changing the number of nearest-neighbors (NN) and the resolution (res)

1. using a Simulated Annealing-based optimization (stochastic search). Retrieve solutions in a user-defined amount of time (2 min) and automatically update `S.obj` by storing the optimal clustering solution in the `seurat_cluster` slot of the metadata: 
```
settings <- list(max.time=120) # max.time must be in s

S.obj <- SAClustering(S.obj=S.obj,
res.range=c(0.1,1),
NN.range=c(3,15),
reduction=TRUE,
control=settings)
```
The default objective function, `type.fun`, in the previous example is the average silhouette, `"mean.silhouette"`, computed across all clusters.
**Remember: the longer the time that you set, the more the number of combinations that are tested, and the higher the chance of obtaining a better solution!"**

2. using a Grid Search approach (deterministic search) span a 10x10 grid of NN and resolution values and output a tibble that enlists the clustering solution corresponding to each combination of the parameters:
```
clustering_GS <- GridSearch_pcs_fast(object = S.obj,
assay.name="RNA",
.resolutions = seq(0.1, 1.9, by = 0.2),
.bootstraps = 1,
.knns = seq(11, 101, by = 10))
```
3. Store the clustering solution obtained with a given set of parameters, e.g. NN and res, and include it in `S.obj`.
```
S.obj <- getFinal(S.obj=pbmc3k.final,
res=1,
NN=30,
reduction=TRUE)
```
The default objective function, `type.fun`, in the previous example is the mean silhouette computed across all clusters, `"mean.silhouette"`.


# Dependencies
1. `Seurat`: `acdc` is fully integrated with the Seurat pipeline for single cell analysis and leverages Seurat functions for graph construction and cell clustering. 
2. `GenSA`: for optimization based on simulated annealing
3. `factoextra`: for silhouette analysis
4. `doParallel`: to enhance parallelization
5. `foreach`: for parallelization of GridSearch
6. `cluster`: for cluster analysis
7. `dplyr`: for easier data manipulation




# References
1. Kiselev, VY, Andrews, TS, Hemberg, M. (2019) Challenges in unsupervised clustering of single-cell RNA-seq data. Nat Rev Genet 20, 273–282.
2. Blondel, V D, Guillaume, J, Lambiotte, R, Lefebvre, E (2008). Fast unfolding of communities in large networks". Journal of Statistical Mechanics: Theory and Experiment. (10) P10008.
3. Satija R, Farrell JA, Gennert D, Schier AF, Regev A (2015). “Spatial reconstruction of single-cell gene expression data.” Nature Biotechnology, 33, 495-502. 
4. Traag, V.A., Waltman, L. & van Eck, N.J. (2019) From Louvain to Leiden: guaranteeing well-connected communities. Sci Rep 9, 5233. 
5. Xiang, Y., Gubian, S., Suomela, B.P., & Hoeng, J. (2013). Generalized Simulated Annealing for Global Optimization: The GenSA Package. R J., 5, 13.


# Contacts

Please contact Alexander Wang, Luca Zanella or Alessandro Vasciaveo for doubts regarding this project (please remember that this project is still under development!).  

Alexander Wang - aw3436@cumc.columbia.edu  

Luca Zanella - lz2841@cumc.columbia.edu  

Alessandro Vasciaveo - av2729@cumc.columbia.edu  

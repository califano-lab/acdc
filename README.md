# acdc

**A**utomated **C**ommunity **D**etection of **Cell** populations

This repo contains the current beta version of ```acdc```, an optimization-based framework to automatize clustering of cell populations from scRNA-seq data using community detection algorithms. 
```acdc``` is currently **under development** and new functionalities will be released, following completion and benchmarking. 
```acdc``` is deployed as an R package and fully compatible with ```Seurat```, the main scRNA-seq analysis pipeline in R.

Several graph-based clustering algorithms are available within ```acdc```, including Leiden and Louvain. 
2 optimization routines for parameter tuning are available, Grid Search and(generalized) Simulated Annealing.
Optimization variables include the number of nearest neighbors, *k*, resolution, *res*, and the number of principal components, *PCs*.
Several objective functions are available, including the Silhouette Score (default).


New releases will expand functionalities to new features, including the possibility to iteratively sub-cluster cell populations to find fine grain and biologically meaningful clustering solutions.

``` 
STAY TUNED FOR UPDATES AND NOVEL DEVELOPMENTS!
```



# Installation 
1. Start R
2. Run the following commands
```
install.packages("devtools")
devtools::install_github("LucaZanella15/silhClust")
```


# Dependencies
1. `Seurat`: the tool is fully integrated with the Seurat pipeline for single cell analysis and leverages Seurat functions for graph construction and cell clustering. 
2. `GenSA`: for optimization based on simulated annealing
3. `factoextra`: for silhouette analysis
4. `foreach`: for parallelization of GridSearch 


# Main functions
1. `SAClustering`: optimizes the clustering solution to find the best set of parameters (k, resolution, etc) using a Simulated Annealing-based optimization 
2. `GridSearch`: optimizes the clustering solution to find best set of parameters (k, resolution, etc) using a Grid Search 
3. `getFinal`: returns the optimal clustering solution with a user-defined set of parameters (k, resolution etc). Useful when a set of optimal parameters has been identified by `SAClustering` or `GridSearch`, and one is willing to store the optimal parameters into the Seurat object without re-running an optimization routine.


# References
1. Blondel, V D, Guillaume, J, Lambiotte, R, Lefebvre, E (2008). Fast unfolding of communities in large networks". Journal of Statistical Mechanics: Theory and Experiment. (10) P10008.
2. Satija R, Farrell JA, Gennert D, Schier AF, Regev A (2015). “Spatial reconstruction of single-cell gene expression data.” Nature Biotechnology, 33, 495-502. 
3. Traag, V.A., Waltman, L. & van Eck, N.J. (2019) From Louvain to Leiden: guaranteeing well-connected communities. Sci Rep 9, 5233. 
4. Xiang, Y., Gubian, S., Suomela, B.P., & Hoeng, J. (2013). Generalized Simulated Annealing for Global Optimization: The GenSA Package. R J., 5, 13.

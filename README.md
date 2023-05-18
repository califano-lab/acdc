# acdc 🤘 - beta version 

**A**utomated **C**ommunity **D**etection of **C**ell populations

This repo contains the current beta version of ```acdc```, an optimization-based framework to automatize clustering of cell populations from scRNA-seq data using community detection algorithms. 
```acdc``` is currently **under development** and new functionalities will be released, following completion and benchmarking. 
```acdc``` is deployed as an R package and fully compatible with ```Seurat```, the main scRNA-seq analysis pipeline in R.

<div align="center">
  <img width="240" alt="image" src="https://github.com/califano-lab/acdc-beta/assets/92543296/09feabaf-d868-48d7-b830-933210db6005">
  <img width="240" alt="image" src="https://github.com/califano-lab/acdc-beta/assets/92543296/336ea80f-04a3-48d0-b96b-79d8b205b436"> 
  <img width="240" alt="image" src="https://github.com/califano-lab/acdc-beta/assets/92543296/d3e8428b-c8c3-4582-bacf-91858f4c34ed">
</div>

- Several graph-based clustering algorithms are available within ```acdc```, including Leiden and Louvain. 
- 2 optimization routines for parameter tuning are available, Grid Search and(generalized) Simulated Annealing.
- Optimization variables include the number of nearest neighbors, *k*, resolution, *res*, and the number of principal components, *PCs*.
- Several objective functions are available, including the Silhouette Score (default).


New releases will expand functionalities to new features, including the possibility to iteratively sub-cluster cell populations to find fine grain and biologically meaningful clustering solutions.

``` 
STAY TUNED FOR UPDATES AND NOVEL DEVELOPMENTS!🤘🏾
```



# Installation 
1. Start R
2. Run the following commands
```
install.packages("devtools")
devtools::install_github("califano-lab/acdc-beta")
```
3. ... Start playing around! 🎸

# Core functions
1. `SAClustering`: optimizes the clustering solution to find the best set of parameters (k, resolution, etc) using a Simulated Annealing-based optimization 
2. `GridSearch`: optimizes the clustering solution to find best set of parameters (k, resolution, etc) using a Grid Search 
3. `getFinal`: returns the optimal clustering solution with a user-defined set of parameters (k, resolution etc). Useful when a set of optimal parameters has been identified by `SAClustering` or `GridSearch`, and one is willing to store the optimal parameters into the Seurat object without re-running an optimization routine.



# Dependencies
1. `Seurat`: the tool is fully integrated with the Seurat pipeline for single cell analysis and leverages Seurat functions for graph construction and cell clustering. 
2. `GenSA`: for optimization based on simulated annealing
3. `factoextra`: for silhouette analysis
4. `foreach`: for parallelization of GridSearch 


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

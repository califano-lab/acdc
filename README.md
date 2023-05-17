# acdc

**A**utomated **C**ommunity **D**etection of **Cell** populations

This repo contains the current beta version of ```acdc```, an optimization-based framework to automatize clustering of cell populations from scRNA-seq data. 
```acdc``` is currently **under development** and new functionalities will be released, upon completion and benchmarking. 
```acdc``` is currently deployed as an R package, and is fully compatible with ```Seurat```, the main scRNA-seq analysis pipeline available inR.

Several graph-based clustering algorithms available within ```acdc``` including Leiden and Louvain. 
2 optimization routines are available for parameter tuning, a Grid Search and a (generalized) Simulated Annealing.
Suggested optimization variables include the number of nearest neighbors, *k*, resolution, *res*, and the number of principal components, *PCs*.
Several objective functions are available, including the Silhouette Score (default).


New releases will include new features, including the possibility to iteratively sub-cluster cell populations to find fine grain, biologically sound, clustering solutions.

``` 
print("\033[32mSTAY TUNED FOR UPDATES AND NOVEL DEVELOPMENTS!\033[0m")
```





# Installation 
1. Start R
2. Run the following commands
```
install.packages("devtools")
devtools::install_github("LucaZanella15/silhClust")
```


# Dependencies
1. `Seurat`: the code actually integrates with the Seurat pipeline for single cell analysis and uses Seurat functions for graph construction and cell clustering. 
2. `GenSA`: for optimization based on simulated annealing
3. `factoextra`: for silhouette analysis
4. `foreach`: for parallelization of GridSearch (In case you don't include it in the final package, remove this line)


-

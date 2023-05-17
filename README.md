# acdc

**A**utomated **C**ommunity **D**etection* of **Cell** populations

This repo contains the current beta version of ```acdc```, an optimization-based framework to highly automatize clustering of cell populations from scRNA-seq data. 
The tool is currently **under development** and new deployed as an R package

Deployed as an R package, fully compatible with Seurat

Several graph-based clustering algorithms available

Optimization framework: 
variables: k, resolution, PCs
several objective functions, e.g. Silhouette and Entropy
routines: Grid Search & (Generalized) Simulated Annealing1 	

Allows iteratively sub-clustering cell populations 						
![image](https://github.com/califano-lab/acdc-beta/assets/92543296/a352924c-296f-4ab2-ae4c-40c54dc0909f)

This 
choice of the number of clusters by optimization of the silhouette score

SAClustering is the function to run Generalized Simulated Annealing Optimization on Silhouette

Remaining functions are not intended to be part of the package





# Installation 
1. Start R
2. Run the following commands
```
install.packages("devtools")
devtools::install_github("LucaZanella15/silhClust")
```
(Currently the repo is still private, therefore please take the following action to install:)
```
devtools::install_github("LucaZanella15/silhClust", auth_token="YourToken")
```


# Dependencies
1. `Seurat`: the code actually integrates with the Seurat pipeline for single cell analysis and uses Seurat functions for graph construction and cell clustering. 
2. `GenSA`: for optimization based on simulated annealing
3. `factoextra`: for silhouette analysis
4. `foreach`: for parallelization of GridSearch (In case you don't include it in the final package, remove this line)


--------------------------------------------------------------------------------------------------

# Code

4. List of negative clusters (Thibshirani)
5. (Add automatic computation of PCA in the tool? )
7. 100,000 -> 500 -> FastClust2 -> kNN 
8. Provide additional means to compute distance
9. Provide additional objective functions

10. For Alec: enlist the new dependencies to be added to the list of Dependencies in this README 



# Others
1. (Add ... to all functions (where needed) )
2. Add output in S.obj[[...]]@misc when the algorithm returns one single cluster (to SAClustering, getFinal. Perhaps, Alec, also for the GridSearch if it is conceived in the same way)
3. Colors silhouette
4. if control = NULL, choose fast action, i.e. default parameters for fast solution (to be decided after some tests) 
5. (RunTime and other parameters added to the output )
6. ( Add progress bar in situations where useful )
7. (Rename control list)
8. (Tweak for span larger space when optmizing NN (possibly useful for other parameters) ) 






# Manuscript
3 datasets: fair comparison between grid search and SAClustering: focus on time and spanned solutions.
  - benchmark (e.g. pbmc3k). in iterClust, Ding used: pbmc (public Chromium 10x Genomics), Tirosh (SKCM), Usoskin 2015. (Mouse sensory neuron dataset), Franti et al. 2016 (Dim1024), Gionis 2007 (Aggregation dataset)
  - large dataset
  - dataset used in two other papers
  
  Figs:
  - UMAP 
  - 3D optimization 
  - Silhouette boxplots




# What comes after the manuscript
1. Add Details in help and vignettes







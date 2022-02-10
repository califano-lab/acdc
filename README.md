# silhClust
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

# Dependencies
1. `Seurat`: the code actually integrates with the Seurat pipeline for single cell analysis and uses Seurat functions for graph construction and cell clustering. 
2. `GenSA`: for optimization based on simulated annealing


--------------------------------------------------------------------------------------------------

# Code

4. List of negative clusters (Thibshirani)
5. (Add automatic computation of PCA in the tool? )
7. 100,000 -> 500 -> FastClust2 -> kNN 
9. Provide additional objective functions


# Others
1. (Add ... to all functions (where needed) )
2. Colors silhouette
3. if control = NULL, choose fast action, i.e. default parameters for fast solution (to be decided after some tests) 
4. (RunTime and other parameters added to the output )
5. ( Add progress bar in situations where useful )
6. (Rename control list)
7. (Tweak for span larger space when optmizing NN (possibly useful for other parameters) ) 






# Paper
3 datasets: fair comparison between grid search and SAClustering: focus on time and spanned solutions.
  - benchmark (e.g. pbmc3k). in iterClust, Ding used: pbmc (public Chromium 10x Genomics), Tirosh (SKCM), Usoskin 2015. (Mouse sensory neuron dataset), Franti et al. 2016 (Dim1024), Gionis 2007 (Aggregation dataset)
  - large dataset
  - dataset used in two other papers
  
  Figs:
  - UMAP 
  - 3D optimization 
  - Silhouette boxplots




# What comes after the paper
1. Add Details in help and vignettes







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

1. Add ... to all functions (where needed)
2. Convert output to Seurat Object
3. RunTime and other parameters added to the output
4. List of negative clusters (Thibshirani)
5. Add automatic computation of PCA in the tool? 
6. Add progress bar in situations where useful
7. 100,000 -> 500 -> FastClust2 -> kNN 
8. Rename control list
9. Provide additional objective functions
10. Tweak for span larger space when optmizing NN (possibly useful for other parameters)  


# Others
1. Colors silhouette
2. if control = NULL, choose fast action, i.e. default parameters for fast solution (to be decided after some tests) 



# Paper
3 datasets: fair comparison between grid search and SAClustering: focus on time and spanned solutions.
  - benchmark (e.g. pbmc3k)
  - large dataset
  - dataset used in two other papers
  
  Figs:
  - UMAP 
  - 3D optimization 
  - Silhouette boxplots




# What comes after the paper
1. Add Details in help and vignettes







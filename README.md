# silhClust
choice of the number of clusters by optimization of the silhouette score

SAClustering is the function to run Generalized Simulated Annealing Optimization on Silhouette

Remaning functions are not intended to be part of the package



# Code

1. Add ... to all functions (where needed)
2. Convert output to Seurat Object
3. RunTime and other parameters added to the output
4. Colors silhouette
5. List of negative clusters (Thibshirani)
6. Add automatic computation of PCA in the tool? 
7. Add progress bar in situations where useful
8. 100,000 -> 500 -> FastClust2 -> kNN 
9. Rename control list
10. if control = NULL, choose fast action (to be decided after some tests) 
11. Provide additional objective functions
12. 


# Paper
3 datasets: fair comparison between grid search and SAClustering: focus on time and spanned solutions.
  - benchmark (e.g. pbmc3k)
  - large dataset
  - dataset used in two other papers
  
  Figs:
  - UMAP 
  - 3D optimization 
  - Silhouette boxplots

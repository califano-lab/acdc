rm(list=ls()) # cleans workspace
#.rs.restartR() # restart memory to free RAM
gc()
cat("\014")


# Preparing folders and load test data
######
Computer <- "iMac" #Dell_Xps15" # "iMac # "CentOS"

switch(Computer,
       "Dell_Xps15"={
         CodesFolder <- paste0("/Users/",Sys.info()[["user"]],"/Desktop/")
         ScratchFolder <- paste0("W:/")
         
       },
       "iMac"={
         CodesFolder <- paste0("/Users/",Sys.info()[["user"]],"/Desktop/")
         ScratchFolder <- paste0("/Volumes/lz2841-2/")
       },
       "CentOS"={
         
         # optparse
         
         option_list = list(
           make_option(c("-p","--p_num"), type="character",default=NULL,
                       help="patient number (p_num)", metavar="chracter")
         );
         opt_parser = OptionParser(option_list=option_list);
         opt = parse_args(opt_parser);
         
         
         if (is.null(opt$p_num)){
           print_help(opt_parser)
           stop("Provide patient number p_num for infercnv analysis.n",
                call.=FALSE)
         } else {
           p_num <- as.integer(opt$p_num)
         }
         
         
         
         CodesFolder <- paste0("/ifs/scratch/c2b2/ac_lab/lz2841/")
         ScratchFolder <- CodesFolder
       }
)



# test dataset
p_num <- 5 # to test
Viper_f <- paste0(ScratchFolder,"NEPCPatients/Viper_analysis/P",p_num,"/")
S.obj <- readRDS(paste0(Viper_f,"Tumor_viper_P",p_num,"_epithelial.rds"))

DefaultAssay(S.obj) <- "Viper"




###### 
#### GenSA options and call to SA_opt function


NN_range <- c(3,30) # min, max number of NN

par <- c(0.011,3.1/NN_range[2]) # resolution, kNN: initial values
lower <- c(0.01, NN_range[1]/NN_range[2]) # LB: resolution, normalized num NN 
upper <- c(2,NN_range[2]/NN_range[2]) # UB

clust_alg <- 1 # for original Louvain algorithm

expected.val <- -1
absTol <- 1e-6

cat("Modify optimization parameters, especially temperature (not included here).\n")
control <- list(maxit=5000, # max number iterations
                #threshold.stop = expected.val + absTol, # 1 (worst) -1 (best) silhouette
                simple.function=FALSE,
                verbose=TRUE,
                smooth=FALSE, # smoothness objective function
                max.call=1e7, # max calls objective function
                max.time=10, # (s) max running time - 1 h
                verbose=TRUE,
                seed=1234
)


##### Inputs to function
# if object is Seurat
assay.name <- "Viper"
dist <- "correlation"
cells.dims <- 2 # 1: cells on row, features on columns
# 2. cells on columns, features on rows

cat("Make it generic.\n")
X <- S.obj[["Viper"]]@scale.data

# Compute distance matrix

if (cells.dims == 1){
  # cells are along rows
  d <- 1 - cor(t(X))
  # cells are along columns
} else if (cells.dims == 2){
  d <- 1 - cor(X)
}




type.fun <- "mean.silhouette" # "mean.silhouette" "median.silhouette" "group.mean.silhouette" "group.median.silhouette"


setwd(paste0(CodesFolder,"silhClust/R/"))
source("SA_tools.R")

clustering_solution <- SAClustering(S.obj=S.obj)



















































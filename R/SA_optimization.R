rm(list=ls()) # cleans workspace
#.rs.restartR() # restart memory to free RAM
gc()
cat("\014")

######
require(GenSA)
require(Seurat)

cat("Then turn it into a function.\n")


#######
# Test section on Rastrigin function


Rastrigin <- function(x){
  fn.call <<- fn.call + 1
  sum(x^2 - 10 * cos(2 * pi * x)) + 10 * length(x)
}

options(digits=10)
dimension <- 2
lower <- rep(-5.12, dimension)
upper <- rep(5.12, dimension)
set.seed(1234)
x.ini <- lower + runif(length(lower)) * (upper - lower)
expected.val <- 0
absTol <- 1e-13
fn.call <- 0
out.GenSA <- GenSA(par=NULL,
                   lower=lower,
                   upper=upper,
                   fn=Rastrigin,
                   control=list(threshold.stop=expected.val+absTol,
                                max.time=1.9,
                                verbose=TRUE))

out.GenSA[c("value","par","trace.mat","counts")]

cat("GenSA call functions", fn.calls, "times.\n")


#######
# Test section on risk allocation

library("quantmod")
tickers <- c("GE", "IBM", "JPM", "MSFT", "WMT")
getSymbols(tickers, from = "2000-12-01", to = "2010-12-31")
P <- NULL
for(ticker in tickers) {
  tmp <- Cl(to.monthly(eval(parse(text = ticker))))
  P <- cbind(P, tmp)
}
colnames(P) <- tickers
R <- diff(log(P))
R <- R[-1,]
mu <- colMeans(R)
sigma <- cov(R)
library("PerformanceAnalytics")
pContribCVaR <- ES(weights = rep(0.2, 5),
                   method = "gaussian", portfolio_method = "component",
                   mu = mu, sigma = sigma)$pct_contrib_ES


obj <- function(w) {
  fn.call <<- fn.call + 1
  if (sum(w) == 0) { w <- w + 1e-2 }
  w <- w / sum(w)
  CVaR <- ES(weights = w,
               method = "gaussian", portfolio_method = "component",
               mu = mu, sigma = sigma)
  tmp1 <- CVaR$ES
  tmp2 <- max(CVaR$pct_contrib_ES - 0.225, 0)
  out <- tmp1 + 1e3 * tmp2
  return(out)
}

set.seed(1234)
fn.call <<-0
out.GenSA <- GenSA(fn=obj,
                   lower=rep(0,5),
                   upper=rep(1,5),
                   control = list(smooth=FALSE,
                   max.call=3000))
fn.call.GenSA <- fn.call
out.GenSA$value
out.GenSA$counts
cat("GenSA call functions", fn.call.GenSA, "times.\n")
wstar.GenSA <- out.GenSA$par
wstar.GenSA <- wstar.GenSA / sum(wstar.GenSA)
rbind(tickers, round(100 * wstar.GenSA, 2))
100 * (sum(wstar.GenSA * mu) - mean(mu))



######
# Preparing folders
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



# Silhouette Optimization

# test dataset
p_num <- 5 # to test
Viper_f <- paste0(ScratchFolder,"NEPCPatients/Viper_analysis/P",p_num,"/")
S.obj <- readRDS(paste0(Viper_f,"Tumor_viper_P",p_num,"_epithelial.rds"))

DefaultAssay(S.obj) <- "Viper"

#### GenSA options

NN_range <- c(3,30) # min, max number of NN

par <- c(0.01,NN_range[1]/NN_range[2]) # resolution, kNN: initial values
lower <- c(0.01, NN_range[1]/NN_range[2]) # LB: resolution, normalized num NN 
upper <- c(2,NN_range[2]/NN_range[2]) # UB

clust_alg <- 1 # for original Louvain algorithm

expected.val <- -1
absTol <- 1e-6

cat("Modify optimization parameters, especially temperature (not included here).\n")
control <- list(
  # maxit=5000, # max number iterations
  threshold.stop = expected.val + absTol, # 1 (worst) -1 (best) silhouette
  smooth=FALSE, # smoothness objective function
  max.call=1e7, # max calls objective function
  max.time=3600 # (s) max running time - 1 h
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





##### Call GenSA
set.seed(1234)
fn.call <<- 0

out.GenSA <- GenSA(fn=obj,
                   par=par
                   lower=lower,
                   upper=upper,
                   control = control)




cat("Make also case for PCA-based obj.fn.\n Should you also add a require for Seurat and cluster and factoextra within obj?\n",
    "Clust alg for the moment is just Louvain, either remove it or consider adding many more.\n",
    "Also consider passing the parameter random.seed=some number in FindClusters.\n",
    "Also be careful with the line with 'return'.\n")

obj <- function(x, d,S.obj,NN_range, assay.name, clust_alg, type.fun){

  
  # first argument: x are the parameters SA is optimizing over
  # first element in x is resolution value, second element is num NN
  # additional aguments: additional parameters needed for computation
    # d = distance matrix
    # S.obj = Seurat object
    # min and max number of NNs
    # assay.name = perhaps needed for PCA-based clustering
    # clust_alg=Louvain, Leiden etc
    # obj.fun=for the switch case statement
  
  fn.call <<- fn.call + 1
  
  x[2] <- floor(x[2]*NN_range[2])
  
  
  
  S.obj@graphs <- Seurat::FindNeighbors(d,
                                    distance.matrix = TRUE,
                                    verbose = TRUE,
                                    k.param = x[2],
                                    annoy.metric = "euclidean",
                                    #dims=NULL,
                                    #assay=assay.name,
                                    compute.SNN = TRUE)
  
  
  names(S.obj@graphs) <- c("nn","snn")
  
  S.obj <- Seurat::FindClusters(S.obj,
                                graph.name="snn",
                                resolution=x[1],
                                verbose=FALSE,
                                modularity.fxn=1,
                                algorithm=clust_alg)
  
  if (nlevels(S.obj$seurat_clusters) == 1) next
  s <- cluster::silhouette( as.integer(S.obj$seurat_clusters) , d)
  
  
  obj.fn <- switch(type.fun,
         "mean.silhouette"={
         },
         "median.silhouette"={},
         "group.mean.silhouette"={
           obj.fn <- sapply( unique(s[,"cluster"]), 
                             function(i) mean( s[ s[,1]==i ,"sil_width"] ) )
           obj.fn <- mean(obj.fn)
         },
         "group.median.silhouette"={
           obj.fn <- sapply( unique(s[,"cluster"]), 
                             function(i) median( s[ s[,1]==i ,"sil_width"] ) )
           obj.fn <- mean(obj.fn)
         })
  

      

}





silhouette.fn <- function(){

}



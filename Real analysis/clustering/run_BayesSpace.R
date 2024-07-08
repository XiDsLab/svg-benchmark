##################################
####### BayesSpace
library(SingleCellExperiment)
library(scater)
library(Seurat)
library(ggplot2)
library(BayesSpace)
library(BiocParallel)
library(BiocGenerics)
library(parallel)
library(scran)
library(tidyverse)

## Load data
slideindex <- c(151507)
n_region <- 7
methods<-c("Binspect","spark","meringue","spatialDE","SOMDE","sparkX","sepal","scGCO","RV","dCor","HSIC")
data <- readRDS(paste0(slideindex,".list.rds"))
rownames(data[[2]]) <- data[[2]]$key
data <- CreateSeuratObject(data[[1]],meta.data = data[[2]])
data <- as.SingleCellExperiment(data)
data@colData@listData <- list(data@colData@listData,
                               row = data@colData@listData$array_row,
                               col = data@colData@listData$array_col)
data <- logNormCounts(data)

## Find HVG
#rankhvg <- modelGeneVar(data, assay.type="logcounts")
#hvg <- getTopHVGs(rankhvg, n=2000)

## BayesSpace
run_BayesSpace <- function(cl_num,genes,data){
  require(scater)
  require(BayesSpace)
  require(SingleCellExperiment)
  data <- runPCA(data, ntop = 2000, subset_row=genes, ncomponents=15, 
                  exprs_values="logcounts", BSPARAM=BiocSingular::ExactParam())
  rowData(data)[["is.HVG"]] <- rownames(data) %in% genes
  set.seed(149)
  # According to the author's instructions, platform and gamma will be tuned based on the datasets
  data <- spatialCluster(data, q=cl_num, platform="Visium", d=15,
                          init.method="mclust", model="t", gamma=3,
                          nrep=50000, burn.in=1000,
                          save.chain=FALSE)
  data@colData@listData$spatial.cluster
}

## Load genes ranked by SVG algorithms
load(paste0(slideindex,"_generank.RData"))
## Run
for (n in 1:11) {
  set.seed(2022)
  tmp<-get(paste0("rank_",methods[n]))
  tmp1<-tmp[order(tmp)]
  tmp1<-names(tmp1[1:2000])
  assign(paste0(methods[n],"info"),run_BayesSpace(n_region,genes =tmp1,data = brain))
  assign(paste0(methods[n],"info"),
         cbind.data.frame(x=brain@colData@listData[["row"]],
                          y=brain@colData@listData[["col"]],
                          cluster = factor(get(paste0(methods[n],"info")))))
  print(n)
}


save.image(paste0(slideindex,"_BayesSpace.RData"))

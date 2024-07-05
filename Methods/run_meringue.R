library(Seurat)
library(parallel)
library(dplyr)
library(tidyverse)
library(peakRAM)
## Load data after quality control
slideindex<-c(151507)
load(paste0(slideindex,".RData"))

## MERINGUE
library(MERINGUE)
tm<-peakRAM({
  set.seed(2022)
  # Remove poor datasets and genes
  mcounts <- cleanCounts(counts = counts, min.reads = 3, min.lib.size = 3, plot=FALSE,verbose=TRUE)
  minfo <- info[colnames(mcounts),]
  mat <- normalizeCounts(counts = mcounts,log=FALSE,verbose=TRUE)
  w <- getSpatialNeighbors(minfo)
  Ameringue <- getSpatialPatterns(mat, w)
  
  sig_genes_meringue <- filterSpatialPatterns(mat = mat,I = Ameringue, w = w,adjustPv = TRUE,alpha = 0.05,minPercentCells = 0.05,verbose = TRUE)
  
})

tchar<-paste0(slideindex,"_result_meringue.RData")
save(Ameringue,sig_genes_meringue,tm,file = tchar)

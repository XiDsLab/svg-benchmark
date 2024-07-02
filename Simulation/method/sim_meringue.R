library(Seurat)
library(parallel)
library(Giotto)
library(ggplot2)
library(cowplot)
library(dplyr)
library(patchwork)
library(tidyverse)
library(peakRAM)
library(SeuratData)

for (k in 1:2) {
  setwd(paste0("./",k))
  slideindex<-c("brainA")
  load("simulation_brainA.RData")
  
  ##meringue
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
  
  rm(list=ls())
  setwd("../")
}
c
save(counts,file='simulation_brainA.RData')
getwd()

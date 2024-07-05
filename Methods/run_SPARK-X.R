library(Seurat)
library(parallel)
library(dplyr)
library(tidyverse)
library(peakRAM)
## Load data after quality control
slideindex<-c(151507)
load(paste0(slideindex,".RData"))

## SPARK-X
library('SPARK')
tm<-peakRAM({
  set.seed(2022)
  ## extract the coordinates from the rawdata
  spark_tstart = Sys.time()
  sinfo<-as.matrix(info)
  
  AsparkX<-sparkx(t(counts),sinfo,numCores=20,option="mixture")
  
  sig_genes_sparkX<-subset(AsparkX[["res_mtest"]],AsparkX[["res_mtest"]]$adjustedPval<0.05)
  sig_genes_sparkX<-sig_genes_sparkX[order(sig_genes_sparkX$adjustedPval),]
})

save(AsparkX,sig_genes_sparkX,tm,file =paste0(slideindex,"_result_sparkX.RData") )


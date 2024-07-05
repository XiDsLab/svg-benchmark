library(Seurat)
library(parallel)
library(Giotto)
library(ggplot2)
library(cowplot)
library(dplyr)
library(patchwork)
library(tidyverse)
library(peakRAM)

for (k in 1:2) {
  setwd(paste0("./",k))
  slideindex<-c("brainA")
  load("simulation_brainA.RData")
  info<-as.data.frame(info)
  Data <- CreateSeuratObject(t(counts),meta.data = info)
  Data<-NormalizeData(object = Data)
  
  counts<-as.matrix(Data@assays[["RNA"]]@data)
  info<-data.frame(x = Data@meta.data[["x"]],
                   y = Data@meta.data[["y"]])
  rm(Data)
  
  ## dCor
  tm<-peakRAM({
    dCor <- function(i,data=counts,coordinate=info){
      # data is a vector of normalized counts
      # coordinate should be a matrix with each column represents an axis
      set.seed(2022)
      library(Rfast)
      tmp <- dcor.ttest(cbind(data[i,]),coordinate)
      return(tmp[4])
    }
    
    dCor_result<-matrix(data=NA,nrow = nrow(counts),ncol = 1)
    rownames(dCor_result)<-rownames(counts)
    
    cl <- makeCluster(getOption("cl.cores", 20))
    dCor_result[,1] <- unlist(parLapply(cl,1:nrow(counts),fun=dCor,data=counts,coordinate=info))
    stopCluster(cl)
  })
  
  save(dCor_result,tm, file = paste0(slideindex,"_result_dCor.RData"))
  
  ##HSIC
  tm<-peakRAM({
    dHSIC <- function(i,data=counts,coordinate=info,replicates = 100,p_value_method="gamma"){
      # p_value_method should be one of the followings:
      # "gamma" (gamma approximation based test), 
      # "permutation" (permutation test (slow)), 
      # "bootstrap" (bootstrap test (slow)) and 
      # "eigenvalue" (eigenvalue based test).
      set.seed(2022)
      library(dHSIC)
      tmp <- dhsic.test(cbind(data[i,]),coordinate,B= replicates,
                        method = p_value_method)
      return(tmp$p.value)
    }
    
    HSIC_result<-matrix(data=NA,nrow = nrow(counts),ncol = 1)
    rownames(HSIC_result)<-rownames(counts)
    
    cl <- makeCluster(getOption("cl.cores", 20))
    HSIC_result[,1] <- unlist(parLapply(cl,1:nrow(counts),fun =dHSIC,data=counts,coordinate=info,replicates = 100))
    stopCluster(cl)
  })
  
  
  save(HSIC_result,tm, file = paste0(slideindex,"_result_HSIC.RData"))
  rm(list=ls())
  setwd("../")
}

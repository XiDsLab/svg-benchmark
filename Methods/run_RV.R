library(Seurat)
library(parallel)
library(dplyr)
library(tidyverse)
library(peakRAM)
## Load data after quality control
slideindex<-c(151507)
load(paste0(slideindex,".RData"))
info<-as.data.frame(info)
Data <- CreateSeuratObject(t(counts),meta.data = info)
Data<-NormalizeData(object = Data)

counts<-as.matrix(Data@assays[["RNA"]]@data)
info<-data.frame(x = Data@meta.data[["x"]],
                 y = Data@meta.data[["y"]])
rm(Data)

## RV-correlation
tm<-peakRAM({
  RV_coef <- function(i,data=counts,coordinate=info){
    # data is a vector of normalized counts
    # coordinate should be a matrix with each column represents an axis
    set.seed(2022)
    library(FactoMineR)
    return(coeffRV(cbind(data[i,]),coordinate)$p.value)
  }
  
  RV_result<-matrix(data=NA,nrow = nrow(counts),ncol = 1)
  rownames(RV_result)<-rownames(counts)
  
  cl <- makeCluster(getOption("cl.cores", 20))
  RV_result[,1] <- unlist(parLapply(cl,1:nrow(counts),fun =RV_coef,data=counts,coordinate=info))
  stopCluster(cl)
})

save(RV_result,tm,file = paste0(slideindex,"_result_RVcor.RData"))



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


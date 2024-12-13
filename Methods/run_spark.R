library(Seurat)
library(parallel)
library(dplyr)
library(tidyverse)
library(peakRAM)
## Load data after quality control
slideindex<-c(151507)
load(paste0(slideindex,".RData"))

## SPARK
library(SPARK)
tm<-peakRAM({
  set.seed(2022)
  sinfo<-as.data.frame(info)
  spark <- CreateSPARKObject(counts = t(counts), location = sinfo[, 1:2], 
                             percentage = 0.1, min_total_counts = 10)
  #caculate total counts
  spark@lib_size <- apply(spark@counts, 2, sum)
  #estimate under null
  spark <- spark.vc(spark, covariates = NULL, lib_size = spark@lib_size, 
                    num_core = 20, verbose = T, fit.maxiter = 500)
  #p-value
  spark <- spark.test(spark, check_positive = T, verbose = T)
  #select spatial varaible genes
  Aspark <- spark@res_mtest[order(spark@res_mtest[["adjusted_pvalue"]]),]
  sig_genes_spark <- subset(spark@res_mtest,spark@res_mtest[["adjusted_pvalue"]]<0.05)
})

tchar <- paste0(slideindex,"_result_spark.RData")
save(Aspark,sig_genes_spark, tm,file = tchar)


library(Seurat)
library(parallel)
library(dplyr)
library(tidyverse)
library(peakRAM)
library(Rcpp)
library(RcppArmadillo)
library(RcppDist)
library(mclust)
library(edgeR)
library(lattice)
source("functions/Boost_Ising_function.R")

## Load data after quality control
slideindex<-c(151507)
load(paste0(slideindex,".RData"))

tm<-peakRAM({
  set.seed(2022)
  boostising_tstart = Sys.time()
  
  bcount <- data.matrix(counts)
  binfo <- data.matrix(info)
  filter_result <- filter_count(bcount, binfo, min_total = 10, min_percentage = 0.1)
  loc_f <- filter_result[[1]]
  count_f <- filter_result[[2]]
  
  #select spatial variable genes
  Aboostising <- Boost_Ising(count_f,loc_f, norm_method = 'tss', clustermethod = 'MGC')
  sig_genes_boostising <- subset(boostising,boostising$BF_neg > 150 | boostising$BF_pos > 150)
})

tchar<-paste0(slideindex[k],"_result_boostising.RData")
save(sig_genes_boostising,Aboostising,tm,file = tchar)
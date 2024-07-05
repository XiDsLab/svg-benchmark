library(Seurat)
library(parallel)
library(dplyr)
library(tidyverse)
library(peakRAM)
library(Rcpp)
library(RcppArmadillo)
library(RcppDist)
library(MASS)
source("R/boost.gp.R")
## Load data after quality control
slideindex<-c(151507)
load(paste0(slideindex,".RData"))

tm<-peakRAM({
  set.seed(2022)
  bcount <- data.matrix(counts)
  binfo <- data.matrix(info)
  colnames(binfo)<-c("x","y")
  
  #select spatial variable genes
  boostGP <- boost.gp(Y = bcount, loc = binfo)
  sig_genes_boostGP <- subset(boostGP,boostGP['pval']<0.05)
  
})

tchar<-paste0(slideindex[k],"_result_boostGP.RData")
save(sig_genes_boostGP,boostGP,tm,file = tchar)

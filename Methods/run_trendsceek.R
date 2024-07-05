library(Seurat)
library(parallel)
library(dplyr)
library(tidyverse)
library(peakRAM)
## Load data after quality control
slideindex<-c(151507)
load(paste0(slideindex,".RData"))

library('trendsceek')
tm<-peakRAM({
  set.seed(2022)
  pp<-pos2pp(info)
  min.ncells.expr = 3
  min.expr = 3
  counts_f<- genefilter_exprmat(as.matrix(t(counts)), min.expr, min.ncells.expr)
  pp<-set_marks(pp,as.matrix(counts_f),log.fcn = log10)
  
  Atrendsceek = trendsceek_test(pp, ncores =20)
  
  ##extract significant genes
  p = 0.05 ##Benjamini-Hochberg
  sig_genes_trendsceek = extract_sig_genes(Atrendsceek, p)
  
  #top <- trendstat_list[["supstats_wide"]][order(trendstat_list[["supstats_wide"]][["p.bh_markcorr"]]),]
  #sig_genes_trendsceek = sig_list[['markcorr']][, 'gene']
  
  ##show spatialDE genes' 4 summary statistic
  #plot_trendstats(trendstat_list, sig_genes_trendsceek)
  ##show spatialDE in 2D
  #pp_sig = pp_select(pp, sig_genes_trendsceek)
  #plot_pp_scatter(pp_sig, log_marks = FALSE, scale_marks = TRUE, pal.direction = -1)
})

tchar<-paste0(slideindex,"_result_trendsceek.RData")
save(Atrendsceek,sig_genes_trendsceek,tm,file = tchar)



library(Seurat)
library(parallel)
library(dplyr)
library(tidyverse)
library(peakRAM)
## Load data after quality control
slideindex<-c(151507)
load(paste0(slideindex,".RData"))

## Binspect
library(Giotto)
tm<-peakRAM({
  set.seed(2022)
  gobject <-
    createGiottoObject(raw_exprs = t(counts), spatial_locs = info)
  
  ## filter and normalize
  gobject <-
    filterGiotto(
      gobject = gobject,
      expression_threshold = 1,
      gene_det_in_min_cells = 5,
      min_det_genes_per_cell = 5,
      expression_values = c('raw'),
      verbose = T
    )
  gobject <-
    normalizeGiotto(gobject = gobject,
                    scalefactor = 6000,
                    verbose = T)
  
  ## create grid and network for prior
  gobject <-
    createSpatialGrid(
      gobject = gobject,
      sdimx_stepsize=400,
      sdimy_stepsize=400
    )
  
  annotated_grid = annotateSpatialGrid(gobject)
  gobject = createSpatialNetwork(
    gobject = gobject,
    minimum_k = 2,
    method = 'kNN',
    k = 10
  )
  
  ## use kmeans-bin
  km_spatialgenes = binSpect(gobject,
                             bin_method = 'kmeans',
                             spatial_network_name = "kNN_network",
                             do_parallel = TRUE,
                             cores = 20)
  ABinspect <- km_spatialgenes[order(km_spatialgenes[["adj.p.value"]]),]
  
  ## select spatial variable genes by adjusted p-values
  sig_genes_Binspect <-
    subset(km_spatialgenes, km_spatialgenes[["adj.p.value"]] < 0.05)
})

tchar <- paste0(slideindex,"_result_Binspect.RData")
save(ABinspect,sig_genes_Binspect, tm, file = tchar)

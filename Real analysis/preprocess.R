## Quality control
PreProcess <- function(data,
                       min_spots_gene = 400,
                       min_gene_spots_r = 0.01) {
  ## data should be list, and data[[1]] is counts matrix, data[[2]] is coordinate and other information
  ## min_spots_genes means that there is no less than this number of gene expression in each spot (the expression amount >0)
  ## min_gene_spots_r means that each gene is expressed in spots of no less than this proportion (the expression amount >0)
  ################################################################################################
  ## load count and info
  counts <- as.matrix(data[[1]])
  counts<-t(counts)
  
  info <- data.frame(x = data[[2]]$array_row,
                     y = data[[2]]$array_col)
  rownames(info) <- data[[2]]$key
  
  ## Quality control
  min_gene_spots = min_gene_spots_r * nrow(counts)
  # Remove all spots with low number of genes detected
  counts = counts[rowSums(counts > 0) > min_spots_gene,]
  # Remove all genes with low number of spots detected
  counts = counts[,colSums(counts > 0) > min_gene_spots]
  info <- info[rownames(counts), ]
  
  tmp <- as.matrix(info)
  info <- matrix(data = NA,
                 nrow = nrow(counts),
                 ncol = 2)
  rownames(info) <- rownames(counts)
  info[, 1] <- as.numeric(tmp[, 1])
  info[, 2] <- as.numeric(tmp[, 2])
  rm(tmp)
  
  return(list(counts, info))
}

## Normalize with Seurat
Normalize<- function(counts, info){
  ## counts matrix: spatial locations * genes, 
  ## info: including spatial coordinates (x,y) and other information about spatial locations
  ################################################################################################
  library(Seurat)
  info<-as.data.frame(info)
  Data <- CreateSeuratObject(t(counts),meta.data = info)
  Data<-NormalizeData(object = Data)
  
  counts<-as.matrix(Data@assays[["RNA"]]@data)
  info<-data.frame(x = Data@meta.data[["x"]],
                   y = Data@meta.data[["y"]])
  return(list(counts, info))
}


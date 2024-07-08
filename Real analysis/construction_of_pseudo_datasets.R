library(Matrix)
for (k in 1:10) {
  sed=k
  dir.create(paste0("./",k))
  ######################################################
  ### Permutation for 10X
  slideindex<-c(151507)
  
  ## Load data
  tchar<-paste0(slideindex,".list.rds")
  data<-readRDS(tchar)
  counts<-as.matrix(data[["count_mat"]])
  counts<-t(counts)
  info<-matrix(NA,nrow = nrow(counts),ncol = 2)
  rownames(info)<-data[["anno_coord"]][["key"]]
  info[,1]<-data[["anno_coord"]][["array_row"]]
  info[,2]<-data[["anno_coord"]][["array_col"]]
  
  setwd(paste0("./",k))
  
  ## Permutation & noise
  set.seed(sed)
  counts <- apply(counts, 2, function(x) {
    sample(x, size = length(x), replace = FALSE)
  })
  rownames(counts) <- rownames(info)
  data[["count_mat"]] <- as(counts, "dgCMatrix")
  data[["anno_coord"]] <- list(key = rownames(info),
                               array_row = info[, 1],
                               array_col = info[, 2])
  
  saveRDS(data,file = paste0("P",k,".list.rds"))
  
  ## Quality control
  min_spots_gene = 400
  min_gene_spots = 0.01*nrow(counts)
  # Remove all genes with low number of spots detected
  counts = counts[rowSums(counts > 0) > min_spots_gene,]
  # Remove all spots with low number of genes detected
  counts = counts[,colSums(counts > 0) > min_gene_spots]
  info<-info[rownames(counts),]
  
  rownames(counts) <- paste(info[, 1], info[, 2], sep = "x")
  tchar <- paste0("P", k, "_counts.csv")
  write.csv(counts, file = tchar)
  
  ### Noisy
  rrange<-c(0.1,0.2,0.3)
  for (r in rrange) {
    setwd("../")
    tchar<-paste0(slideindex,".list.rds")
    data<-readRDS(tchar)
    counts<-as.matrix(data[["count_mat"]])
    counts<-t(counts)
    info<-matrix(NA,nrow = nrow(counts),ncol = 2)
    rownames(info)<-data[["anno_coord"]][["key"]]
    info[,1]<-data[["anno_coord"]][["array_row"]]
    info[,2]<-data[["anno_coord"]][["array_col"]]
    
    setwd(paste0("./",k))
    
    set.seed(sed)
    n = round(nrow(counts)*r)
    for (i in 1:ncol(counts)) {
      set.seed(sed+i)
      index<-sample(1:nrow(counts), n, replace = FALSE)
      counts[index,i]<-sample(counts[index,i], size = n, replace = FALSE)
      print(i)
    }
    
    data[["count_mat"]]<-as(counts, "dgCMatrix")
    
    data[["anno_coord"]]<-list(
      key=rownames(info),
      array_row=info[,1],
      array_col=info[,2]
    )
    
    tchar <- paste0("N",r,"_",k, ".list.rds")
    saveRDS(data, file = tchar)
    
    ## Quality control
    min_spots_gene = 400
    min_gene_spots = 0.01*nrow(counts)
    # Remove all genes with low number of spots detected
    counts = counts[rowSums(counts > 0) > min_spots_gene,]
    # Remove all spots with low number of genes detected
    counts = counts[,colSums(counts > 0) > min_gene_spots]
    info<-info[rownames(counts),]
    
    rownames(counts)<-paste(info[,1],info[,2],sep="x")
    tchar<-paste0("N",r,"_",k,"_counts.csv")
    write.csv(counts,file = tchar)
  }
  
  setwd("../")
}



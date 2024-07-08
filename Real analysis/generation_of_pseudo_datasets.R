library(Matrix)
for (k in 1:10) {
  seed=k
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
  
  ## Permutation
  set.seed(seed)
  counts <- apply(counts, 2, function(x) {
    sample(x, size = length(x), replace = FALSE)
  })
  rownames(counts) <- rownames(info)
  data[["count_mat"]] <- as(counts, "dgCMatrix")
  data[["anno_coord"]] <- list(key = rownames(info),
                               array_row = info[, 1],
                               array_col = info[, 2])
  
  saveRDS(data,file = paste0("P",k,".list.rds"))

  setwd("../")
}



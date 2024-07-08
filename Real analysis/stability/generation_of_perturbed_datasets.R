library(Matrix)
for (k in 1:10) {
  seed=k
  dir.create(paste0("./",k))
  ######################################################
  ### swapping with r ratio
  rrange<-c(0.1,0.2,0.3)
  for (r in rrange) {
    slideindex<-c(151507)
    tchar<-paste0(slideindex,".list.rds")
    data<-readRDS(tchar)
    counts<-as.matrix(data[["count_mat"]])
    counts<-t(counts)
    info<-matrix(NA,nrow = nrow(counts),ncol = 2)
    rownames(info)<-data[["anno_coord"]][["key"]]
    info[,1]<-data[["anno_coord"]][["array_row"]]
    info[,2]<-data[["anno_coord"]][["array_col"]]
    
    setwd(paste0("./",k))
    set.seed(seed)
    n = round(nrow(counts)*r)
    for (i in 1:ncol(counts)) {
      set.seed(seed+i)
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
  }
  setwd("../")
}



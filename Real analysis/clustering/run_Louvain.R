library(Seurat)
slideindex<-c(151507)
brain <- readRDS(paste0(slideindex,".list.rds"))
rownames(brain[[2]]) <- brain[[2]]$key
brain <- CreateSeuratObject(brain[[1]],meta.data = brain[[2]])
brain <- SCTransform(brain, assay = "RNA", verbose = FALSE,return.only.var.genes = FALSE)

## load SV genes detected by methods
load(paste0(slideindex,"_generank.RData"))
methods<-c("Binspect","spark","meringue","spatialDE","SOMDE","sparkX","sepal","scGCO","RV","dCor","HSIC")

#region information for the data
Iregion<-as.matrix(brain@meta.data[["layer_guess_reordered"]])
rownames(Iregion)<-rownames(brain@meta.data[["key"]])
nregion=nlevels(brain@meta.data[["layer_guess_reordered"]])

#create ARI matrix
ARI<-matrix(data = NA,nrow = 12,ncol = 4)
row.names(ARI)<-c(methods,"region")
colnames(ARI)<-c("methods","n-region","resolution","ARI")

for (n in 1:11) {
  r1=0.2
  r2=2.8
  tmp<-get(paste0("rank_",methods[n]))
  tmp1<-tmp[order(tmp)]
  svg<-names(tmp1[1:2000])
  while (abs(r1 - r2) > 1e-6) {
    mid <- (r1 + r2) / 2
    brain2 <- ScaleData(brain, features = svg)
    brain2 <- RunPCA(brain2, verbose = FALSE, features = svg)
    brain2 <- FindNeighbors(object = brain2, dims = 1:30)
    brain2 <- FindClusters(object = brain2, resolution = mid)
    brain2 <- RunUMAP(brain2, reduction = "pca", dims = 1:30)
    Isvg <- as.matrix(list(brain2@active.ident)[[1]])
    if(nlevels(factor(Isvg))==nregion){
      ARI[n,1]<-methods[n];ARI[n,2]<-nlevels(factor(Isvg));ARI[n,3]<-mid
      assign(paste0(methods[n],"info"),
             cbind.data.frame(x=brain@meta.data[["array_row"]],
                              y=brain@meta.data[["array_col"]],
                              cluster = factor(Isvg)))
      library(mclust)
      x<-which(brain@meta.data$layer_guess_reordered!="NA's")
      ARI[n,4] <- adjustedRandIndex(Iregion[x], Isvg[x])
      break
    }else if (nlevels(factor(Isvg))<nregion){
      r1<-mid
    }else{
      r2<-mid
    }
  }
  
  if(abs(r1 - r2) < 1e-6){
    ARI[n,1]<-methods[n];ARI[n,2]<-nlevels(factor(Isvg));ARI[n,3]<-mid
    assign(paste0(methods[n],"info"),
           cbind.data.frame(x=brain@meta.data[["array_row"]],
                            y=brain@meta.data[["array_col"]],
                            cluster = factor(Isvg)))
    library(mclust)
    x<-which(brain@meta.data$layer_guess_reordered!="NA's")
    ARI[n,4] <- adjustedRandIndex(Iregion[x], Isvg[x])
  }
}

ARI[12,1]<-c("true");ARI[12,2]<-nregion;ARI[12,4]<-1
write.csv(ARI,file=paste0(slideindex[k],"_Louvain.csv"))
save.image(paste0(slideindex[k],"_Louvain.RData"))

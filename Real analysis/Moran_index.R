library(Seurat)
slideindex<-c(151507)
fchar<- paste0(slideindex,".list.rds")
brain <- readRDS(fchar)
counts<-as.matrix(brain[["count_mat"]])
info<-brain[["anno_coord"]]
info<-as.data.frame(info)
Data <- CreateSeuratObject(counts,meta.data = info)
Data<-NormalizeData(object = Data)

counts<-as.matrix(Data@assays[["RNA"]]@data)
rm(Data,info)
load(paste0(slideindex[k],"_svg.RData"))

library(spdep)
x <- brain[["anno_coord"]][["array_row"]]
y <- brain[["anno_coord"]][["array_col"]]
coords <- cbind(x, y)
dk<-knn2nb(knearneigh(coords))
r <- max(unlist(nbdists(dk, coords)))
my.nb <- dnearneigh(coords,0,r)
pM<-c()
SM<-c()
for (i in seekindex) {
  obs <- counts[i,]
  mytest<-moran.test(obs, nb2listw(my.nb, style = "W",zero.policy = TRUE),zero.policy=TRUE)
  SM<-c(SM,as.numeric(mytest[["estimate"]][["Moran I statistic"]]))
  pM<-c(pM,as.numeric(mytest[["p.value"]]))
  print(which(seekindex==i))
}

GMoran<-data.frame(seekindex=seekindex,moran.index=SM,p.value=pM)
save(GMoran,file = paste0(slideindex,"_Mgene.RData"))

p<-p.adjust(GMoran$p.value,method="BH",n=length(GMoran$p.value))
sliver<-unique(GMoran$seekindex[which(p<0.05)])

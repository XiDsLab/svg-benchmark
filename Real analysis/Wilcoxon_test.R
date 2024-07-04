library(Seurat)

slideindex<-c(151507)
## data input
fchar<- paste0(slideindex,".list.rds")
brain <- readRDS(fchar)

## delete spot whose region annotation is "NA" 
eindex<-which(brain[["anno_coord"]][["layer_guess_reordered_short"]]!="NA's")
info<-brain[["anno_coord"]][eindex,]
counts<-as.matrix(brain[["count_mat"]])[,eindex]

info<-as.data.frame(info)
Data <- CreateSeuratObject(counts,meta.data = info)
Data<-NormalizeData(object = Data)

counts<-as.matrix(Data@assays[["RNA"]]@data)
rm(Data,info)

## perform wilcox-test to find genes
bcounts<-as.data.frame(t(counts))
bcounts[,"region"]<-as.factor(brain[["anno_coord"]][["layer_guess_reordered_short"]][eindex])
comparisons <- combn(unique(as.character(bcounts[,"region"])), 2)

load(paste0(slideindex[k],"_svg.RData"))
p=matrix(data=NA,nrow = length(seekindex),ncol = ncol(comparisons))
rownames(p)<-seekindex
for (i in seekindex) {
  p[i,]<-apply(comparisons, 2, function(x) {
    wilcox.test(bcounts[,i ] ~ region, data = bcounts, subset = region %in% x)$p.value
  })
  print(which(seekindex==i))
}
tmp <- p.adjust(as.vector(p), method = "bonferroni")
p_adjusted <- matrix(tmp, nrow = length(seekindex))
rownames(p_adjusted)<-seekindex

save(seekindex,p,p_adjusted,file = paste0(slideindex,"_Wgenes.RData"))

library(ICSKAT)
p_combined<-c()
for (km in seekindex) {
  p_combined<-c(p_combined,ACAT(Pvals=p[km,][!is.na(p[km,])]))
}
p_adjusted<-p.adjust(as.vector(p_combined), method = "bonferroni")
sliver<-unique(seekindex[which(p_adjusted<0.05)])

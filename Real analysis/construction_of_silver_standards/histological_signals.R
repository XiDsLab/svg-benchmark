library(Seurat)
slideindex<-c("FFPE_Human_Breast_Cancer")

#Lchannel obtained from Lchannel.py
Lchannel <- read.csv(paste0(slideindex,"_Lchannel.csv"),header=T,row.names = 1)
Lchannel <-Lchannel$Lchannel
brain<-readRDS(paste0(slideindex,".list.rds"))
load(paste0(slideindex[k],"_svg.RData"))
counts<-as.matrix(brain[["count_mat"]])
info<-data.frame(x=brain[["anno_coord"]][["array_row"]],
                 y=brain[["anno_coord"]][["array_col"]])
Data <- CreateSeuratObject(counts,meta.data = info)
Data<-NormalizeData(object = Data)

counts<-as.matrix(Data@assays[["RNA"]]@data)
rm(Data,info)


cor=c()
p=c()
for (i in seekindex) {
  tmp<-cor.test(counts[i,],Lchannel)
  cor<-c(cor,tmp[["estimate"]][["cor"]])
  p<-c(p,tmp[["p.value"]])
  print(which(seekindex==i))
}
qvalue<-p.adjust(p,method="BH",n=length(p))
Lresult<-data.frame(gene = seekindex,
                    cor = cor,
                    p.value = p,
                    p.adj = qvalue)

save(Lresult,file = paste0(slideindex,"_Lcor.RData"))

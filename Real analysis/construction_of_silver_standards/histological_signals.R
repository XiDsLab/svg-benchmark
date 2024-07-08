library(Seurat)
slideindex<-c("FFPE_Human_Breast_Cancer")

#Lchannel obtained from LABTrans.py
Lchannel <- read.csv(paste0(slideindex,"_Lchannel.csv"),header=T,row.names = 1)
Lchannel <-Lchannel$Lchannel
data<-readRDS(paste0(slideindex,".list.rds"))

## load genes after quality control (seekindex)
load(paste0(slideindex[k],"_svg.RData"))
counts<-as.matrix(data[["count_mat"]])[seekindex,]
info<-data.frame(x=data[["anno_coord"]][["array_row"]],
                 y=data[["anno_coord"]][["array_col"]])
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
sliver<-data.frame(gene = seekindex,
                    cor = cor,
                    p.value = p,
                    p.adj = qvalue)

library(Seurat)
slideindex<-c("FFPE_Human_Breast_Cancer")

#Lchannel obtained from LABTrans.py
Lchannel <- read.csv(paste0(slideindex,"_Lchannel.csv"),header=T,row.names = 1)
Lchannel <-Lchannel$Lchannel
data<-readRDS(paste0(slideindex,".list.rds"))

## load genes after quality control (seekindex) or perform quality control with PreProcess
#load(paste0(slideindex[k],"_svg.RData"))
#counts<-as.matrix(data[["count_mat"]])[seekindex,]
tmp<-PreProcess(data=data,mat_form = T)
counts<-as.matrix(tmp[[1]])
seekindex<-rownames(counts)
info<-data[["anno_coord"]]
info<-as.data.frame(info)
Data <- CreateSeuratObject(t(counts),meta.data = info)
Data<-NormalizeData(object = Data)
counts<-as.matrix(Data@assays[["RNA"]]@data)
rm(Data,info,tmp)

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
sliver<-seekindex[which(qvalue<0.05)]

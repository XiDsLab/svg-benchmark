library(Seurat)
source("../preprocess.R")
## Load data
slideindex<-c(151507)
fchar<- paste0(slideindex,".list.rds")
data <- readRDS(fchar)

##load genes after quality control (seekindex) or perform quality control with PreProcess
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
rm(Data,info)

#############################################################################
## For datasets with less than 20,000 spatial locations
library(spdep)
x <- data[["anno_coord"]][["array_row"]]
y <- data[["anno_coord"]][["array_col"]]
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
p.adj<-p.adjust(GMoran$p.value,method="BH",n=length(GMoran$p.value))
sliver<-GMoran$seekindex[which(p.adj<0.05)]

#############################################################################
## For datesets with more than 20,000 spatial locations
library(parallel)
x <- as.numeric(data[[2]][,1])
y <- as.numeric(data[[2]][,2])
rm(data)
#library(devtools)
#install_github('mcooper/moranfast')

pM <- function(i,data=counts,co1=x,co2=y){
  # data is a vector of normalized counts
  # coordinate should be a matrix with each column represents an axis
  library(moranfast)
  mytest<-moranfast(data[i,], co1, co2)
  S<-as.numeric(mytest[["observed"]])
  p<-as.numeric(mytest[["p.value"]])
  return(c(S,p))
}

GMoran<-matrix(data=NA,nrow = length(seekindex),ncol = 3)
rownames(GMoran)<-seekindex
colnames(GMoran)<-c("moran.index","p.value","p.adj")

cl <- makeCluster(getOption("cl.cores", 20))
res<- unlist(parLapply(cl,1:length(seekindex),fun=pM,data=counts,co1=x,co2=y))
stopCluster(cl)
GMoran[,1]<-res[seq(1,length(res),2)]
GMoran[,2]<-res[seq(0,length(res),2)]
GMoran[,3]<-p.adjust(GMoran[,2],method="BH",n=length(seekindex))
sliver<-rownames(GMoran)[which(GMoran[,3]<0.05)]


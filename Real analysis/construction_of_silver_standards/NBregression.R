library(Seurat)
library(car)
library(scran)
library(caret)
library(MASS)

## slide ID
slideindex<-c(151507)

## data input
data <- readRDS(paste0(slideindex,".list.rds"))

## delete spot whose region annotation is "NA" 
eindex<-which(data[["anno_coord"]][["layer_guess_reordered_short"]]!="NA's")
data[["anno_coord"]]<-data[["anno_coord"]][eindex,]
counts<-as.matrix(data[["count_mat"]])[,eindex]
load(paste0(slideindex[k],"_svg.RData"))

## perform ANOVA to find genes
bcounts<-as.data.frame(t(counts+1))
bcounts[,"region"]<-brain[["anno_coord"]][["layer_guess_reordered_short"]]
dmy<-dummyVars(~region,data=bcounts)
trfs<-data.frame(predict(dmy,newdata = bcounts))
bcounts<-cbind(bcounts,trfs)
rm(dmy)
rm(counts)

p=c()
for (i in seekindex) {
  g1 = glm.nb(
    bcounts[,i ] ~ region.L1 + region.L2 + region.L3 + region.L4 + region.L5 + region.L6 +region.WM,
    data = bcounts,
    link = "log",
  )
  g2 = glm.nb(
    bcounts[,i ] ~ 1,
    data = bcounts,
    link = "log",
  )
  kr<-anova(g1, g2)$`Pr(Chi)`[2]
  p<-c(p,kr)
  print(which(seekindex==i))
}

save(seekindex,p,file = paste0(slideindex,"_Agenes.RData"))

p<-p.adjust(p,method="BH",n=length(p))
sliver<-unique(seekindex[which(p<0.05)])

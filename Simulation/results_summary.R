#install.packages('PRROC')
library(cpm)
library(Seurat)
library(ggplot2)
library(parallel)
library(PRROC)
setwd('result_path')

## Load results of the SVG methods--------------------------------------------
slideindex<-c("brainA")
k=1
methods<-c("Binspect","spark","meringue","spatialDE","SOMDE","sparkX","sepal","scGCO","RVcor","dCor","HSIC")
load('simulation_brainA.Rdata')
for (llindex in 1:3) {
  load(paste0(slideindex[k],"_result_",methods[llindex],".RData"))
}
sig_genes_spatialDE<-read.csv(paste0(slideindex[k],"_result_",methods[4],".csv"), header=T,row.names = 1)
sig_genes_SOMDE<-read.csv(paste0(slideindex[k],"_result_",methods[5],".csv"), header=T,row.names = 1)
load(paste0(slideindex[k],"_result_",methods[6],".RData"))
sig_genes_sepal<-read.table(paste0(slideindex[k],"_result_",methods[7],".tsv"), sep=",",header=T,row.names = 1)
sig_genes_scGCO<-read.csv(paste0(slideindex[k],"_results_",methods[8],".csv"), header=T,row.names = 1)
for (llindex in 9:11) {
  load(paste0(slideindex[k],"_result_",methods[llindex],".RData"))
}

## data process
sig_genes_Binspect<-sig_genes_Binspect[order(sig_genes_Binspect$adj.p.value),]
sig_genes_spark<-sig_genes_spark[order(sig_genes_spark$adjusted_pvalue),]
sig_genes_spatialDE<-sig_genes_spatialDE[sig_genes_spatialDE$g!="log_total_count",]
sig_genes_spatialDE<-subset(sig_genes_spatialDE,sig_genes_spatialDE$qval<0.05)
sig_genes_spatialDE<-sig_genes_spatialDE[order(sig_genes_spatialDE$qval),]
sig_genes_SOMDE<-sig_genes_SOMDE[sig_genes_SOMDE$g!="log_total_count",]
sig_genes_SOMDE<-subset(sig_genes_SOMDE,sig_genes_SOMDE$qval<0.05)
sig_genes_SOMDE<-sig_genes_SOMDE[order(sig_genes_SOMDE$qval),]
sig_genes_sparkX<-rownames(sig_genes_sparkX)
positive<-subset(sig_genes_sepal,sig_genes_sepal$average>0)
y<-sort(positive$average)
x<-c(1:length(y))*0.01
cut<-processStream(y, cpmType = "Exponential")$changePoints
point<-y[quantile(cut, 0.9)]

#plot(x,y,cex=0.9)
#abline(v=quantile(cut, 0.9)*0.01,col='blue')

sig_genes_sepal<-subset(sig_genes_sepal,sig_genes_sepal$average>point)
sig_genes_scGCO<-sig_genes_scGCO[,1:4]
sig_genes_scGCO<-subset(sig_genes_scGCO,sig_genes_scGCO$fdr<0.05)
sig_genes_scGCO<-sig_genes_scGCO[order(sig_genes_scGCO$fdr),]

## overall SV genes
spark<-rownames(sig_genes_spark)
Binspect<-c(sig_genes_Binspect[,1])
Binspect<-Binspect[["genes"]]
meringue<-sig_genes_meringue
spatialDE<-c(sig_genes_spatialDE["g"])
spatialDE<-spatialDE[["g"]]
SOMDE<-c(sig_genes_SOMDE["g"])
SOMDE<-SOMDE[["g"]]
sparkX<-sig_genes_sparkX
sepal<-rownames(sig_genes_sepal)
scGCO<-rownames(sig_genes_scGCO)

nameindex<-c("RV","dCor","HSIC")
for (n in 1:3) {
  x<-p.adjust(get(paste0(nameindex[n],"_result")),method="BH",n=length(get(paste0(nameindex[n],"_result"))))
  x<-as.matrix(x)
  rownames(x)<-rownames(get(paste0(nameindex[n],"_result")))
  x<-subset(x,x<0.05)
  x<-x[order(x[,1]),]
  x<-names(x)
  assign(nameindex[n],x)
}

## seekindex including all genes after QC
save(Binspect,spark,meringue,spatialDE,SOMDE,sparkX,sepal,scGCO,RV,dCor,HSIC,seekindex,
     file = paste0(slideindex[k],"_svg.RData"))

methods<-c("Binspect","spark","meringue","spatialDE","SOMDE","sparkX","sepal","scGCO","RVcor","dCor","HSIC")

for (llindex in 1:3) {
  load(paste0(slideindex[k],"_result_",methods[llindex],".RData"))
}
sig_genes_spatialDE<-read.csv(paste0(slideindex[k],"_result_",methods[4],".csv"), header=T,row.names = 1)
sig_genes_SOMDE<-read.csv(paste0(slideindex[k],"_result_",methods[5],".csv"), header=T,row.names = 1)
load(paste0(slideindex[k],"_result_",methods[6],".RData"))
sig_genes_sepal<-read.table(paste0(slideindex[k],"_result_",methods[7],".tsv"), sep=",",header=T,row.names = 1)
sig_genes_scGCO<-read.csv(paste0(slideindex[k],"_results_",methods[8],".csv"), header=T,row.names = 1)
for (llindex in 9:11) {
  load(paste0(slideindex[k],"_result_",methods[llindex],".RData"))
}

## data process
for (nm in c(1:3,6)) {
  tmp<-get(paste0("A",methods[nm]))
  if(nm==1){
    tmp1<-as.vector(tmp$adj.p.value)
    names(tmp1)<-tmp$genes
    dif<-setdiff(seekindex,tmp$genes)
  }
  if(nm==2){
    tmp1<-as.vector(tmp$adjusted_pvalue)
    names(tmp1)<-rownames(tmp)
    dif<-setdiff(seekindex,rownames(tmp))
  }
  if(nm==3){
    tmp1<-as.vector(tmp$p.adj)
    names(tmp1)<-rownames(tmp)
    dif<-setdiff(seekindex,rownames(tmp))
  }
  if(nm==6){
    tmp1<-as.vector(tmp[["res_mtest"]][["adjustedPval"]])
    names(tmp1)<-rownames(tmp[["res_mtest"]])
    dif<-setdiff(seekindex,rownames(tmp[["res_mtest"]]))
  }
  
  if(length(dif)!=0){
    tmp2<-as.vector(rep(2,length(dif)))
    names(tmp2)<-dif
  }else{
    tmp2<-c()
  }
  tmpc<-c(tmp1,tmp2)
  tmpc<-rank(tmpc)
  assign(paste0("rank_",methods[nm]),tmpc)
}

for (nm in c(4,5,7,8)) {
  tmp<-get(paste0("sig_genes_",methods[nm]))
  if(nm==4){
    tmp<-tmp[tmp$g!="log_total_count",]
    tmp1<-as.vector(tmp$qval)
    names(tmp1)<-tmp$g
    dif<-setdiff(seekindex,tmp$g)
  }
  if(nm==5){
    tmp1<-as.vector(tmp$qval)
    names(tmp1)<-tmp$g
    dif<-setdiff(seekindex,tmp$g)
  }
  if(nm==7){
    tmp1<-as.vector(tmp$average)
    names(tmp1)<-rownames(tmp)
    dif<-setdiff(seekindex,rownames(tmp))
  }
  if(nm==8){
    tmp1<-as.vector(tmp$fdr)
    names(tmp1)<-rownames(tmp)
    dif<-setdiff(seekindex,rownames(tmp))
  }
  
  if(length(dif)!=0){
    if(nm!=7){
      tmp2<-as.vector(rep(2,length(dif)))
      names(tmp2)<-dif
    }else{
      tmp2<-as.vector(rep(-1,length(dif)))
      names(tmp2)<-dif
    }
  }else{
    tmp2<-c()
  }
  tmpc<-c(tmp1,tmp2)
  if(nm!=7){
    tmpc<-rank(tmpc)
  }else{
    tmpc<-rank(-tmpc)
  }
  assign(paste0("rank_",methods[nm]),tmpc)
}

nameindex<-c("RV","dCor","HSIC")
for (nm in 1:3) {
  tmp<-as.vector(p.adjust(get(paste0(nameindex[nm],"_result")),method="BH",n=length(get(paste0(nameindex[nm],"_result")))))
  names(tmp)<-rownames(get(paste0(nameindex[nm],"_result")))
  tmp<-rank(tmp)
  assign(paste0("rank_",nameindex[nm]),tmp)
}
#setwd("../")
save(rank_Binspect,rank_spark,rank_meringue,rank_spatialDE,rank_SOMDE,rank_sparkX,
     rank_sepal,rank_scGCO,rank_RV,rank_dCor,rank_HSIC,seekindex,
     file = paste0(slideindex[k],"_generank.RData"))


## Calculate indicators based on simulated gold standards-------------------
methods<-c("Binspect","spark","meringue","spatialDE","SOMDE","sparkX","sepal","scGCO","RV","dCor","HSIC")
Information<-matrix(data=NA,nrow=11,ncol=10)
rownames(Information)<-methods
colnames(Information)<-c("precision","sensitivity","accuracy","similarity","F1",
                         "num_sliver","AUPR","AUROC","EP","EPratio")

#####Load true svg (in simulation) to L
ncol(counts)
# G<-colnames(counts)
# hah<-sapply(1:ncol(counts), function(t){
#   grepl("DE",G[t])
# })
L<-DE_gen
load(paste0(slideindex[k],"_svg.RData"))
for (tt in 1:11) {
  cat(methods[tt],length(get(methods[tt])),"\n")
  Information[tt,1]=length(intersect(L,get(methods[tt])))/length(get(methods[tt]))
  Information[tt,2]=length(intersect(L,get(methods[tt])))/length(L)
  Information[tt,3]=(length(seekindex)-length(L)-length(get(methods[tt]))+length(intersect(L,get(methods[tt]))))/length(seekindex)
  Information[tt,4]=length(intersect(L,get(methods[tt])))/length(union(L,get(methods[tt])))
  Information[tt,5]=2*Information[tt,1]*Information[tt,2]/(Information[tt,1]+Information[tt,2])
  Information[tt,6]=length(L)
}

load(paste0(slideindex[k],"_generank.RData"))

nonsvg<-setdiff(seekindex,L)
for (tt in 1:11) {
  tmp<-get(paste0("rank_",methods[tt]))
  pr <- pr.curve(scores.class0 = -tmp[L], scores.class1 = -tmp[nonsvg], curve = T)
  Information[tt,7]=pr$auc.integral
  #plot(pr)
  roc <- roc.curve(scores.class0 = -tmp[L], scores.class1 = -tmp[nonsvg], curve = T)
  Information[tt,8]=roc$auc
  #plot(roc)
  tmp1<-tmp[order(tmp)]
  tmp1<-names(tmp1[1:length(L)])
  Information[tt,9]=length(intersect(L,tmp1))/length(tmp1)
  Information[tt,10]=length(intersect(L,tmp1))/length(tmp1)/(length(L)/length(seekindex))
}

Information
write.csv(Information,paste0(slideindex[k],"_sim.csv"))

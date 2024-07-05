## load reference data------------------------------------------------------
library(ggplot2)
library(SRTsim)
load('reference/reference2.RData')

## Estimate model parameters for data generation----------------------------
ngen<-10000
set.seed(1)
counts<-t(counts)
de_potion<-0.3
non_potion<-0.15
colnames(info) <- c("x","y","label")

simSRT  <- createSRT(count_in=counts,loc_in =info)
simSRT2 <- srtsim_fit(simSRT,sim_schem="domain")
save(simSRT2,file='sim_result.RData')
load('sim_result.RData')
simSRT3<-simSRT2
gene<-rownames(counts)

## Calculate indicators-----------------------------------------------------
### Kruskal-Wallis test statistic
zhi<-c()
counts<-t(counts)
KW<-sapply(1:ncol(counts),function(t){
  data<-counts[,t]
  label<-info$kind
  R<-0
  indx<-order(data)
  zhi[indx]<-c(1:length(data))
  for(i in 0:(length(unique(label))-1)){
    n_i<-length(which(label==as.character(i)))
    re<-sum(data[which(label==as.character(i))])
    R<-R+re^2/n_i
  }
  return(R)
})
names(KW)<-colnames(counts)
info<-as.data.frame(info)
names(info)<-c('x','y','kind')
domain<-as.character(unique(info$kind))

### The proportion of spots in each domain to the overall spot
domain_pro<-sapply(1:length(domain), function(t){
  length(which(info$kind==domain[t]))/nrow(info)
})


### For the SRT method, parameter estimation is not performed in regions with low gene expression, and it is necessary to identify which genes do not have parameter estimation in a certain region
DE_gen<-c()
nonexp_gen<-c()
for (i in 1:length(domain)){
  nonexp_gen<-c(nonexp_gen,names(simSRT2@EstParam[[domain[i]]]@listData[["gene_sel2"]]))
}
nonexp_gen<-unique(nonexp_gen)
gener_gene<-setdiff(gene,nonexp_gen)

### The mean of each gene in each region
domain_mean<-matrix(0,ncol=length(gene),nrow=length(domain))
for (i in 1:length(domain)){
  indx<-match(names(simSRT2@EstParam[[domain[i]]]@listData[["gene_sel1"]]),gene)
  domain_mean[i,indx]<-simSRT2@EstParam[[domain[i]]]@listData[["marginal_param1"]]@listData[["mu"]]
}
colnames(domain_mean)<-gene
rownames(domain_mean)<-domain

### Caculate differences in mean expression between regions for each gene
differ<-sapply(1:ncol(domain_mean),function(t){
  mu<-domain_mean[,t]%*%domain_pro
  re<-(domain_mean[,t]-mu)/mu
  #re<-re*domain_pro
  return(re%*%re)
})
names(differ)<-colnames(domain_mean)
differ_non<-differ[match(nonexp_gen,names(differ))]
names(differ_non)<-nonexp_gen
differ_ger<-differ[match(gener_gene,names(differ))]
names(differ_ger)<-gener_gene

### Caculate global mean for each gene
global_mean<-sapply(1:ncol(domain_mean),function(t){
  mu<-domain_mean[,t]%*%domain_pro
  return(mu)
})
names(global_mean)<-colnames(domain_mean)

### Calculate propotion points with non zero expression in each region
domain_exp<-matrix(0,nrow(domain_mean),ncol(domain_mean))
for (i in 1:length(domain)){
  indx<-which(info$kind==domain[i])
  counts_use<-counts[,indx]
  domain_exp[i,]<-rowSums(counts_use>0)/ncol(counts_use)
}
colnames(domain_exp)<-colnames(domain_mean)
rownames(domain_exp)<-rownames(domain_mean)

###Calculate the maximum proportion of regional expression for each gene
svg_pot_max<-sapply(1:ncol(domain_exp), function(t){
  return(max(domain_exp[,t]))
})
names(svg_pot_max)<-colnames(domain_mean)

###Calculate the minimum proportion of regional expression for each gene
svg_pot_min<-sapply(1:ncol(domain_exp), function(t){
  return(min(domain_exp[,t]))
})
names(svg_pot_min)<-colnames(domain_mean)

info_gene<-matrix(0,ncol(domain_mean),9)
colnames(info_gene)<-c('gener_gen','differ','differ_rank','max_spot','min_spot','global_mean','DE_gen',"final_gen",'KW')
rownames(info_gene)<-colnames(domain_mean)
info_gene<-as.data.frame(info_gene)
info_gene$gener_gen[match(gener_gene,rownames(info_gene))]<-1
info_gene$differ<-differ[match(names(differ),rownames(info_gene))]
info_gene$differ_rank[rev(order(differ))]<-c(1:nrow(info_gene))
info_gene$max_spot<-svg_pot_max
info_gene$min_spot<-svg_pot_min
info_gene$global_mean<-global_mean
info_gene$KW<-KW

## Filter the required SVG and non SVG based on indicators------------------
DE_gen<-rownames(info_gene)[which((info_gene$gener_gen==1&info_gene$differ>18&info_gene$max_spot>0.05&info_gene$global_mean>0.005)|(info_gene$gener_gen==0&info_gene$max_spot>0.05))]
normal_gen_cho<-setdiff(rownames(info_gene)[which((info_gene$gener_gen==1&info_gene$differ<18&info_gene$global_mean>0.001)|(info_gene$gener_gen==1&info_gene$differ>18&info_gene$max_spot<0.01&info_gene$global_mean>0.005)|(info_gene$gener_gen==0&info_gene$max_spot<0.05&info_gene$max_spot>0.01))],DE_gen)
DE_gen<-sample(DE_gen,size=de_potion*ngen,replace=FALSE)
normal_gene<-sample(normal_gen_cho,size=ngen-length(DE_gen),replace = FALSE)

gen_final<-c(DE_gen,normal_gene)
gen_final<-sample(gen_final,length(gen_final),replace = FALSE)

info_gene$DE_gen[match(DE_gen,rownames(info_gene))]<-1
info_gene$final_gen[match(gen_final,rownames(info_gene))]<-1

## Modify non SVG parameters to make their means equal
for (i in 1:length(normal_gene)){
  indx<-which(colnames(domain_mean)==normal_gene[i])
  domain_mean[,indx]<-rep(sample(domain_mean[which(domain_mean[,indx]>0),indx],1),length(domain))
}

for(i in 1:length(domain)){
  indx<-match(gener_gene,names(simSRT2@EstParam[[domain[i]]]@listData[["gene_sel1"]]))
  simSRT2@EstParam[[domain[i]]]@listData[["marginal_param1"]]@listData[["mu"]][indx]<-domain_mean[i,match(gener_gene,colnames(domain_mean))]
}



## Generate synthetic data with estimated parameters------------------------
simSRTe <- srtsim_count(simSRT2)
refcounts_all<-counts
counts<-simSRTe@simCounts[match(gen_final,rownames(simSRTe@simCounts)),]
refcounts<-simSRTe@refCounts[match(gen_final,rownames(simSRTe@refCounts)),]
counts<-as.matrix(counts)
counts<-t(counts)
kind<-info$kind
info<-data.frame(x=info$x,y=info$y)

rownames(info)<-rownames(counts)
save(counts,DE_gen,info,info_gene,kind,file="simcounts_mh.RData")
write.csv(counts,file = "simcounts_mh.csv")
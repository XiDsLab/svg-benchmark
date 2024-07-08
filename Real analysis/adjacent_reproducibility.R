setwd("/Users/xwchen/Desktop/spatial analysis/similarity_generank")

for (k in 1:15) {
  slideindex<-c(151507:151510,151669:151676,
                "01A","01B","09A","09B",
                "17A","17B","19A","19B",
                "29A","29B","36A","36B",
                "F1","F2","F2","F3")
  topN=2000
  methods<-c("Binspect","spark","meringue","spatialDE","SOMDE","sparkX","sepal","scGCO","RV","dCor","HSIC")
  env1 <- new.env()
  env2 <- new.env()
  load(paste0(slideindex[(2*k-1)],"_generank.RData"), envir = env1)
  load(paste0(slideindex[(2*k)],"_generank.RData"), envir = env2)
  Information<-matrix(data=NA,nrow=11,ncol=2)
  rownames(Information)<-methods
  colnames(Information)<-c("num_initial","Jaccard_similarity")
  for (nm in methods) {
    df1<-get(paste0("rank_",nm), envir = env1)
    df2<-get(paste0("rank_",nm), envir = env2)
    common<-intersect(names(df1),names(df2))
    Information[nm,1]<-length(common)
    
    df1<-df1[common]
    df2<-df2[common]
    df1<-df1[order(df1)]
    df2<-df2[order(df2)]
    df1<-names(df1[1:topN])
    df2<-names(df2[1:topN])
    Information[nm,2]<-length(intersect(df1,df2))/length(union(df1,df2))
  }
  
  
  Information
  write.csv(Information,paste0("pair",k,"_reproducibility.csv"))
  
}

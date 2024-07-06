load('D:/Data/slide-seq/final2/1/simulation_brainA.RData')

counts1<-counts
info1<-info
kind1<-kind
library(ggplot2)
df<-data.frame(info,kind)
ggplot(df)+
  geom_point(aes(x=x,y=y,color=kind),cex=0.8)

num<-c(100,250,500,1000,1500,5000,10000,15000,50000)##Number of downsampled spots
num2<-c(1,2)
for(i in 1:2){
  for(j in 1:9){
    filename<-paste0('D:/Data/slide-seq/final2/',num2[i],'/simulation_brainA.RData')
    load(filename)
    n<-num[j]
    x<-info$x
    y<-info$y
    # library(ggplot2)
    # ggplot(info)+geom_point(aes(x=x,y=y),cex=0.7)
    center<-which(info$x>3000&info$x<3001)#&info$y>-3001&info$y<(-3000))
    center<-center[5]
    dist<-sapply(1:nrow(counts),function(t){
      di<-sqrt((info$x[t]-info$x[center])^2+(info$y[t]-info$y[center])^2)
    })
    #which(dist==0)
    
    indx<-order(dist)[1:n]
    counts<-counts[indx,]
    info<-info[indx,]
    ggplot(info)+geom_point(aes(x=x,y=y))
    file<-paste0('D:/Data/xiacaiyang/',num[j],'/',num2[i])
    setwd(file)
    save(counts,DE_gen,info,info_gene,kind,file="simulation_brainA.RData")
    write.csv(counts,file = "brainA_counts.csv")
  }
}

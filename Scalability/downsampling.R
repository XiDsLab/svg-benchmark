load('../Simulation/reference/reference3.RData')
num<-c(100,250,500,1000,1500,5000,10000,15000,50000) ## Number of downsampled spots
num2<-c(1:2) ## Number of times to repeat down-sampling
for(i in 1:2){
  for(j in 1:9){
    n<-num[j]
    # library(ggplot2)
    # ggplot(info)+geom_point(aes(x=info$x,y=info$y),cex=0.7)
    
    ## Please note that implementing completely random down-sampling will completely 
    ## destroy the spatial structure of the data, causing the calculation of memory 
    ## and time to deviate from the real situation of the SVG algorithms.
    
    ## Select n spatial locations' coordinates
    center<-which(info$x>3000&info$x<3001) # &info$y>-3001&info$y<(-3000))
    center<-sample(center,1)
    dist<-sapply(1:nrow(counts),function(t){
      di<-sqrt((info$x[t]-info$x[center])^2+(info$y[t]-info$y[center])^2)
    })
    index<-order(dist)[1:n]
    counts<-counts[index,]
    info<-info[index,]
    
    ## Select 10,000 genes
    index<-sample(colnames(counts),10000)
    counts<-counts[,index]
    
    ## Save down-sampling results
    setwd(paste0('D:/Data/downsampling/',num[j],'/',num2[i]))
    save(counts,info,file="scalability_test.RData")
    write.csv(counts,file ="scalability_test.csv")
  }
}

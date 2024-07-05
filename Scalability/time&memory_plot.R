library(ggplot2)
library(reticulate)
library(ggpubr)
method1<-c('Binspect','dCor','RV','meringue','HSIC','spark','sparkX','sepal','SOMDE','spatialDE','scGCO')
method2<-c('Binspect','dCor','RV','meringue','HSIC','spark','sparkX','sepal','SOMDE','spatialDE','scGCO')
method3<-c('Binspect','dCor','HSIC','sparkX','sepal','SOMDE','scGCO')
method4<-c('Binspect','dCor','sparkX','sepal','SOMDE')
method<-c()
spotn<-c()
t_clock<-c()
t_cpu<-c()
m<-c()
seednum<-c()

pd <- import("pandas")
spotnum<-c(500,1000,2000,4000,8000,16000,32000,50000)
seed<-c(1,2)

for(i in 1:5){
  for(j in 1){
    file<-paste0('D:/Data/downsampling/',spotnum[i],'/',seed[j])
    setwd(file)
    
    load('brainA_result_Binspect.RData')
    t_clock<-c(t_clock,tm$Elapsed_Time_sec)
    t_cpu<-c(t_cpu,tm$Elapsed_Time_sec*20)
    m<-c(m,tm$Peak_RAM_Used_MiB)
    seednum<-c(seednum,seed[j])
    
    load('brainA_result_Dcor.RData')
    t_clock<-c(t_clock,tm$Elapsed_Time_sec)
    t_cpu<-c(t_cpu,tm$Elapsed_Time_sec*20)
    load('brainA_result_Dcor_mem.RData')
    m<-c(m,tm$Peak_RAM_Used_MiB)
    seednum<-c(seednum,seed[j])
    
    load('brainA_result_RVcor.RData')
    t_clock<-c(t_clock,tm$Elapsed_Time_sec)
    t_cpu<-c(t_cpu,tm$Elapsed_Time_sec*20)
    load('brainA_result_RVcor_mem.RData')
    m<-c(m,tm$Peak_RAM_Used_MiB)
    seednum<-c(seednum,seed[j])
    
    load('brainA_result_meringue.RData')
    t_clock<-c(t_clock,tm$Elapsed_Time_sec)
    t_cpu<-c(t_cpu,tm$Elapsed_Time_sec)
    m<-c(m,tm$Peak_RAM_Used_MiB)
    seednum<-c(seednum,seed[j])
    
    load('brainA_result_HSIC.RData')
    t_clock<-c(t_clock,tm$Elapsed_Time_sec)
    t_cpu<-c(t_cpu,tm$Elapsed_Time_sec*20)
    load('brainA_result_HSIC_mem.RData')
    m<-c(m,tm$Peak_RAM_Used_MiB)
    seednum<-c(seednum,seed[j])
    
    load('brainA_result_spark.RData')
    t_clock<-c(t_clock,tm$Elapsed_Time_sec*20/20)
    t_cpu<-c(t_cpu,tm$Elapsed_Time_sec*20)
    m<-c(m,tm$Peak_RAM_Used_MiB)
    seednum<-c(seednum,seed[j])
    
    load('brainA_result_sparkX.RData')
    t_clock<-c(t_clock,tm$Elapsed_Time_sec)
    t_cpu<-c(t_cpu,tm$Elapsed_Time_sec*20)
    m<-c(m,tm$Peak_RAM_Used_MiB)
    seednum<-c(seednum,seed[j])
    
    mem_r<-read.table('result.txt',sep=' ',header = FALSE)
    sepal<-pd$read_pickle('brainA_result_sepal.data')
    t_clock<-c(t_clock,as.numeric(sepal[1]))
    t_cpu<-c(t_cpu,as.numeric(sepal[1])*24)
    #m<-c(m,as.numeric(sepal[2])/(1024*1024))
    m<-c(m,as.numeric(mem_r[1]))
    seednum<-c(seednum,seed[j])
    
    SOMDE<-pd$read_pickle("brainA_result_SOMDE.data")
    t_clock<-c(t_clock,as.numeric(SOMDE[1]))
    t_cpu<-c(t_cpu,as.numeric(SOMDE[1]))
    #m<-c(m,as.numeric(SOMDE[2])/(1024*1024))
    m<-c(m,as.numeric(mem_r[2]))
    seednum<-c(seednum,seed[j])
    
    spatial <- pd$read_pickle("brainA_result_spatialDE.data")
    t_clock<-c(t_clock,as.numeric(spatial[1]))
    t_cpu<-c(t_cpu,as.numeric(spatial[1]))
    #m<-c(m,as.numeric(spatial[2])/(1024*1024))
    m<-c(m,as.numeric(mem_r[3]))
    seednum<-c(seednum,seed[j])
    
    scgo<-pd$read_pickle("brainA_results_scGCO.data")
    t_clock<-c(t_clock,as.numeric(scgo[1])*as.numeric(scgo[3])/20)
    t_cpu<-c(t_cpu,as.numeric(scgo[1])*as.numeric(scgo[3]))
    #m<-c(m,as.numeric(scgo[2])/(1024*1024))
    m<-c(m,as.numeric(mem_r[4]))
    seednum<-c(seednum,seed[j])
    
    method<-c(method,method1)
    spotn<-c(spotn,rep(spotnum[i],length(method1)))
    
  }
}


for(i in 6){
  for(j in 1){
    file<-paste0('D:/Data/downsampling/',spotnum[i],'/',seed[j])
    setwd(file)
    
    load('brainA_result_Binspect.RData')
    t_clock<-c(t_clock,tm$Elapsed_Time_sec)
    t_cpu<-c(t_cpu,tm$Elapsed_Time_sec*20)
    m<-c(m,tm$Peak_RAM_Used_MiB)
    seednum<-c(seednum,seed[j])
    
    load('brainA_result_Dcor.RData')
    t_clock<-c(t_clock,tm$Elapsed_Time_sec)
    t_cpu<-c(t_cpu,tm$Elapsed_Time_sec*20)
    load('brainA_result_Dcor_mem.RData')
    m<-c(m,tm$Peak_RAM_Used_MiB)
    seednum<-c(seednum,seed[j])
    
    load('brainA_result_RVcor.RData')
    t_clock<-c(t_clock,tm$Elapsed_Time_sec)
    t_cpu<-c(t_cpu,tm$Elapsed_Time_sec*20)
    load('brainA_result_RVcor_mem.RData')
    m<-c(m,tm$Peak_RAM_Used_MiB)
    seednum<-c(seednum,seed[j])
    
    load('brainA_result_meringue.RData')
    t_clock<-c(t_clock,tm$Elapsed_Time_sec)
    t_cpu<-c(t_cpu,tm$Elapsed_Time_sec)
    m<-c(m,tm$Peak_RAM_Used_MiB)
    seednum<-c(seednum,seed[j])
    
    load('brainA_result_HSIC.RData')
    t_clock<-c(t_clock,tm$Elapsed_Time_sec)
    t_cpu<-c(t_cpu,tm$Elapsed_Time_sec*20)
    load('brainA_result_HSIC_mem.RData')
    m<-c(m,tm$Peak_RAM_Used_MiB)
    seednum<-c(seednum,seed[j])
    
    load('brainA_result_spark.RData')
    t_clock<-c(t_clock,tm$Elapsed_Time_sec*20/20)
    t_cpu<-c(t_cpu,tm$Elapsed_Time_sec*20)
    m<-c(m,tm$Peak_RAM_Used_MiB)
    seednum<-c(seednum,seed[j])
    
    load('brainA_result_sparkX.RData')
    t_clock<-c(t_clock,tm$Elapsed_Time_sec)
    t_cpu<-c(t_cpu,tm$Elapsed_Time_sec*20)
    m<-c(m,tm$Peak_RAM_Used_MiB)
    seednum<-c(seednum,seed[j])
    
    mem_r<-read.table('result.txt',sep=' ',header = FALSE)
    sepal<-pd$read_pickle('brainA_result_sepal.data')
    t_clock<-c(t_clock,as.numeric(sepal[1]))
    t_cpu<-c(t_cpu,as.numeric(sepal[1])*24)
    #m<-c(m,as.numeric(sepal[2])/(1024*1024))
    m<-c(m,as.numeric(mem_r[1]))
    seednum<-c(seednum,seed[j])
    
    SOMDE<-pd$read_pickle("brainA_result_SOMDE.data")
    t_clock<-c(t_clock,as.numeric(SOMDE[1]))
    t_cpu<-c(t_cpu,as.numeric(SOMDE[1]))
    #m<-c(m,as.numeric(SOMDE[2])/(1024*1024))
    m<-c(m,as.numeric(mem_r[2]))
    seednum<-c(seednum,seed[j])
    
    spatial <- pd$read_pickle("brainA_result_spatialDE.data")
    t_clock<-c(t_clock,as.numeric(spatial[1]))
    t_cpu<-c(t_cpu,as.numeric(spatial[1]))
    #m<-c(m,as.numeric(spatial[2])/(1024*1024))
    m<-c(m,as.numeric(mem_r[3]))
    seednum<-c(seednum,seed[j])
    
    scgo<-pd$read_pickle("brainA_results_scGCO.data")
    t_clock<-c(t_clock,as.numeric(scgo[1])*as.numeric(scgo[3])/20)
    t_cpu<-c(t_cpu,as.numeric(scgo[1])*as.numeric(scgo[3]))
    #m<-c(m,as.numeric(scgo[2])/(1024*1024))
    m<-c(m,as.numeric(mem_r[4]))
    seednum<-c(seednum,seed[j])
    
    
    method<-c(method,method2)
    spotn<-c(spotn,rep(spotnum[i],length(method2)))
    
  }
}

for(i in 7){
  for(j in 1){
    file<-paste0('D:/Data/downsampling/',spotnum[i],'/',seed[j])
    setwd(file)
    
    load('brainA_result_Binspect.RData')
    t_clock<-c(t_clock,tm$Elapsed_Time_sec)
    t_cpu<-c(t_cpu,tm$Elapsed_Time_sec*20)
    m<-c(m,tm$Peak_RAM_Used_MiB)
    seednum<-c(seednum,seed[j])
    
    load('brainA_result_Dcor.RData')
    t_clock<-c(t_clock,tm$Elapsed_Time_sec)
    t_cpu<-c(t_cpu,tm$Elapsed_Time_sec*20)
    load('brainA_result_Dcor_mem.RData')
    m<-c(m,tm$Peak_RAM_Used_MiB)
    seednum<-c(seednum,seed[j])
    # 
    # load('brainA_result_RVcor.RData')
    # t_clock<-c(t_clock,tm$Elapsed_Time_sec)
    # t_cpu<-c(t_cpu,tm$Elapsed_Time_sec*20)
    # m<-c(m,tm$Peak_RAM_Used_MiB)
    # seednum<-c(seednum,seed[j])
    
    # load('brainA_result_meringue.RData')
    # t_clock<-c(t_clock,tm$Elapsed_Time_sec)
    # t_cpu<-c(t_cpu,tm$Elapsed_Time_sec)
    # m<-c(m,tm$Peak_RAM_Used_MiB)
    # seednum<-c(seednum,seed[j])
    
    load('brainA_result_HSIC.RData')
    t_clock<-c(t_clock,tm$Elapsed_Time_sec)
    t_cpu<-c(t_cpu,tm$Elapsed_Time_sec*20)
    load('brainA_result_HSIC_mem.RData')
    m<-c(m,tm$Peak_RAM_Used_MiB)
    seednum<-c(seednum,seed[j])
    
    # load('brainA_result_spark.RData')
    # t_clock<-c(t_clock,tm$Elapsed_Time_sec*10/20)
    # t_cpu<-c(t_cpu,tm$Elapsed_Time_sec*10)
    # m<-c(m,tm$Peak_RAM_Used_MiB)
    # seednum<-c(seednum,seed[j])
    
    load('brainA_result_sparkX.RData')
    t_clock<-c(t_clock,tm$Elapsed_Time_sec)
    t_cpu<-c(t_cpu,tm$Elapsed_Time_sec*20)
    m<-c(m,tm$Peak_RAM_Used_MiB)
    seednum<-c(seednum,seed[j])
    
    mem_r<-read.table('result.txt',sep=' ',header = FALSE)
    sepal<-pd$read_pickle('brainA_result_sepal.data')
    t_clock<-c(t_clock,as.numeric(sepal[1]))
    t_cpu<-c(t_cpu,as.numeric(sepal[1])*24)
    #m<-c(m,as.numeric(sepal[2])/(1024*1024))
    m<-c(m,as.numeric(mem_r[1]))
    seednum<-c(seednum,seed[j])
    
    SOMDE<-pd$read_pickle("brainA_result_SOMDE.data")
    t_clock<-c(t_clock,as.numeric(SOMDE[1]))
    t_cpu<-c(t_cpu,as.numeric(SOMDE[1]))
    #m<-c(m,as.numeric(SOMDE[2])/(1024*1024))
    m<-c(m,as.numeric(mem_r[2]))
    seednum<-c(seednum,seed[j])
    
    # spatial <- pd$read_pickle("brainA_result_spatialDE.data")
    # t_clock<-c(t_clock,as.numeric(spatial[1]))
    # t_cpu<-c(t_cpu,as.numeric(spatial[1]))
    # m<-c(m,as.numeric(spatial[2])/(1024*1024))
    # seednum<-c(seednum,seed[j])
    
    scgo<-pd$read_pickle("brainA_results_scGCO.data")
    t_clock<-c(t_clock,as.numeric(scgo[1])*as.numeric(scgo[3])/20)
    t_cpu<-c(t_cpu,as.numeric(scgo[1])*as.numeric(scgo[3]))
    #m<-c(m,as.numeric(scgo[2])/(1024*1024))
    m<-c(m,as.numeric(mem_r[3]))
    seednum<-c(seednum,seed[j])
    
    method<-c(method,method3)
    spotn<-c(spotn,rep(spotnum[i],length(method3)))
    
  }
}

for(i in 8){
  for(j in 1){
    file<-paste0('D:/Data/downsampling/',spotnum[i],'/',seed[j])
    setwd(file)
    
    load('brainA_result_Binspect.RData')
    t_clock<-c(t_clock,tm$Elapsed_Time_sec)
    t_cpu<-c(t_cpu,tm$Elapsed_Time_sec*20)
    m<-c(m,tm$Peak_RAM_Used_MiB)
    seednum<-c(seednum,seed[j])
    
    load('brainA_result_Dcor.RData')
    t_clock<-c(t_clock,tm$Elapsed_Time_sec/4)
    t_cpu<-c(t_cpu,tm$Elapsed_Time_sec*20/4)
    load('brainA_result_Dcor_mem.RData')
    m<-c(m,tm$Peak_RAM_Used_MiB)
    seednum<-c(seednum,seed[j])
    # 
    # load('brainA_result_RVcor.RData')
    # t_clock<-c(t_clock,tm$Elapsed_Time_sec)
    # t_cpu<-c(t_cpu,tm$Elapsed_Time_sec*20)
    # m<-c(m,tm$Peak_RAM_Used_MiB)
    # seednum<-c(seednum,seed[j])
    
    # load('brainA_result_meringue.RData')
    # t_clock<-c(t_clock,tm$Elapsed_Time_sec)
    # t_cpu<-c(t_cpu,tm$Elapsed_Time_sec)
    # m<-c(m,tm$Peak_RAM_Used_MiB)
    # seednum<-c(seednum,seed[j])
    
    # load('brainA_result_HSIC.RData')
    # t_clock<-c(t_clock,tm$Elapsed_Time_sec/2)
    # t_cpu<-c(t_cpu,tm$Elapsed_Time_sec*20/2)
    # load('brainA_result_HSIC_mem.RData')
    # m<-c(m,tm$Peak_RAM_Used_MiB)
    # seednum<-c(seednum,seed[j])
    
    # load('brainA_result_spark.RData')
    # t_clock<-c(t_clock,tm$Elapsed_Time_sec*10/20)
    # t_cpu<-c(t_cpu,tm$Elapsed_Time_sec*10)
    # m<-c(m,tm$Peak_RAM_Used_MiB)
    # seednum<-c(seednum,seed[j])
    
    load('brainA_result_sparkX.RData')
    t_clock<-c(t_clock,tm$Elapsed_Time_sec)
    t_cpu<-c(t_cpu,tm$Elapsed_Time_sec*20)
    m<-c(m,tm$Peak_RAM_Used_MiB)
    seednum<-c(seednum,seed[j])
    
    mem_r<-read.table('result.txt',sep=' ',header = FALSE)
    sepal<-pd$read_pickle('brainA_result_sepal.data')
    t_clock<-c(t_clock,as.numeric(sepal[1]))
    t_cpu<-c(t_cpu,as.numeric(sepal[1])*24)
    #m<-c(m,as.numeric(sepal[2])/(1024*1024))
    m<-c(m,as.numeric(mem_r[1]))
    seednum<-c(seednum,seed[j])
    
    SOMDE<-pd$read_pickle("brainA_result_SOMDE.data")
    t_clock<-c(t_clock,as.numeric(SOMDE[1]))
    t_cpu<-c(t_cpu,as.numeric(SOMDE[1]))
    #m<-c(m,as.numeric(SOMDE[2])/(1024*1024))
    m<-c(m,as.numeric(mem_r[2]))
    seednum<-c(seednum,seed[j])
    
    # spatial <- pd$read_pickle("brainA_result_spatialDE.data")
    # t_clock<-c(t_clock,as.numeric(spatial[1]))
    # t_cpu<-c(t_cpu,as.numeric(spatial[1]))
    # m<-c(m,as.numeric(spatial[2])/(1024*1024))
    # seednum<-c(seednum,seed[j])
    
    # scgo<-pd$read_pickle("brainA_results_scGCO.data")
    # t_clock<-c(t_clock,as.numeric(scgo[1])*as.numeric(scgo[3])/20)
    # t_cpu<-c(t_cpu,as.numeric(scgo[1])*as.numeric(scgo[3]))
    # m<-c(m,as.numeric(scgo[2])/(1024*1024))
    # seednum<-c(seednum,seed[j])
    
    method<-c(method,method4)
    spotn<-c(spotn,rep(spotnum[i],length(method4)))
    
  }
}


df<-data.frame(memory=m,method=method,time_clock=t_clock,time_cpu=t_cpu,seednum=seednum,spotnum=spotn)

df$time_clock<-df$time_clock/60
df$time_cpu<-df$time_cpu/60


df2<-data.frame(memory=(df$memory[which(df$seednum==1)]),
                method=df$method[which(df$seednum==1)],
                time_clock=(df$time_clock[which(df$seednum==1)]),
                time_cpu=(df$time_cpu[which(df$seednum==1)]),
                spotnum=df$spotnum[which(df$seednum==1)])
mycolor<-c("#4da0a0","#ffbc14","#156077","#2e409a","#942d8d",
           "#00bdcd","#9999CC","#e4ce00", "#9ec417","#CC6666","#56B4E9")
names(mycolor)<-c("Binspect","spark","meringue",
                  "spatialDE","SOMDE","sparkX","sepal",
                  "scGCO","RV","dCor","HSIC")
myshapes <- c(22,21,24,25,21,21,23,24,22,22,23)
names(myshapes)<-c("Binspect","spark","meringue",
                   "spatialDE","SOMDE","sparkX","sepal",
                   "scGCO","RV","dCor","HSIC")

df3<-df2
indx<-c(rep(c(65),each=2))
df4<-df2[indx,]
df4[2,3]<-48*60
df4[2,5]<-2*df4[2,5]
df5<-df4
df5$time_clock<-log10(df5$time_clock)
df5$spotnum<-log10(df5$spotnum)
df3$time_clock<-log10(df3$time_clock)
df3$spotnum<-log10(df3$spotnum)
df3$time_cpu<-log10(df3$time_cpu)

library(latex2exp)


df3<-df2

df3$memory<-log10(df3$memory)
df3$spotnum<-log10(df3$spotnum)

## memory plot
ggplot(df3)+
  geom_line(aes(x=spotnum,y=memory,color=method),linewidth=3)+geom_point(aes(x=spotnum,y=memory,shape=method,fill=method),size=8)+
  scale_shape_manual(values = myshapes)+
  scale_color_manual(values = mycolor)+
  scale_fill_manual(values = mycolor)+
  scale_y_continuous(breaks=c(1024,1024*5,1024*10,1024*20,1024*40,1024*60), labels = c('1G','5G','10G','20G','40G','60G'))+
  labs(x='Number of cells',y='Memory',color='Method',fill='Method',shape='Method')+
  geom_hline(yintercept=5*1024, linetype='dashed', col = 'red')+
  annotate("text", x = log10(700), y = 5*1024, label = "5 GB",size=12, vjust = -0.5)+
  geom_hline(yintercept=50.4*1024, linetype='dashed', col = 'red')+
  annotate("text", x = log10(700), y = 46*1024, label = "50.4 GB",size=12, vjust = -0.5)+
  scale_x_continuous(breaks=unique(df3$spotnum), labels = c('500','1000','2000','4000','8000','16k','32k','50k'))+
  theme_pubr()+theme(legend.position = "right")+theme(text = element_text(size=40,hjust=0.5))+
  theme(
    legend.key.size = unit(40, "pt")
  )
df4$spotnum<-log10(df4$spotnum)

## clock time plot
ggplot()+
  geom_line(data=df3,aes(x=spotnum,y=time_clock,color=method),linewidth=3)+
  geom_point(data=df3,aes(x=spotnum,y=time_clock,shape=method,fill=method),size=8)+
  scale_shape_manual(values = myshapes)+
  scale_color_manual(values = mycolor)+
  scale_fill_manual(values = mycolor)+
  geom_hline(yintercept=24*60, linetype='dashed', col = 'red',linewidth=0.5)+
  annotate("text", x = log10(700), y = 24*60, label = "1 day", size=12,vjust = -0.5)+
  geom_hline(yintercept=48*60, linetype='dashed', col = 'red')+
  annotate("text", x = log10(700), y = 48*60, label = "2 days", size=12,vjust = 1.4)+
  geom_hline(yintercept=10, linetype='dashed', col = 'red')+
  annotate("text", x = log10(700), y = 10, label = "10 mins", size=12,vjust = -0.5)+
  scale_x_continuous(breaks=unique(df3$spotnum), labels = c('500','1000','2000','4000','8000','16k','32k','50k'))+
  scale_y_continuous(breaks = c(60,60*5,60*10,60*15,60*24,60*48),labels = c('1h','5h','10h','15h','24h','48h'))+
  labs(x='Number of cells',y='Time(min)',color='Method',shape='Method',fill='Method')+
  theme_pubr()+theme(legend.position = "right")+theme(text = element_text(size=40,hjust=0.5))+theme(
    legend.key.size = unit(40, "pt")
  )


## cpu time plot
ggplot(df3)+
  geom_line(aes(x=spotnum,y=time_cpu,color=method),linewidth=3)+
  geom_point(aes(x=spotnum,y=time_cpu,shape=method,fill=method),size=8)+
  scale_shape_manual(values = myshapes)+
  scale_color_manual(values = mycolor)+
  scale_fill_manual(values = mycolor)+
  geom_hline(yintercept=24*60, linetype='dashed', col = 'red',linewidth=0.5)+
  annotate("text", x = log10(750), y = log10(24*60), label = "1 day", size=12,vjust = 0.5)+
  geom_hline(yintercept=60, linetype='dashed', col = 'red')+
  annotate("text", x = log10(750), y = 60, label = "1 hour", size=12,vjust = -0.5)+
  geom_hline(yintercept=24*60*20, linetype='dashed', col = 'red')+
  annotate("text", x = log10(750), y =24*60*20, label = "20 days", size=12,vjust = -0.5)+
  geom_hline(yintercept=24*60*10, linetype='dotted', col = 'red',linewidth=0.5)+
  annotate("text", x = log10(750), y = 24*60*10, label = "10 days",size=12, vjust = -0.5)+
  scale_x_continuous(breaks=unique(df3$spotnum), labels = c('500','1000','2000','4000','8000','16k','32k','50k'))+
  scale_y_continuous(breaks = c(60,24*60,24*60*5,24*60*10,24*60*20),labels = c('1 h','1 day','5 days','10 days','20 days'))+
  labs(x='Number of cells',y='Time(min)',color='Method',shape='Method',fill='Method')+
  theme_pubr()+theme(legend.position = "right")+theme(text = element_text(size=40,hjust=0.5))+theme(
    legend.key.size = unit(40, "pt")
  )



library(ggplot2)
library(reticulate)
# use_python('D:/software/Anaconda3')
pd<- import("pandas")
method_R<-list(c('Binspect','dCor','RV','meringue','HSIC','spark','sparkX'),
                c('Binspect','dCor','HSIC','sparkX'),
               c('Binspect','dCor','sparkX'))
method_py<-list(c('sepal','SOMDE','spatialDE','scGCO'),
                c('sepal','SOMDE','scGCO'),
                c('sepal','SOMDE'))
clock2cpu<-c(20,20,20,1,20,20,20,24,1,1,1)
names(clock2cpu)<-c('Binspect','dCor','RV','meringue','HSIC','spark','sparkX','sepal','SOMDE','spatialDE','scGCO')
df<-matrix(NA,ncol=5)
colnames(df)<-c('method','memory','time_cpu','time_clock','spotnum')

spotnum<-c(500,1000,2000,4000,8000,16000,32000,50000)
for (i in 1:length(spotnum)){
  if(i==7){
    tmp_methodR<-method_R[[2]]
    tmp_methodpy<-method_py[[2]]
  }else if(i==8){
    tmp_methodR<-method_R[[3]]
    tmp_methodpy<-method_py[[3]]
  } else{
    tmp_methodR<-method_R[[1]]
    tmp_methodpy<-method_py[[1]]
  }
  path<-paste0('D:/Data/downsampling/',spotnum[i],'/1')
  setwd(path)
  filename<-list.files()
  r_result<-filename[grep('.RData',filename)]
  r_result<-r_result[grep('result',r_result)]
  for(j in 1:length(tmp_methodR)){
    tmp_method<-tmp_methodR[j]
    tmp_path<-r_result[grep(tmp_method,r_result)]
    if(length(tmp_path)==1){
      load(tmp_path)
      a<-clock2cpu[which(names(clock2cpu)==tmp_method)]
      new_row=matrix(c(tmp_method,tm$Peak_RAM_Used_MiB,
                         a*tm$Elapsed_Time_sec,
                         tm$Elapsed_Time_sec,spotnum[i]),ncol=5)
      df<-rbind(df,new_row)
    }else{
      load(tmp_path[1])
      a<-clock2cpu[which(names(clock2cpu)==tmp_method)]
      new_row=matrix(c(tmp_method,tm$Peak_RAM_Used_MiB,
                         a*tm$Elapsed_Time_sec,
                         tm$Elapsed_Time_sec,spotnum[i]),ncol=5)
      load(tmp_path[2])
      new_row[2]<-tm$Peak_RAM_Used_MiB
      df<-rbind(df,new_row)
    }
  }
  mem_py<-read.table('result.txt',sep=' ',header = FALSE) ## Manually summarize by using .out file  
  py_result<-filename[grep('.data',filename)]
  for (j in 1:length(tmp_methodpy)){
    tmp_method<-tmp_methodpy[j]
    tmp_path<-py_result[grep(tmp_method,py_result)]
    tmp_mem<-mem_py[2,which(mem_py[1,]==tmp_method)]
    a<-clock2cpu[which(names(clock2cpu)==tmp_method)]
    tmp_clotime<- pd$read_pickle(tmp_path)[[1]]
    new_row=matrix(c(tmp_method,tmp_mem,
                       a*tmp_clotime,
                       tmp_clotime,spotnum[i]),ncol=5)
    df<-rbind(df,new_row)
  }
  
}
df<-df[-1,]
df<-as.data.frame(df)

df$method<-as.character(df$method)
df$memory<-as.numeric(df$memory)
df$time_clock<-as.numeric(df$time_clock)
df$time_cpu<-as.numeric(df$time_cpu)
df$spotnum<-as.numeric(df$spotnum)
df$time_clock<-df$time_clock/60
df$time_cpu<-df$time_cpu/60
df$time_cpu[75]<-df$time_cpu[75]/4
df$time_clock[75]<-df$time_clock[75]/4


df$spotnum<-log10(df$spotnum)
library(ggpubr)
mycolor<-c("#4da0a0","#ffbc14","#156077","#2e409a","#942d8d",
           "#00bdcd","#9999CC","#e4ce00", "#9ec417","#CC6666","#56B4E9")
names(mycolor)<-c("Binspect","spark","meringue",
                  "spatialDE","SOMDE","sparkX","sepal",
                  "scGCO","RV","dCor","HSIC")
myshapes <- c(22,21,24,25,21,21,23,24,22,22,23)
names(myshapes)<-c("Binspect","spark","meringue",
                   "spatialDE","SOMDE","sparkX","sepal",
                   "scGCO","RV","dCor","HSIC")

library(latex2exp)
## memory plot
ggplot(df)+
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
  scale_x_continuous(breaks=unique(df$spotnum), labels = c('500','1000','2000','4000','8000','16k','32k','50k'))+
  theme_pubr()+theme(legend.position = "right")+theme(text = element_text(size=40,hjust=0.5))+
  theme(
    legend.key.size = unit(40, "pt")
  )


##clocktime plot
ggplot()+
  geom_line(data=df,aes(x=spotnum,y=time_clock,color=method),linewidth=3)+
  geom_point(data=df,aes(x=spotnum,y=time_clock,shape=method,fill=method),size=8)+
  #geom_line(data=df4,aes(x=spotnum,y=time_clock,color=method),linewidth=2,linetype='dotted')+
  #geom_point(data=df4,aes(x=spotnum,y=time_clock,shape=method,fill=method),size=8)+
  scale_shape_manual(values = myshapes)+
  scale_color_manual(values = mycolor)+
  scale_fill_manual(values = mycolor)+
  geom_hline(yintercept=24*60, linetype='dashed', col = 'red',linewidth=0.5)+
  annotate("text", x = log10(700), y = 24*60, label = "1 day", size=12,vjust = -0.5)+
  #geom_hline(yintercept=log10(60), linetype='dashed', col = 'red')+
  #annotate("text", x = log10(700), y = 60, label = "1 hour", size=12,vjust = -1)+
  geom_hline(yintercept=48*60, linetype='dashed', col = 'red')+
  annotate("text", x = log10(700), y = 48*60, label = "2 days", size=12,vjust = 1.4)+
  geom_hline(yintercept=10, linetype='dashed', col = 'red')+
  annotate("text", x = log10(700), y = 10, label = "10 mins", size=12,vjust = -0.5)+
  scale_x_continuous(breaks=unique(df$spotnum), labels = c('500','1000','2000','4000','8000','16k','32k','50k'))+
  scale_y_continuous(breaks = c(60,60*5,60*10,60*15,60*24,60*48),labels = c('1h','5h','10h','15h','24h','48h'))+
  labs(x='Number of cells',y='Time(min)',color='Method',shape='Method',fill='Method')+
  theme_pubr()+theme(legend.position = "right")+theme(text = element_text(size=40,hjust=0.5))+theme(
    legend.key.size = unit(40, "pt")
  )


##cputime plot
ggplot(df)+
  geom_line(aes(x=spotnum,y=time_cpu,color=method),linewidth=3)+
  geom_point(aes(x=spotnum,y=time_cpu,shape=method,fill=method),size=8)+
  scale_shape_manual(values = myshapes)+
  scale_color_manual(values = mycolor)+
  scale_fill_manual(values = mycolor)+
  geom_hline(yintercept=24*60, linetype='dashed', col = 'red',linewidth=0.5)+
  annotate("text", x = log10(750), y = log10(24*60), label = "1 day", size=12,vjust = 0.5)+
  geom_hline(yintercept=60, linetype='dashed', col = 'red')+
  annotate("text", x = log10(750), y = 60, label = "1 hour", size=12,vjust = -0.5)+
  #geom_hline(yintercept=1, linetype='dashed', col = 'red')+
  #annotate("text", x = log10(32000), y = 1, label = "1 min", size=12,vjust = -0.5)+
  geom_hline(yintercept=24*60*20, linetype='dashed', col = 'red')+
  annotate("text", x = log10(750), y =24*60*20, label = "20 days", size=12,vjust = -0.5)+
  geom_hline(yintercept=24*60*10, linetype='dotted', col = 'red',linewidth=0.5)+
  annotate("text", x = log10(750), y = 24*60*10, label = "10 days",size=12, vjust = -0.5)+
  scale_x_continuous(breaks=unique(df$spotnum), labels = c('500','1000','2000','4000','8000','16k','32k','50k'))+
  scale_y_continuous(breaks = c(60,24*60,24*60*5,24*60*10,24*60*20),labels = c('1 h','1 day','5 days','10 days','20 days'))+
  labs(x='Number of cells',y='Time(min)',color='Method',shape='Method',fill='Method')+
  theme_pubr()+theme(legend.position = "right")+theme(text = element_text(size=40,hjust=0.5))+theme(
    legend.key.size = unit(40, "pt")
  )




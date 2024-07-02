tech<-c('10X','stereoseq','slideseq')
seed<-1:10
data<-as.data.frame(matrix(nrow=0,ncol=13))
for (i in 1:3){
  for(j in 1:10){
    filename<-paste0('D:/Data/summary/',tech[i],'/',seed[j])
    setwd(filename)
    
    #filename<-paste0('D:/Data/summary/',tech,'/',seed)
    #setwd(filename)
    data1<-read.csv('brainA_sim.csv')
    data1<-data.frame(data1,tech=rep(tech[i],nrow(data1)),seed=rep(seed[j],nrow(data1)))
    data<-rbind(data,data1)
  }
}
library(ggplot2)
library(ggpubr)
col_fun1 = colorRamp2(c(0, 0.5, 1), c("#0f86a9", "white", "#FC8452"))
col_fun2 = colorRamp2(c(0, 0.5, 1), c("#A5CC26", "white", "#FF7BAC"))
col_fun3 = colorRamp2(c(0, 0.5, 1), c("#3FA9F5", "white", "#FF931E"))
col_fun4 = colorRamp2(c(0, 0.5, 1), c("#ffa500", "white", "#B3A9EB"))
data$AUPR[which(is.na(data$AUPR))]<-0
data$AUROC[which(is.na(data$AUROC))]<-0
data$EP[which(is.na(data$EP))]<-0

## boxplot
df<-data.frame(AUPR=data$AUPR,method=data$X,AUROC=data$AUROC,EP=data$EP,seed=data$seed,tech=data$tech,args=(data$AUPR+data$AUROC+data$EP)/3)  
indx<-which(df$EP==0)
df<-df[-indx,]
df$method<-factor(df$method,levels = c('Binspect','spatialDE','dCor','sparkX','HSIC','spark','meringue','RV','SOMDE','scGCO','sepal'))
df$tech[which(df$tech=='10X')]<-'Set1'
df$tech[which(df$tech=='stereoseq')]<-'Set2'
df$tech[which(df$tech=='slideseq')]<-'Set3'
df$tech<-factor(df$tech,levels = c('Set1','Set2','Set3'))
df2<-reshape2::melt(df,id.vars=c('tech','seed','method','AUPR',"AUROC","EP"))


ggplot(df2)+
  geom_boxplot(aes(x=method,y=value,color=method),cex=0.8)+
  geom_point(aes(x=method,y=value))+geom_jitter(aes(x=method,y=value),shape=16, position=position_jitter(0.2),cex = 0.8)+
  scale_color_manual(values = mycolor)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  labs(title = "Simulation",x="",y="",color='Method')+theme_pubr()+
  theme(legend.position = "right",,plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90,size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size=18))+
  facet_wrap(.~tech,nrow=1)

ggplot(df)+
  geom_boxplot(aes(x=method,y=arg,color=method))+geom_jitter(aes(x=method,y=arg),alpha = 0.7)+
  facet_grid(tech~.,scales = 'free')+scale_color_manual(values = mycolor)+
  theme_bw()+theme(legend.position = "right")+labs(y='')+
  theme(panel.grid=element_blank())+theme(text = element_text(size=40,hjust=0.5),
                                          axis.text.x = element_text(angle = 30,vjust = 0.85,hjust = 0.75))+theme(
    legend.key.size = unit(40, "pt")
  )

  
  
  ggplot(df)+
    geom_boxplot(aes(x=method,y=arg,color=tech))+#geom_jitter(aes(x=method,y=arg),alpha = 0.7)+
    theme_bw()+theme(legend.position = "right")+labs(y='',color='')+
    theme(panel.grid=element_blank())+theme(text = element_text(size=40,hjust=0.5),
                                            axis.text.x = element_text(angle = 30,vjust = 0.85,hjust = 0.75))+theme(
                                              legend.key.size = unit(40, "pt")
                                            )  
  
  ggplot(df)+
    geom_boxplot(aes(x=tech,y=arg,color=tech))+#geom_jitter(aes(x=method,y=arg),alpha = 0.7)+
    facet_wrap(~method,strip.position = "bottom",ncol=11)+
    theme_bw()+theme(legend.position = "right")+labs(y='',color='',x='')+
    theme(panel.grid=element_blank())+theme(text = element_text(size=40,hjust=0.5),
                                            axis.text.x = element_blank())+theme(
                                              legend.key.size = unit(40, "pt")
                                            )+theme(
                                              strip.text.x = element_text(
                                                size = 20, color = "black"                                              ) # 这里设置x轴方向的字体类型，

                                            )    
## heatmap
df<-data.frame(method=aggregate(data$AUPR,by=list(data$X,data$tech),mean)[,1],
               tech=aggregate(data$AUPR,by=list(data$X,data$tech),mean)[,2],
               AUPR=aggregate(data$AUPR,by=list(data$X,data$tech),mean)[,3],
               AUROC=aggregate(data$AUROC,by=list(data$X,data$tech),mean)[,3],
               EP=aggregate(data$EP,by=list(data$X,data$tech),mean)[,3])

X<-(df$AUPR+df$AUROC+df$EP)/3
df<-data.frame(df,X=X)
df<-df[,c(1,2,6)]
df1<-reshape2::melt(df,id.vars=c('method','tech'))
myfun<-function(x){
  return((x-min(x))/(max(x)-min(x)))
}
a<-aggregate(df1$value,by=list(df1$tech,df1$variable),myfun)
#a1<-myfun(df1$value[which(df1$tech=='10X')])
b<-a[,3]
s<-matrix(t(b),ncol=1,nrow=33)
df2<-data.frame(df1,scal_value=s)

df1<-reshape2::melt(df)
tmp1<-(order(df$AUPR+df$AUROC+df$EP))

tmp1<-(order((df2$value[which(df2$tech=='10X'&df2$variable=='X')])+
               (df2$value[which(df2$tech=='stereoseq'&df2$variable=='X')])))
tmp1<-(order((df2$value[which(df2$tech=='10X'&df2$variable=='AUPR')])+
               (df2$value[which(df2$tech=='10X'&df2$variable=='EP')])+
               (df2$value[which(df2$tech=='10X'&df2$variable=='AUROC')])+
               (df2$value[which(df2$tech=='stereoseq'&df2$variable=='AUPR')])+
               (df2$value[which(df2$tech=='stereoseq'&df2$variable=='EP')])+
               (df2$value[which(df2$tech=='stereoseq'&df2$variable=='AUROC')])))
df2$method<-factor(df2$method,levels = df2$method[tmp1])

df2$scal_value[which(df2$tech=='slideseq'&df2$scal_value==0)]<--1
df2$tech[which(df$tech=='10X')]<-'Setting1'
df2$tech[which(df$tech=='stereoseq')]<-'Setting2'
df2$tech[which(df$tech=='slideseq')]<-'Setting3'
df2$tech<-factor(df2$tech,levels = c('Setting1','Setting2','Setting3'))
ggplot(df2, aes(x = variable, y = method, fill = scal_value)) +
  facet_grid(.~tech)+
  geom_tile(color = "white", size=1) +
  scale_fill_gradientn(limits = c(-1,1),colours = c('gray','#0f86a9',"#0f86a9", "white", "#FC8452"),breaks=c(1,0,-1),
                       labels=c("high","low","NULL"))+
  coord_fixed(rati=0.3) + 
  theme_minimal() +
  labs(title = "",x="",y="")+
  theme(legend.position = "right",
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size=18))+
  scale_x_discrete(position = "top")+geom_text(aes(label = round(value,2)),data=df2[which(df2$value>0.1),],color = 'black', size = 4)+
  labs(fill='')+guides(fill=guide_colorbar(label=TRUE))+theme(strip.text.x = element_text(
    size = 30
  ))#+scale_color_continuous(limits=c(-1,-0.9))




data<-as.data.frame(matrix(nrow=0,ncol=13))
for (i in 3){
  for(j in 1){
    filename<-paste0('D:/Data/summary/',tech[i],'/',seed[j])
    setwd(filename)
    data1<-read.csv('brainA_sim.csv')
    data1<-data.frame(data1,tech=rep(tech[i],nrow(data1)),seed=rep(seed[j],nrow(data1)))
    data<-rbind(data,data1)
  }
}
data$AUPR[which(is.na(data$AUPR))]<-0
data$AUROC[which(is.na(data$AUROC))]<-0
data$EP[which(is.na(data$EP))]<-0
library(ggplot2)
df<-data.frame(method=aggregate(data$AUPR,by=list(type=data$X),mean)[,1],
               AUPR=aggregate(data$AUPR,by=list(type=data$X),mean)[,2],
               AUROC=aggregate(data$AUROC,by=list(type=data$X),mean)[,2],
               EP=aggregate(data$EP,by=list(type=data$X),mean)[,2])

df1<-reshape2::melt(df)
tmp1<-(order(df$AUPR+df$AUROC+df$EP))
df1$method<-factor(df1$method,levels = df1$method[tmp1])
p2<-ggplot(df1, aes(x = variable, y = method, fill = value)) +
  geom_tile(color = "white", size=1) +
  scale_fill_gradientn(limits = c(0.13,0.95),colours = c("#0f86a9", "white", "#FC8452"))+
  coord_fixed(rati=0.3) + 
  theme_minimal() +
  labs(title = "",x="Large-Scale",y="")+
  theme(legend.position = "right",
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size=18))+
  scale_x_discrete(position = "top")+geom_text(aes(label = round(value,2)),data=df1[which(df1$value>0.1),],color = 'black', size = 4)

p1+geom_text(aes(x=variable[which(value>0)],y=value[which(value>0)],label = round(value[which(value>0)],2)), color = 'gray', size = 4)

library(ggpubr)
library(cowplot)
library(ggforce)
h<-ggdraw() +
  draw_plot(p1, x = 0, y = 0, width = .5, height = 1) + 
  draw_plot(p2, x = 0.5, y = 0, width = .5, height = 1) 


coord_cartesian(ylim = c(0.95, 1))
df<-reshape2::melt(data)
df<-df[which(df$variable=='AUPR'|df$variable=='AUROC'|df$variable=='EP'),]
group<-paste0(df$variable,'x',df$X)
df$group<-factor(df$group,levels=c("AUPRxBinspect","AUROCxBinspect","EPxBinspect",
                                   "AUPRxspatialDE","AUROCxspatialDE","EPxspatialDE",
                                   "AUPRxdCor","AUROCxdCor","EPxdCor",
                                   "AUPRxsparkX","AUROCxsparkX","EPxsparkX",
                                   "AUPRxHSIC","AUROCxHSIC","EPxHSIC",
                                   "AUPRxspark","AUROCxspark","EPxspark",
                                   "AUPRxmeringue","AUROCxmeringue","EPxmeringue",
                                   "AUPRxRV","AUROCxRV","EPxRV",
                                   "AUPRxSOMDE","AUROCxSOMDE","EPxSOMDE",
                                   "AUPRxscGCO","AUROCxscGCO","EPxscGCO",
                                   "AUPRxsepal","AUROCxsepal","EPxsepal"
                                   ))
df<-data.frame(df,group=group)
df$X<-factor(df$X,levels = c('Binspect','spatialDE','dCor','sparkX','HSIC','spark','meringue','RV','SOMDE','scGCO','sepal'))
ggplot(df,aes(x=variable,y=value))+
  geom_boxplot(aes(x=variable,y=value,color=variable),size=1)+geom_jitter(aes(color=variable),shape=15, position = position_jitter(0.2),size=3)+
  facet_grid(. ~ X,scales = "free")+theme_minimal()+labs(x='')+theme(axis.text.x = element_blank())+theme(text = element_text(size = 40))+
  geom_vline(xintercept = 14)+
  geom_vline(xintercept = 20)


ggplot(df[which(df$tech=='10X'),],aes(x=variable,y=value))+
  geom_boxplot(aes(x=variable,y=value,color=variable),size=1)+geom_jitter(aes(color=variable),shape=15, position = position_jitter(0.2),size=3)+
  facet_wrap(. ~ X,ncol=3,nrow=4,scales = "free")+theme_minimal()+labs(x='')+theme(axis.text.x = element_blank())+theme(text = element_text(size = 40))

  
  

write.csv(df,file='sim_result.csv')

##propotion of svg & FDR
tech<-c('10X','stereoseq','slideseq')
seed<-1:10
data<-as.data.frame(matrix(nrow=0,ncol=2))
for (i in 1:3){
  for(j in 1:10){
    filename<-paste0('D:/Data/summary/',tech[i],'/',seed[j])
    setwd(filename)
    
    #filename<-paste0('D:/Data/summary/',tech,'/',seed)
    #setwd(filename)
    load('result2.RData')
    data1<-data.frame(result,method=rownames(result))
    data1<-data.frame(data1,tech=rep(tech[i],nrow(data1)),seed=rep(seed[j],nrow(data1)))
    data<-rbind(data,data1)
  }
}
indx<-which(data$FDR==0)
data<-data[-indx,]
library(ggplot2)
library(ggpubr)
df<-data
df$method<-factor(df$method,levels = c('Binspect','spatialDE','dCor','sparkX','HSIC','spark','meringue','RVcor','SOMDE','scGCO','sepal'))
df$method<-factor(df$method,levels = c('Binspect','HSIC','dCor','sparkX','SOMDE','spatialDE','spark','sepal','RVcor','meringue','scGCO'))
df$method<-factor(df$method,levels = c('sepal','Binspect','HSIC','sparkX','dCor','SOMDE','spark','spatialDE','scGCO','RVcor','meringue'))

ggplot(df)+
  geom_boxplot(aes(x=method,y=svg.propotion,color=method))+
  geom_point(aes(x=method,y=svg.propotion))+geom_jitter(aes(x=method,y=svg.propotion),shape=16, position=position_jitter(0.2),cex = 0.8)+
  scale_color_manual(values = mycolor)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  labs(title = "Simulation",x="",y="",color='Method')+theme_pubr()+
  theme(legend.position = "right",,plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90,size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size=18))

indx<-which(df$method=='sepal')
df1<-df[-indx,]
ggplot(df1)+
  geom_boxplot(aes(x=method,y=FDR,color=method))+
  geom_point(aes(x=method,y=FDR))+geom_jitter(aes(x=method,y=FDR),shape=16, position=position_jitter(0.2),cex = 0.8)+
  scale_color_manual(values = mycolor)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  labs(title = "Simulation",x="",y="FDR",color='Method')+theme_pubr()+
  theme(legend.position = "right",,plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90,size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size=18))
df$FDR[which(df$method=='sparkX')]
write.csv(df,file='sim_data.csv')

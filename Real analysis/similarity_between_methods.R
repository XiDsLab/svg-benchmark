library(iCOBRA)
library(rjson)
library(SummarizedExperiment)
library(MultiAssayExperiment)
library(dplyr)
library(tidyr)

fn <- list.files(pattern = "\\.RData$")
methods<-c("Binspect","spark","meringue",
           "spatialDE","SOMDE","sparkX","sepal",
           "scGCO","RV","dCor","HSIC")

D <- array(dim = c(11,11,74))

for (rfn in fn) {   
  load(rfn)
  for (nm1 in methods) { 
    if(exists(paste0("rank_",nm1))){
      for (nm2 in methods) {
        if(exists(paste0("rank_",nm2))){
          if(length(seekindex)>2000){
            topN=2000
          }else{
            topN=200
          }
          tmp <- get(paste0("rank_",nm1))
          tmp<-tmp[order(tmp)]
          df1<-names(tmp[1:topN])
          
          tmp <- get(paste0("rank_",nm2))
          tmp<-tmp[order(tmp)]
          df2<-names(tmp[1:topN])
          
          tmp<-length(intersect(df1,df2))/length(union(df1,df2))
          D[which(methods==nm1),which(methods==nm2),which(fn==rfn)]<-tmp 
        }
      }
    }
  }
  print(which(fn==rfn))
  rm(rank_Binspect,rank_spark,rank_meringue,rank_spatialDE,rank_SOMDE,rank_sparkX,rank_sepal,rank_scGCO,rank_RV,rank_dCor,rank_HSIC)
}

res <- apply(D, 1:2, mean, na.rm = TRUE)
rownames(res)=colnames(res)=methods
require(graphics)
tmpdist <- 1 - res
hc <- hclust(as.dist(tmpdist),"single")
plot(hc)


library(psych)
library(BiocManager)
#install("ComplexHeatmap")
library(ComplexHeatmap)
library(circlize)
res<-res[c("RV","dCor","HSIC","sparkX","Binspect","meringue","SOMDE",
           "spark","spatialDE","sepal","scGCO"),
         c("RV","dCor","HSIC","sparkX","Binspect","meringue","SOMDE",
           "spark","spatialDE","sepal","scGCO")]

col_fun5 = colorRamp2(c(0, 0.5, 1), c("white", "#FF5D5D", "#FF0000"))

cellwidth = 0.7
cellheight = 0.7
cn = dim(res)[2]
rn = dim(res)[1]
res<-res[c("spark","spatialDE","SOMDE","meringue","Binspect",
           "sparkX","dCor","HSIC","RV","scGCO","sepal"),]
res<-res[,c("spark","spatialDE","SOMDE","meringue","Binspect",
            "sparkX","dCor","HSIC","RV","scGCO","sepal")]
res2<-(res-min(res))/(max(res)-min(res))

Heatmap(res2,name ="r", col = col_fun5,width = unit(cellwidth*cn, "cm"),height = unit(cellheight*rn, "cm"),
        rect_gp = gpar(col = "white", lwd = 4),cluster_rows = F,cluster_columns  = F,
        row_title = NULL,column_names_side = c("top"),column_names_gp = gpar(fontsize = 7.5),
        column_names_rot = 90,column_names_centered = TRUE,row_names_centered = TRUE,row_names_side = c("left"),
        row_names_gp = gpar(fontsize = 8),heatmap_legend_param = list(title = "",
                                    legend_height = unit(3, "cm"), 
                                    grid_width = unit(0.4, "cm"),
                                    labels_gp = gpar(col = "gray20", 
                                                     fontsize = 8)),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.2f", res[i, j]), x, y, 
                    gp = gpar(fontsize = 6))}
)

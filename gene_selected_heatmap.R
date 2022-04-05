Sys.setenv(LANGUAGE = "en") #显示英文报错信息
gc()
memory.limit(9999999999)
set.seed(123)
rm(list = ls())  
options(stringsAsFactors = F)
library(Seurat)
library(tidyverse)
library(patchwork)
library(SCENIC)
library(plyr)
library(permute)
library(data.table)
library(SCopeLoomR)
library(pheatmap)
library("biomaRt")
library("org.Hs.eg.db")
library("org.Mm.eg.db")
library("DESeq2")
library(edgeR)
library(limma)
library("Rgraphviz")
setwd("D:/anxiiuli_data/count_mouse_new")


setwd("D:/课题总结/ptp4a3转录组分析/figure1-安秀丽转录组文章")




genes<-read.table("D:/课题总结/ptp4a3转录组分析/figure3-kd转录组分析/细胞周期基因.txt")
gene<-genes$V1
gene<-c("Nfya","Nfyb","Nfyc")

fpkm<-read.csv(file="ternimalE_fpkm.csv",row.names = 1)
head(fpkm)
fpkm$PRO<-apply(fpkm[,1:3],1, mean, na.rm = T) 
fpkm$BASO<-apply(fpkm[,4:6],1, mean, na.rm = T) 
fpkm$POLY<-apply(fpkm[,7:9],1, mean, na.rm = T) 
fpkm$ORTHO<-apply(fpkm[,10:12],1, mean, na.rm = T) 

data <-fpkm[,13:16]
data<-data[which(rowSums(data) > 0),]#去掉全为零的行 情况

newdata<-data[c(gene),]
newdata<-na.omit(newdata)
write.csv(newdata,"磷酸酶的表达情况.csv")
data<-newdata
data<-data[which(rowSums(data) > 4),]#去掉全为零的行 情况
newdata<-data
library(pheatmap)
library(ggplot2)
library(colorRamps)
library(RColorBrewer)
library(viridis)
library(cowplot)

pheatmap::pheatmap(newdata,scale = "row",
                   cluster_rows = T,
                   cluster_cols = F)

pheatmap::pheatmap(newdata,scale = "row",clustering_method = "ward.D2",
                   cluster_rows = T,show_rownames = F,
                   cluster_cols = F,color = colorRampPalette(colors = c("blue","white","red"))(100))
pheatmap(
  mat   = newdata,color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
  scale = "row",
  cluster_rows = T,
  cluster_cols = F,
  border_color = NA,
  #annotation_col = annotation_col,
  #annotation_row = annotation_row,
  show_colnames = TRUE,
  show_rownames = TRUE,
  drop_levels   = TRUE,
  fontsize  = 8#,filename = "fig2h.pdf"
)


pheatmap(
  mat   = newdata,color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
  scale = "row",clustering_method = "average",
  cluster_rows = T,
  cluster_cols = F,
  border_color = NA,
  #annotation_col = annotation_col,
  #annotation_row = annotation_row,
  show_colnames = TRUE,
  show_rownames = TRUE,
  drop_levels   = TRUE,
  fontsize  = 8#,filename = "fig2h.pdf"
)



pheatmap(newdata,fontsize = 7,
         cluster_rows = T,border_color = NA,cluster_cols = F,color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))



fpkm<-read.csv(file="ternimalE_fpkm.csv",row.names = 1)
head(fpkm)
fpkm$PRO<-apply(fpkm[,1:3],1, mean, na.rm = T) 
fpkm$BASO<-apply(fpkm[,4:6],1, mean, na.rm = T) 
fpkm$POLY<-apply(fpkm[,7:9],1, mean, na.rm = T) 
fpkm$ORTHO<-apply(fpkm[,10:12],1, mean, na.rm = T) 

data <-fpkm[,13:16]
head(data)

genes<-read.csv("D:/anxiiuli_data/ptp4a3/count/KEGG_up_kk2.csv",row.names = 1)
genes<-genes[c("mmu04151"),]
genes<-genes[c("mmu04110"),]

genes<-genes$geneID
vec <- unlist(strsplit(genes,split='/'))

#删掉所有列上都重复的
newdata<-data[c(vec),]
newdata<-na.omit(newdata)
#低值为蓝色，高值为红色，中间值为白色：
#pdf("fig1g.pdf")
pheatmap(newdata,fontsize = 7,
         cluster_rows = F,border_color = NA,cluster_cols = F,color = colorRampPalette(colors = c("white","red","darkred"))(100))

pheatmap(newdata,scale = "row",fontsize = 7,
         cluster_rows = F,border_color = NA,cluster_cols = F,color = colorRampPalette(colors = c("navyblue","white","red"))(100))
pheatmap(newdata,scale = "row",fontsize = 7,
         cluster_rows = F,border_color = NA,cluster_cols = F,color = colorRampPalette(colors = c("blue","white","red"))(100))

pheatmap(newdata,fontsize = 7,
         cluster_rows = T,border_color = NA,cluster_cols = T,color = colorRampPalette(colors = c("white","red","darkred"))(100))

pheatmap(newdata,scale = "row",fontsize = 7,
         cluster_rows = T,border_color = NA,cluster_cols = T,color = colorRampPalette(colors = c("navyblue","white","red"))(100))
pheatmap(newdata,scale = "row",fontsize = 7,
         cluster_rows = F,border_color = NA,cluster_cols = F,color = colorRampPalette(colors = c("blue","white","red"))(100))




# 加载包
library(Seurat)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
pheatmap(newdata,scale = "row",fontsize = 7,
         cluster_rows = F,border_color = NA,cluster_cols = F,color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))
pheatmap(newdata,scale = "row",fontsize = 7,
         cluster_rows = T,border_color = NA,cluster_cols = T,color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))

pheatmap(newdata,scale = "row",fontsize = 7,clustering_method = "ward.D2",
         cluster_rows = T,border_color = NA,cluster_cols = F,color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))




pheatmap(newdata,scale = "row",fontsize = 6,filename = "new2.pdf",width = 10,height = 100,
         cluster_rows = F,border_color = NA,cluster_cols = F,color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))
dev.off()


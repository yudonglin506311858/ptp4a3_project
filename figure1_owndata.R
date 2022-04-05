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
setwd("D:/课题总结/ptp4a3转录组分析/figure1-安秀丽转录组文章")


cts<-read.csv(file="ternimalE_fpkm.csv",row.names = 1)


ProgCounts<-cts[,1:12]
head(ProgCounts)
ProgCounts$PRO<-apply(ProgCounts[,1:3],1, mean, na.rm = T) 
ProgCounts$BASO<-apply(ProgCounts[,4:6],1, mean, na.rm = T) 
ProgCounts$POLY<-apply(ProgCounts[,7:9],1, mean, na.rm = T) 
ProgCounts$ORTHO<-apply(ProgCounts[,10:12],1, mean, na.rm = T) 

ProgCounts <-ProgCounts[,13:16]



gene<-read.table("D:/课题总结/ptp4a3转录组分析/GO_0004725_protein tyrosine phosphatase activity.txt",header = T)
gene<-unique(gene$Symbol)

newdata<-ProgCounts[c(gene),]
newdata<-na.omit(newdata)

newdata<-newdata[which(rowSums(newdata) > 0),]#去掉全为零的行 情况
newdata<-newdata[which(rowSums(newdata) > 4),]#去掉全为零的行 情况

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
pheatmap(newdata,scale = "row",fontsize = 7,clustering_method = "average",
         cluster_rows = T,border_color = NA,cluster_cols = F,color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))

pheatmap(newdata,scale = "row",fontsize = 7,clustering_method = "centroid",
         cluster_rows = T,border_color = NA,cluster_cols = F,color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))
pheatmap(newdata,scale = "row",fontsize = 7,clustering_method = "complete",
         cluster_rows = T,border_color = NA,cluster_cols = F,color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))

pheatmap(newdata,scale = "row",fontsize = 7,clustering_method = "single",
         cluster_rows = T,border_color = NA,cluster_cols = F,color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))
pheatmap(newdata,scale = "row",fontsize = 7,clustering_method = "median",
         cluster_rows = T,border_color = NA,cluster_cols = F,color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))
pheatmap(newdata,scale = "row",fontsize = 7,clustering_method = "mcquitty",
         cluster_rows = T,border_color = NA,cluster_cols = F,color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))
pheatmap(newdata,scale = "row",fontsize = 7,clustering_method = "median",
         cluster_rows = T,border_color = NA,cluster_cols = F,color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))


pheatmap(newdata,fontsize = 7,clustering_method = "complete",
         cluster_rows = T,border_color = NA,cluster_cols = F,color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))

pheatmap(newdata,fontsize = 7,clustering_method = "complete",
         cluster_rows = T,border_color = NA,cluster_cols = F,color =colorRampPalette(colors = c("white","blue","red"))(100))




pheatmap(newdata,scale = "row",fontsize = 6,filename = "new2.pdf",width = 10,height = 100,
         cluster_rows = F,border_color = NA,cluster_cols = F,color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))
dev.off()


cts<-read.csv("D:/anxiiuli_data/count_own/fivesample/ternimalE_fpkm.csv",row.names = 1)

# colnames(cts)<-c("gene","Length","PRO_1","PRO_2","PRO_3","PRO_4","PRO_5",
#                  "BASO_1","BASO_2","BASO_3","BASO_4","BASO_5",
#                  "POLY_1","POLY_2","POLY_3","POLY_4","POLY_5",
#                  "ORTHO_1","ORTHO_2","ORTHO_3","ORTHO_4","ORTHO_5")

ProgCounts<-cts
head(ProgCounts)
ProgCounts$PRO<-apply(ProgCounts[,1:5],1, mean, na.rm = T) 
ProgCounts$BASO<-apply(ProgCounts[,6:10],1, mean, na.rm = T) 
ProgCounts$POLY<-apply(ProgCounts[,11:15],1, mean, na.rm = T) 
ProgCounts$ORTHO<-apply(ProgCounts[,16:20],1, mean, na.rm = T) 

ProgCounts <-ProgCounts[,21:24]





cts<-read.csv(file="ternimalE_fpkm.csv",row.names = 1)
ptp4a3<-cts[c("Ptp4a3"),]
ptp4a3[c(3:7)]
a<-as.data.frame(t(ptp4a3[c(1:3)]))
b<-as.data.frame(t(ptp4a3[c(4:6)]))
c<-as.data.frame(t(ptp4a3[c(7:9)]))
d<-as.data.frame(t(ptp4a3[c(10:12)]))
a<-as.data.frame(a$Ptp4a3)
b<-as.data.frame(b$Ptp4a3)
c<-as.data.frame(c$Ptp4a3)
d<-as.data.frame(d$Ptp4a3)




data<-dplyr::bind_rows(a,b)
data<-dplyr::bind_rows(data,c)
data<-dplyr::bind_rows(data,d)

head(data)
colnames(data)<-c("pro","baso","poly","ortho")
#data<-log2(data)
library(ggplot2)
library(ggsignif)
data %>% 
  gather(key="Type", value="Expression") %>%
  ggplot( aes(x=Type, y=Expression, fill=Type)) +
  geom_violin()+
  theme(legend.position = 'none')+ geom_boxplot(width=0.4,outlier.shape = NA,fill = "white", colour=c( "black"))+ theme_bw()+ggtitle("G2M.Score in baso&poly&ortho cell") +
  
  theme(plot.title = element_text(hjust = 0.1))+ #设置标题居中
  geom_signif(
    comparisons = list(c("pro","baso","poly","ortho")),
    map_signif_level = TRUE
  )

pdf("ptp4a3表达情况.pdf")
data1<-data %>% 
  gather(key="cellType", value="Expression")

data1$cellType<-factor(data1$cellType,levels = c("pro","baso","poly","ortho"))
ggplot(data1, aes(x=cellType, y=Expression, fill=cellType)) +
  geom_violin()+
  theme(legend.position = 'none')+ geom_boxplot(width=0.05,outlier.shape = NA,fill = "white", colour=c( "black"))+ theme_bw()+ggtitle("Ptp4a3 expression in pro&baso&poly&ortho cell") +
  
  theme(plot.title = element_text(hjust = 0.5))+ #设置标题居中
  geom_signif(
    comparisons = list(c("pro","baso"),c("pro","poly"),c("pro","ortho"),c("baso","poly"),c("baso","ortho",c("poly","ortho"))),
    map_signif_level =F,step_increase = 0.1,test = t.test
  )
ggplot(data1, aes(x=cellType, y=Expression, fill=cellType)) +
  geom_bar(stat="identity")+ theme_bw()+ggtitle("Ptp4a3 expression in pro&baso&poly&ortho cell")  #设置标题居中

dev.off()








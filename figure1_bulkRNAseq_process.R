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




#读入数据并合并
#安秀丽第一批wt数据：人类
WT_pro_1<-read.table("mm_proerythroblast_1.count",header=T)
WT_baso_1<-read.table("mm_basophilic_1.count",header=T)
WT_poly_1<-read.table("mm_polychromatic_1.count",header=T)
WT_ortho_1<-read.table("mm_orthochromatic_1.count",header=T)

#安秀丽第二批wt数据（新测序）
WT_pro_2<-read.table("mm_proerythroblast_2.count",header=T)
WT_baso_2<-read.table("mm_basophilic_2.count",header=T)
WT_poly_2<-read.table("mm_polychromatic_2.count",header=T)
WT_ortho_2<-read.table("mm_orthochromatic_2.count",header=T)

WT_pro_3<-read.table("mm_proerythroblast_3.count",header=T)
WT_baso_3<-read.table("mm_basophilic_3.count",header=T)
WT_poly_3<-read.table("mm_polychromatic_3.count",header=T)
WT_ortho_3<-read.table("mm_orthochromatic_3.count",header=T)



#合并pro的数据
pro<-merge(WT_pro_1,WT_pro_2,by=c("Geneid","Chr","Start","End","Strand","Length"))
pro<-merge(pro,WT_pro_3,by=c("Geneid","Chr","Start","End","Strand","Length"))
#合并baso的数据
baso<-merge(WT_baso_1,WT_baso_2,by=c("Geneid","Chr","Start","End","Strand","Length"))
baso<-merge(baso,WT_baso_3,by=c("Geneid","Chr","Start","End","Strand","Length"))
#合并poly的数据
poly<-merge(WT_poly_1,WT_poly_2,by=c("Geneid","Chr","Start","End","Strand","Length"))
poly<-merge(poly,WT_poly_3,by=c("Geneid","Chr","Start","End","Strand","Length"))
#合并ortho的数据
ortho<-merge(WT_ortho_1,WT_ortho_2,by=c("Geneid","Chr","Start","End","Strand","Length"))
ortho<-merge(ortho,WT_ortho_3,by=c("Geneid","Chr","Start","End","Strand","Length"))
##合并五个时期的数据
pro_baso<-merge(pro,baso,by=c("Geneid","Chr","Start","End","Strand","Length"))
poly_ortho<-merge(poly,ortho,by=c("Geneid","Chr","Start","End","Strand","Length"))
pro_baso_poly_ortho<-merge(pro_baso,poly_ortho,by=c("Geneid","Chr","Start","End","Strand","Length"))


WT_count<-pro_baso_poly_ortho

head(WT_count)
dim(WT_count)
colnames(WT_count)


#删除ensmusg的版本号

#WT_count$Geneid<-gsub("\\.*","",WT_count$Geneid)
#WT_count$Geneid <- gsub("\\.[0-9]*$", "", WT_count$Geneid)
rownames(WT_count)<-WT_count$Geneid
#https://www.nhooo.com/note/qa02b6.html
#https://www.biostars.org/p/178726/
WT_count_all<-WT_count[,c(1,6:18)]
head(WT_count_all)
# colnames(WT_count_all)<-c("Geneid","Length","mm_proerythroblast_1","mm_proerythroblast_2","mm_proerythroblast_3",
#                           "mm_basophilic_1","mm_basophilic_2","mm_basophilic_3",
#                           "mm_polychromatic_1","mm_polychromatic_2","mm_polychromatic_3",
#                           "mm_orthochromatic_1","mm_orthochromatic_2","mm_orthochromatic_3")
colnames(WT_count_all)<-c("Geneid","Length","PRO_1","PRO_2","PRO_3",
                          "BASO_1","BASO_2","BASO_3",
                          "POLY_1","POLY_2","POLY_3",
                          "ORTHO_1","ORTHO_2","ORTHO_3")


cts<-WT_count_all
head(cts)
dim(cts)

write.csv(cts,"ensembl_gene_expression_count_anxiuli.csv")



#基因ID转换
library('biomaRt')
library("curl")
library(ensembldb)
library(dplyr)
library(AnnotationHub)

mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
saveRDS(mart,"mart_mmusculus.rds")
mart<-readRDS("mart_mmusculus.rds")

gene<-read.csv("ensembl_gene_expression_count_anxiuli.csv")
gene<-as.matrix(gene$X)
head(gene)
colnames(gene)[1]<-"ensembl_gene_id"
#listAttributes(mart)
id_con<-getBM(attributes=c('ensembl_gene_id','external_gene_name',"description"),filters = 'ensembl_gene_id', values = gene, mart = mart)
head(id_con)
write.csv(id_con,"mouse_gene_ensembl_transition.csv")






library(stringr)
#cts$Geneid<-str_sub(cts$Geneid,1,str_locate(cts$Geneid,"\\.")[1]-1)
id_con<-read.csv("mouse_gene_ensembl_transition.csv")
id_con<-id_con[,c(2,3)]
colnames(cts)[1]<-"ensembl_gene_id"
head(cts)
head(id_con)
cts<-merge(id_con,cts,by=c("ensembl_gene_id"))
dim(cts)
cts<-cts[which(rowSums(cts[c(4:15)]) > 0),]
write.csv(cts,file = "ternimalE_count_genesymbol.csv",row.names = F)
cts<-cts[,-1]
cts$external_gene_name<-make.names(cts$external_gene_name, unique = TRUE)
rownames(cts)<-cts$external_gene_name
write.table(cts,"data.txt", sep = "\t")


#读取count文件：
write.csv(cts,file="ternimalE_count.csv")

#读取tpm文件：
head(cts)
kb <- cts$Length / 1000
head(kb)
countdata <- cts[,3:14]
rpk <- countdata / kb
head(rpk)
tpm <- t(t(rpk)/colSums(rpk) * 1000000)
head(tpm)
write.csv(tpm,file="ternimalE_tpm.csv")

#读取fpkm文件：
fpkm <- t(t(rpk)/colSums(countdata) * 10^6) 
head(fpkm)
write.csv(fpkm,file="ternimalE_fpkm.csv")


#PCA
condition<-factor(c("PRO","PRO","PRO",
                    "BASO","BASO","BASO",
                    "POLY","POLY","POLY",
                    "ORTHO","ORTHO","ORTHO"),
                  levels = c("PRO","BASO","POLY","ORTHO"))

head(cts)
tmp<-cts[,3:14]
head(tmp)
colData <- data.frame(row.names=colnames(cts[,3:14]), condition)
head(colData,10)
dds_all<- DESeqDataSetFromMatrix(countData = cts[,3:14],colData = colData,design= ~condition)
head(dds_all)
dds_all<- DESeq(dds_all)
vsd_all<-vst(dds_all,blind=FALSE)
head(vsd_all)
dist(t(assay(vsd_all)))
plotPCA(vsd_all,intgroup="condition")


#样本的聚类图
sampleDists <- dist(t(assay(vsd_all)))
library("RColorBrewer")
library(pheatmap)
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd_all$condition, vsd_all$type, sep="")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)



gene<-read.table("D:/anxiiuli_data/ptp4a3/磷酸化酶.txt",header = T)
gene<-unique(gene$Gene)
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




write.csv(newdata,"pi3k基因表达fpkm.csv")
write.csv(newdata,"细胞周期基因表达fpkm.csv")


#做一下差异表达分析
res_all <- results(dds_all)
head(res_all)
res_all <- res_all[order(res_all$padj),]
head(res_all)
diff_gene <- subset(res_all, padj < 0.05 & (log2FoldChange < -1|log2FoldChange > 1))
head(diff_gene)
dim(diff_gene)
head(res_all)
write.csv(diff_gene,file = "diff_gene_all.csv",row.names = T)
diff_gene_up <- subset(res_all, padj < 0.05 & (log2FoldChange > 1))
write.csv(diff_gene_up,file = "diff_gene_up.csv",row.names = T)
diff_gene_down <- subset(res_all, padj < 0.05 & (log2FoldChange < -1))
write.csv(diff_gene_down,file = "diff_gene_down.csv",row.names = T)
dim(diff_gene)
dim(diff_gene_up)
dim(diff_gene_down)


resdata <- merge(as.data.frame(res_all), as.data.frame(counts(dds_all, normalized=TRUE)),by="row.names",sort=FALSE)
write.csv(resdata,file = "DEG_all.csv",row.names = F)
head(resdata)
resdata = res_all[order(resdata$pvalue),]
summary(resdata)
table(resdata$padj<0.05)#number of true 小于0.05 的基因个数


data<-cts
data<-data[,-1]
data[1:4,1:4]#就是普通的表达矩阵
dim(data)
class(data)
group=c("PRO","BASO","POLY","ORTHO")#就是分成多少组
num=c(3,3,3,3) #这几组每组分别有多少个
num
library(DESeq2)
#dir.create("Alldiffenece")
Batch_Deseq<-function(exprSet,group,num){
  ## creat a group
  group_list= factor(rep(group,num))
  group_list
  #exprSet<-data
  colData=data.frame(row.names = colnames(exprSet),
                     group=group_list)
  
  
  for (i in 1:length(group)){
    name=unique(group)[i]
    print(name)
    colData$group<-relevel(colData$group,ref=name)
    dds=DESeq2::DESeqDataSetFromMatrix(countData = exprSet,
                                       colData = colData,
                                       design = ~group) 
    dds <- dds[ rowSums(DESeq2::counts(dds)) > 10, ]
    dds <- DESeq2::DESeq(dds)
    for (j in 2:length(DESeq2::resultsNames(dds))){
      
      resname=DESeq2::resultsNames(dds)[j]
      
      res=DESeq2::results(dds, name=resname)
      print(resname)
      res_lfc <- lfcShrink(dds, coef=j, res=res, type="apeglm")
      res_lfc
      res = res[order(res$pvalue),]
      res<-na.omit(as.data.frame(res))
      print(head(res))
      summary(res)
      
      res_lfc = res_lfc[order(res_lfc$pvalue),]
      res_lfc<-na.omit(as.data.frame(res_lfc))
      print(head(res_lfc))
      summary(res_lfc)
      
      res<-res_lfc
      #resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)),by="row.names",sort=FALSE)
      #write.csv(resdata,paste0(resname,".csv"),row.names = F)
      #resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)),by="row.names",sort=FALSE)
      resdata <- merge(as.data.frame(res_lfc), as.data.frame(counts(dds, normalized=TRUE)),by="row.names",sort=FALSE)
      write.csv(resdata,paste0(resname,".csv"),row.names = F)
      
      head(resdata)
      summary(res[order(resdata$pvalue),])
      table(resdata$padj<0.05)#number of true 小于0.05 的基因个数
      diff_gene <- subset(res, padj < 0.05 & (log2FoldChange < -1|log2FoldChange > 1))
      head(diff_gene)
      dim(diff_gene)
      diff_gene <- merge(as.data.frame(diff_gene), as.data.frame(counts(dds, normalized=TRUE)),by="row.names",sort=FALSE)
      write.csv(diff_gene,file = paste0("diff_gene_all_",resname,".csv"),row.names = F)
      diff_gene_up <- subset(res, padj < 0.05 & (log2FoldChange > 1))
      diff_gene_up <- merge(as.data.frame(diff_gene_up), as.data.frame(counts(dds, normalized=TRUE)),by="row.names",sort=FALSE)
      write.csv(diff_gene_up,file = paste0("diff_gene_up_",resname,".csv"),row.names = F)
      diff_gene_down <- subset(res, padj < 0.05 & (log2FoldChange < -1))
      diff_gene_down <- merge(as.data.frame(diff_gene_down), as.data.frame(counts(dds, normalized=TRUE)),by="row.names",sort=FALSE)
      write.csv(diff_gene_down,file = paste0("diff_gene_down_",resname,".csv"),row.names = F)
      
      #pdf(paste0(resname,".pdf"))
      library(ggplot2)
      dataset <-read.csv(paste0(resname,".csv"),header = TRUE)
      dim(dataset)
      head(dataset)
      dataset <-na.omit(dataset)
      cut_off_pvalue = 0.01
      cut_off_log2FoldChange = 1
      dataset$change = ifelse(dataset$pvalue < cut_off_pvalue & abs(dataset$log2FoldChange) >= cut_off_log2FoldChange,ifelse(dataset$log2FoldChange>cut_off_log2FoldChange ,'Up','Down'), 'Stable')
      p<- ggplot(dataset,aes(x = log2FoldChange, y = -log10(pvalue),colour=change))+labs(title= paste0(resname,"_vocano plot")) +geom_point(alpha=0.4, size=3.5) +scale_color_manual(values=c("#546de5","#d2dae2","#ff4757"))+geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +geom_hline(yintercept =-log10(cut_off_pvalue),lty=4,col="black",lwd=0.8) +labs(x="log2(fold change)", y="-log10 (p-value)")+ theme_bw()+theme(plot.title = element_text(hjust = 0.5), legend.position="right",        legend.title = element_blank() )
      
      
      #options(ggrepel.max.overlaps = Inf)
      library(ggrepel)
      dataset$change = ifelse(dataset$pvalue < cut_off_pvalue & abs(dataset$log2FoldChange) >= cut_off_log2FoldChange,ifelse(dataset$log2FoldChange> cut_off_log2FoldChange ,'Up','Down'), 'Stable')
      
      dataset$label = ifelse(dataset$pvalue <cut_off_pvalue & dataset$log2FoldChange <= -1| dataset$log2FoldChange >= 1,as.character(dataset$Row.names),"")
      
      p+geom_text_repel(data = dataset, aes(x =log2FoldChange, y =-log10(pvalue), label =label), size = 3,box.padding =unit(0.5, "lines"),point.padding = unit(0.8,"lines"),segment.color ="black",show.legend = FALSE)
      ggsave(paste0(resname,"_volcano1.pdf"))
      
      
      # grammar
      library(tidyverse)
      library(magrittr)
      library(glue)
      library(data.table)
      
      # analysis
      library(DESeq2)
      
      # graphics
      library(ggplot2)
      library(ggrepel)
      library(ggsci)
      library(scales)
      library(latex2exp)
      diffData <- fread(paste0(resname,".csv"))
      
      colnames(diffData)[1] <- "gene"
      
      diffData[is.na(padj), padj := 1][]
      diffData[, p := -log10(padj)][]
      
      
      diffData[, type := "ns"][]
      diffData[log2FoldChange > 1 & padj < 0.05, type := "up"][log2FoldChange < -1 & padj < 0.05, type := "down"][]
      
      labelGene1 <- diffData[order(p, decreasing = T)][type == "up"][1:10]
      labelGene2 <- diffData[order(p, decreasing = T)][type == "down"][1:10]
      labelGene <- rbind(labelGene1,labelGene2)
      
      #pal_nejm()(8) %>% show_col()
      typeColor <- structure(
        c(pal_nejm()(2), "gray80"),
        names = c("up", "down", "ns")
      )
      
      ggplot(diffData, aes(x = log2FoldChange, y = p)) +
        geom_point(aes(color = type, size = p), show.legend = F) +
        geom_hline(yintercept = -log10(0.05), color = "gray60", linetype = "dashed") +
        geom_vline(xintercept = 1, color = "gray60", linetype = "dashed") +
        geom_vline(xintercept = -1, color = "gray60", linetype = "dashed") +
        geom_text_repel(
          data = labelGene, aes(label = gene),
          size = 3, fontface = 3,
          nudge_x = .5, nudge_y = .5) +
        scale_radius(range = c(.1, 2)) +
        scale_color_manual(values = typeColor) +
        scale_y_continuous(expand = expansion(c(0, 0.05))) +
        labs(
          x = TeX("$log_{2}(Fold\\,Change)$"),
          y = TeX("$-log_{10}(\\textit{P}\\,value)$")) +
        theme(
          aspect.ratio = 1,
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line())+labs(title= paste0(resname,"_vocano plot"))
      ggsave(paste0(resname,"_volcano2.pdf"))
      
    }
    
  }
}
Batch_Deseq(data,group,num)

#https://www.jianshu.com/p/c8dcb3348040






dir.create("Batch_Enrichment")
setwd("./Batch_Enrichment")



Batch_Enrichment<-function(exprSet,group,num){
  ## creat a group
  group_list= factor(rep(group,num))
  group_list
  #exprSet<-data
  colData=data.frame(row.names = colnames(exprSet),
                     group=group_list)
  
  
  for (i in 1:length(group)){
    name=unique(group)[i]
    print(name)
    colData$group<-relevel(colData$group,ref=name)
    dds=DESeq2::DESeqDataSetFromMatrix(countData = exprSet,
                                       colData = colData,
                                       design = ~group) 
    dds <- dds[ rowSums(DESeq2::counts(dds)) > 10, ]
    dds <- DESeq2::DESeq(dds)
    for (j in 2:length(DESeq2::resultsNames(dds))){
      
      resname=DESeq2::resultsNames(dds)[j]
      
      res=DESeq2::results(dds, name=resname)
      print(resname)
      res_lfc <- lfcShrink(dds, coef=j, res=res, type="apeglm")
      res_lfc
      res = res[order(res$pvalue),]
      res<-na.omit(as.data.frame(res))
      print(head(res))
      summary(res)
      
      res_lfc = res_lfc[order(res_lfc$pvalue),]
      res_lfc<-na.omit(as.data.frame(res_lfc))
      print(head(res_lfc))
      summary(res_lfc)
      
      res<-res_lfc
      #resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)),by="row.names",sort=FALSE)
      #write.csv(resdata,paste0(resname,".csv"),row.names = F)
      #resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)),by="row.names",sort=FALSE)
      resdata <- merge(as.data.frame(res_lfc), as.data.frame(counts(dds, normalized=TRUE)),by="row.names",sort=FALSE)
      write.csv(resdata,paste0(resname,".csv"),row.names = F)
      
      head(resdata)
      summary(res[order(resdata$pvalue),])
      table(resdata$padj<0.05)#number of true 小于0.05 的基因个数
      diff_gene <- subset(res, padj < 0.05 & (log2FoldChange < -1|log2FoldChange > 1))
      head(diff_gene)
      dim(diff_gene)
      diff_gene <- merge(as.data.frame(diff_gene), as.data.frame(counts(dds, normalized=TRUE)),by="row.names",sort=FALSE)
      write.csv(diff_gene,file = paste0("diff_gene_all_",resname,".csv"),row.names = F)
      diff_gene_up <- subset(res, padj < 0.05 & (log2FoldChange > 1))
      diff_gene_up <- merge(as.data.frame(diff_gene_up), as.data.frame(counts(dds, normalized=TRUE)),by="row.names",sort=FALSE)
      write.csv(diff_gene_up,file = paste0("diff_gene_up_",resname,".csv"),row.names = F)
      diff_gene_down <- subset(res, padj < 0.05 & (log2FoldChange < -1))
      diff_gene_down <- merge(as.data.frame(diff_gene_down), as.data.frame(counts(dds, normalized=TRUE)),by="row.names",sort=FALSE)
      write.csv(diff_gene_down,file = paste0("diff_gene_down_",resname,".csv"),row.names = F)
      
      #pdf(paste0(resname,".pdf"))
      library(ggplot2)
      dataset <-read.csv(paste0(resname,".csv"),header = TRUE)
      dim(dataset)
      head(dataset)
      dataset <-na.omit(dataset)
      cut_off_pvalue = 0.01
      cut_off_log2FoldChange = 1
      dataset$change = ifelse(dataset$pvalue < cut_off_pvalue & abs(dataset$log2FoldChange) >= cut_off_log2FoldChange,ifelse(dataset$log2FoldChange>cut_off_log2FoldChange ,'Up','Down'), 'Stable')
      p<- ggplot(dataset,aes(x = log2FoldChange, y = -log10(pvalue),colour=change))+labs(title= paste0(resname,"_volcano plot")) +geom_point(alpha=0.4, size=3.5) +scale_color_manual(values=c("#546de5","#d2dae2","#ff4757"))+geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +geom_hline(yintercept =-log10(cut_off_pvalue),lty=4,col="black",lwd=0.8) +labs(x="log2(fold change)", y="-log10 (p-value)")+ theme_bw()+theme(plot.title = element_text(hjust = 0.5), legend.position="right",        legend.title = element_blank() )
      
      
      #options(ggrepel.max.overlaps = Inf)
      library(ggrepel)
      dataset$change = ifelse(dataset$pvalue < cut_off_pvalue & abs(dataset$log2FoldChange) >= cut_off_log2FoldChange,ifelse(dataset$log2FoldChange> cut_off_log2FoldChange ,'Up','Down'), 'Stable')
      
      dataset$label = ifelse(dataset$pvalue <cut_off_pvalue & dataset$log2FoldChange <= -1| dataset$log2FoldChange >= 1,as.character(dataset$Row.names),"")
      
      p+geom_text_repel(data = dataset, aes(x =log2FoldChange, y =-log10(pvalue), label =label), size = 3,box.padding =unit(0.5, "lines"),point.padding = unit(0.8,"lines"),segment.color ="black",show.legend = FALSE)
      ggsave(paste0(resname,"_volcano1.pdf"))
      
      
      # grammar
      library(tidyverse)
      library(magrittr)
      library(glue)
      library(data.table)
      
      # analysis
      library(DESeq2)
      
      # graphics
      library(ggplot2)
      library(ggrepel)
      library(ggsci)
      library(scales)
      library(latex2exp)
      diffData <- fread(paste0(resname,".csv"))
      
      colnames(diffData)[1] <- "gene"
      
      diffData[is.na(padj), padj := 1][]
      diffData[, p := -log10(padj)][]
      
      
      diffData[, type := "ns"][]
      diffData[log2FoldChange > 1 & padj < 0.05, type := "up"][log2FoldChange < -1 & padj < 0.05, type := "down"][]
      
      labelGene1 <- diffData[order(p, decreasing = T)][type == "up"][1:10]
      labelGene2 <- diffData[order(p, decreasing = T)][type == "down"][1:10]
      labelGene <- rbind(labelGene1,labelGene2)
      
      #pal_nejm()(8) %>% show_col()
      typeColor <- structure(
        c(pal_nejm()(2), "gray80"),
        names = c("up", "down", "ns")
      )
      
      ggplot(diffData, aes(x = log2FoldChange, y = p)) +
        geom_point(aes(color = type, size = p), show.legend = F) +
        geom_hline(yintercept = -log10(0.05), color = "gray60", linetype = "dashed") +
        geom_vline(xintercept = 1, color = "gray60", linetype = "dashed") +
        geom_vline(xintercept = -1, color = "gray60", linetype = "dashed") +
        geom_text_repel(
          data = labelGene, aes(label = gene),
          size = 3, fontface = 3,
          nudge_x = .5, nudge_y = .5) +
        scale_radius(range = c(.1, 2)) +
        scale_color_manual(values = typeColor) +
        scale_y_continuous(expand = expansion(c(0, 0.05))) +
        labs(
          x = TeX("$log_{2}(Fold\\,Change)$"),
          y = TeX("$-log_{10}(\\textit{P}\\,value)$")) +
        theme(
          aspect.ratio = 1,
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line())+labs(title= paste0(resname,"_vocano plot"))
      ggsave(paste0(resname,"_volcano2.pdf"))
      
      try({
        # 
        # #不区分上下调基因：all
        library(clusterProfiler)
        library("org.Mm.eg.db")
        library(ggplot2)
        b<-as.data.frame(diff_gene$Row.names)
        eg = bitr(b$`diff_gene$Row.names`, fromType="SYMBOL", toType="ENTREZID", OrgDb ="org.Mm.eg.db")
        gene <- eg[,2]
        head(gene)
        #分子功能(MolecularFunction)
        
        ego <- enrichGO(
          
          gene          = gene,
          
          keyType = "ENTREZID",
          
          OrgDb         = org.Mm.eg.db,
          
          ont           = "MF",
          
          pAdjustMethod = "BH",
          
          pvalueCutoff  = 0.05,
          
          qvalueCutoff  = 0.05,
          
          readable      = TRUE)
        
        barplot(ego, showCategory =50)
        ggsave(paste0("barplot_GO_MF_all_",resname,".pdf"),width=20,height=20)
        
        dotplot(ego, showCategory =50)
        ggsave(paste0("dotplot_GO_MF_all_",resname,".pdf"),width=20,height=20)
        write.csv(ego,file=paste0("GO_MF_all_",resname,".csv"))
        
        
        
        #生物过程(biologicalprocess)
        
        ego <- enrichGO(
          
          gene          = gene,
          
          keyType = "ENTREZID",
          
          OrgDb         = org.Mm.eg.db,
          
          ont           = "BP",
          
          pAdjustMethod = "BH",
          
          pvalueCutoff  = 0.05,
          
          qvalueCutoff  = 0.05,
          
          readable      = TRUE)
        
        barplot(ego, showCategory =50)
        ggsave(paste0("barplot_GO_BP_all_",resname,".pdf"),width=20,height=20)
        
        dotplot(ego, showCategory =50)
        ggsave(paste0("dotplot_GO_BP_all_",resname,".pdf"),width=20,height=20)
        write.csv(ego,file=paste0("GO_BP_all_",resname,".csv"))
        
        
        #细胞组成(cellularcomponent)
        
        ego <- enrichGO(
          
          gene          = gene,
          
          keyType = "ENTREZID",
          
          OrgDb         = org.Mm.eg.db,
          
          ont           = "CC",
          
          pAdjustMethod = "BH",
          
          pvalueCutoff  = 0.05,
          
          qvalueCutoff  = 0.05,
          
          readable      = TRUE)
        barplot(ego, showCategory =50)
        ggsave(paste0("barplot_GO_CC_all_",resname,".pdf"),width=20,height=20)
        
        dotplot(ego, showCategory =50)
        ggsave(paste0("dotplot_GO_CC_all_",resname,".pdf"),width=20,height=20)
        write.csv(ego,file=paste0("GO_CC_all_",resname,".csv"))
        
        
        
        ekegg <- enrichKEGG(
          
          gene          = gene,
          
          keyType     = "kegg",
          
          organism   = "mmu",
          
          pvalueCutoff      = 0.05,
          
          pAdjustMethod     = "BH",
          
          qvalueCutoff  = 0.05
          
        )
        
        
        barplot(ekegg, showCategory =50)
        ggsave(paste0("barplot_KEGG_all_",resname,".pdf"),width=20,height=20)
        
        dotplot(ekegg, showCategory =50)
        ggsave(paste0("dotplot_KEGG_all_",resname,".pdf"),width=20,height=20)
        write.csv(ekegg,file=paste0("KEGG_all_",resname,".csv"))
        
        
        ego <- enrichGO(gene = gene,
                        OrgDb = org.Mm.eg.db, 
                        pvalueCutoff =0.05, 
                        qvalueCutoff = 0.05,
                        ont="all",
                        readable =T)
        
        
        barplot(ego, showCategory =50)
        ggsave(paste0("barplot_GO_MFBPCC_all_",resname,".pdf"),width=20,height=20)
        barplot(ego,drop = TRUE, showCategory =20,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
        ggsave(paste0("barplot_GO_MFBPCC_all_facet_grid_",resname,".pdf"),width=20,height=20)
        dotplot(ego, showCategory =50)
        ggsave(paste0("dotplot_GO_MFBPCC_all_",resname,".pdf"),width=20,height=20)
        dotplot(ego,drop = TRUE, showCategory =20,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
        ggsave(paste0("dotplot_GO_MFBPCC_all_facet_grid_",resname,".pdf"),width=20,height=20)
        write.csv(ego,file=paste0("GO_MFBPCC_all_",resname,".csv"))
        # 
        
        
        #下调基因
        library(clusterProfiler)
        library("org.Mm.eg.db")
        library(ggplot2)
        b<-as.data.frame(diff_gene_down$Row.names)
        eg = bitr(b$`diff_gene_down$Row.names`, fromType="SYMBOL", toType="ENTREZID", OrgDb ="org.Mm.eg.db")
        gene <- eg[,2]
        head(gene)
        #分子功能(MolecularFunction)
        
        ego <- enrichGO(
          
          gene          = gene,
          
          keyType = "ENTREZID",
          
          OrgDb         = org.Mm.eg.db,
          
          ont           = "MF",
          
          pAdjustMethod = "BH",
          
          pvalueCutoff  = 0.05,
          
          qvalueCutoff  = 0.05,
          
          readable      = TRUE)
        
        barplot(ego, showCategory =50)
        ggsave(paste0("barplot_GO_MF_down_",resname,".pdf"),width=20,height=20)
        
        dotplot(ego, showCategory =50)
        ggsave(paste0("dotplot_GO_MF_down_",resname,".pdf"),width=20,height=20)
        write.csv(ego,file=paste0("GO_MF_down_",resname,".csv"))
        
        
        
        #生物过程(biologicalprocess)
        
        ego <- enrichGO(
          
          gene          = gene,
          
          keyType = "ENTREZID",
          
          OrgDb         = org.Mm.eg.db,
          
          ont           = "BP",
          
          pAdjustMethod = "BH",
          
          pvalueCutoff  = 0.05,
          
          qvalueCutoff  = 0.05,
          
          readable      = TRUE)
        
        barplot(ego, showCategory =50)
        ggsave(paste0("barplot_GO_BP_down_",resname,".pdf"),width=20,height=20)
        
        dotplot(ego, showCategory =50)
        ggsave(paste0("dotplot_GO_BP_down_",resname,".pdf"),width=20,height=20)
        write.csv(ego,file=paste0("GO_BP_down_",resname,".csv"))
        
        
        #细胞组成(cellularcomponent)
        
        ego <- enrichGO(
          
          gene          = gene,
          
          keyType = "ENTREZID",
          
          OrgDb         = org.Mm.eg.db,
          
          ont           = "CC",
          
          pAdjustMethod = "BH",
          
          pvalueCutoff  = 0.05,
          
          qvalueCutoff  = 0.05,
          
          readable      = TRUE)
        barplot(ego, showCategory =50)
        ggsave(paste0("barplot_GO_CC_down_",resname,".pdf"),width=20,height=20)
        
        dotplot(ego, showCategory =50)
        ggsave(paste0("dotplot_GO_CC_down_",resname,".pdf"),width=20,height=20)
        write.csv(ego,file=paste0("GO_CC_down_",resname,".csv"))
        
        
        
        ekegg <- enrichKEGG(
          
          gene          = gene,
          
          keyType     = "kegg",
          
          organism   = "mmu",
          
          pvalueCutoff      = 0.05,
          
          pAdjustMethod     = "BH",
          
          qvalueCutoff  = 0.05
          
        )
        
        
        barplot(ekegg, showCategory =50)
        ggsave(paste0("barplot_KEGG_down_",resname,".pdf"),width=20,height=20)
        
        dotplot(ekegg, showCategory =50)
        ggsave(paste0("dotplot_KEGG_down_",resname,".pdf"),width=20,height=20)
        write.csv(ekegg,file=paste0("KEGG_down_",resname,".csv"))
        
        
        ego <- enrichGO(gene = gene,
                        OrgDb = org.Mm.eg.db, 
                        pvalueCutoff =0.05, 
                        qvalueCutoff = 0.05,
                        ont="all",
                        readable =T)
        
        
        barplot(ego, showCategory =50)
        ggsave(paste0("barplot_GO_MFBPCC_down_",resname,".pdf"),width=20,height=20)
        barplot(ego,drop = TRUE, showCategory =20,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
        ggsave(paste0("barplot_GO_MFBPCC_down_facet_grid_",resname,".pdf"),width=20,height=20)
        dotplot(ego, showCategory =50)
        ggsave(paste0("dotplot_GO_MFBPCC_down_",resname,".pdf"),width=20,height=20)
        dotplot(ego,drop = TRUE, showCategory =20,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
        ggsave(paste0("dotplot_GO_MFBPCC_down_facet_grid_",resname,".pdf"),width=20,height=20)
        write.csv(ego,file=paste0("GO_MFBPCC_down_",resname,".csv"))
        # 
        
        #上调基因
        library(clusterProfiler)
        library("org.Mm.eg.db")
        library(ggplot2)
        b<-as.data.frame(diff_gene_up$Row.names)
        eg = bitr(b$`diff_gene_up$Row.names`, fromType="SYMBOL", toType="ENTREZID", OrgDb ="org.Mm.eg.db")
        gene <- eg[,2]
        head(gene)
        #分子功能(MolecularFunction)
        
        ego <- enrichGO(
          
          gene          = gene,
          
          keyType = "ENTREZID",
          
          OrgDb         = org.Mm.eg.db,
          
          ont           = "MF",
          
          pAdjustMethod = "BH",
          
          pvalueCutoff  = 0.05,
          
          qvalueCutoff  = 0.05,
          
          readable      = TRUE)
        
        barplot(ego, showCategory =50)
        ggsave(paste0("barplot_GO_MF_up_",resname,".pdf"),width=20,height=20)
        
        dotplot(ego, showCategory =50)
        ggsave(paste0("dotplot_GO_MF_up_",resname,".pdf"),width=20,height=20)
        write.csv(ego,file=paste0("GO_MF_up_",resname,".csv"))
        
        
        
        #生物过程(biologicalprocess)
        
        ego <- enrichGO(
          
          gene          = gene,
          
          keyType = "ENTREZID",
          
          OrgDb         = org.Mm.eg.db,
          
          ont           = "BP",
          
          pAdjustMethod = "BH",
          
          pvalueCutoff  = 0.05,
          
          qvalueCutoff  = 0.05,
          
          readable      = TRUE)
        
        barplot(ego, showCategory =50)
        ggsave(paste0("barplot_GO_BP_up_",resname,".pdf"),width=20,height=20)
        
        dotplot(ego, showCategory =50)
        ggsave(paste0("dotplot_GO_BP_up_",resname,".pdf"),width=20,height=20)
        write.csv(ego,file=paste0("GO_BP_up_",resname,".csv"))
        
        
        #细胞组成(cellularcomponent)
        
        ego <- enrichGO(
          
          gene          = gene,
          
          keyType = "ENTREZID",
          
          OrgDb         = org.Mm.eg.db,
          
          ont           = "CC",
          
          pAdjustMethod = "BH",
          
          pvalueCutoff  = 0.05,
          
          qvalueCutoff  = 0.05,
          
          readable      = TRUE)
        barplot(ego, showCategory =50)
        ggsave(paste0("barplot_GO_CC_up_",resname,".pdf"),width=20,height=20)
        
        dotplot(ego, showCategory =50)
        ggsave(paste0("dotplot_GO_CC_up_",resname,".pdf"),width=20,height=20)
        write.csv(ego,file=paste0("GO_CC_up_",resname,".csv"))
        
        
        
        ekegg <- enrichKEGG(
          
          gene          = gene,
          
          keyType     = "kegg",
          
          organism   = "mmu",
          
          pvalueCutoff      = 0.05,
          
          pAdjustMethod     = "BH",
          
          qvalueCutoff  = 0.05
          
        )
        
        
        barplot(ekegg, showCategory =50)
        ggsave(paste0("barplot_KEGG_up_",resname,".pdf"),width=20,height=20)
        
        dotplot(ekegg, showCategory =50)
        ggsave(paste0("dotplot_KEGG_up_",resname,".pdf"),width=20,height=20)
        write.csv(ekegg,file=paste0("KEGG_up_",resname,".csv"))
        
        
        ego <- enrichGO(gene = gene,
                        OrgDb = org.Mm.eg.db, 
                        pvalueCutoff =0.05, 
                        qvalueCutoff = 0.05,
                        ont="all",
                        readable =T)
        
        
        barplot(ego, showCategory =50)
        ggsave(paste0("barplot_GO_MFBPCC_up_",resname,".pdf"),width=20,height=20)
        barplot(ego,drop = TRUE, showCategory =20,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
        ggsave(paste0("barplot_GO_MFBPCC_up_facet_grid_",resname,".pdf"),width=20,height=20)
        dotplot(ego, showCategory =50)
        ggsave(paste0("dotplot_GO_MFBPCC_up_",resname,".pdf"),width=20,height=20)
        dotplot(ego,drop = TRUE, showCategory =20,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
        ggsave(paste0("dotplot_GO_MFBPCC_up_facet_grid_",resname,".pdf"),width=20,height=20)
        write.csv(ego,file=paste0("GO_MFBPCC_up_",resname,".csv"))
        
      })  
    }
    
  }
}
Batch_Enrichment(data,group,num)





data<-read.csv("ternimalE_count.csv",row.names = 1)
data<-data[c(3:14)]
data<-data[which(rowSums(data) > 0),]#去掉全为零的行 情况
head(data)
YDL <- CreateSeuratObject(counts = data,project = "ternimalE")#min.cells=10,min.genes=200,
##计算每个细胞的线粒体基因转录本数的百分比（%）,使用[[ ]] 操作符存放到metadata中，mit-开头的为线粒体基因
YDL[["percent.mt"]] <- PercentageFeatureSet(object = YDL, pattern = "^mt-")
###展示基因及线粒体百分比（这里将其进行标记并统计其分布频率，"nFeature_RNA"为基因数，"nCount_RNA"为UMI数，"percent.mt"为线粒体占比）
VlnPlot(object = YDL, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(object = YDL, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)



#对数据进行标准化
##表达量数据标准化,LogNormalize的算法：A = log( 1 + ( UMIA ÷ UMITotal ) × 10000
YDL <- NormalizeData(object = YDL, normalization.method = "LogNormalize", scale.factor = 10000)
#提取那些在细胞间变异系数较大的基因
##鉴定表达高变基因(2000个）,用于下游分析,如PCA；
YDL <- FindVariableFeatures(object = YDL, selection.method = "vst", nfeatures = 2000)
YDL <- ScaleData(object = YDL, features = rownames(YDL))
#线性降维（PCA）,默认用高变基因集,但也可通过features参数自己指定；
YDL=RunPCA(object= YDL,npcs = 2,pc.genes=VariableFeatures(object = YDL))     #PCA分析
ElbowPlot(YDL)#选择top20个PC
pcSelect=2
YDL <- FindNeighbors(object = YDL, dims = 1:pcSelect)                #计算邻接距离
##接着优化模型,resolution参数决定下游聚类分析得到的分群数,对于3K左右的细胞,设为0.4-1.2 能得到较好的结果(官方说明)；如果数据量增大,该参数也应该适当增大。
YDL <- FindClusters(object = YDL, resolution =1)                  #对细胞分组,优化标准模块化
##使用Idents（）函数可查看不同细胞的分群；
head(Idents(YDL), 5)
DimPlot(YDL,reduction = "pca",label = TRUE,pt.size = 1.5)

YDL@meta.data$celltype<-factor(c("PRO","PRO","PRO",
                                 "BASO","BASO","BASO",
                                 "POLY","POLY","POLY",
                                 "ORTHO","ORTHO","ORTHO"),
                               levels = c("PRO","BASO","POLY","ORTHO"))

rownames(YDL@meta.data)
DimPlot(YDL,reduction = "pca",label = TRUE,pt.size = 1.5,group.by = "celltype")
DimPlot(YDL,reduction = "pca",label = F,pt.size = 1.5,group.by = "celltype")


Idents(YDL) <- YDL@meta.data$celltype
YDL.markers <- FindAllMarkers(YDL, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = 'roc')
# install.packages("magrittr") # package installations are only needed the first time you use it
# install.packages("dplyr")    # alternative installation of the %>%
library(magrittr) # needs to be run evYDL time you start R and want to use %>%
library(dplyr)    # alternatively, this also loads %>%
YDL.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
##存储marker

write.csv(YDL.markers,file="allmarker_ternimalE.csv")




fpkm<-read.csv(file="ternimalE_fpkm.csv",row.names = 1)
head(fpkm)
fpkm$PRO<-apply(fpkm[,1:3],1, mean, na.rm = T) 
fpkm$BASO<-apply(fpkm[,4:6],1, mean, na.rm = T) 
fpkm$POLY<-apply(fpkm[,7:9],1, mean, na.rm = T) 
fpkm$ORTHO<-apply(fpkm[,10:12],1, mean, na.rm = T) 

data <-fpkm[,13:16]
head(data)


#删掉所有列上都重复的
newdata<-data[c(YDL.markers$gene),]
newdata<-na.omit(newdata)
#低值为蓝色，高值为红色，中间值为白色：
#pdf("fig1g.pdf")
pheatmap(newdata,fontsize = 7,
         cluster_rows = F,border_color = NA,cluster_cols = F,color = colorRampPalette(colors = c("white","red","darkred"))(100))

pheatmap(newdata,scale = "row",fontsize = 7,
         cluster_rows = F,border_color = NA,cluster_cols = F,color = colorRampPalette(colors = c("navyblue","white","red"))(100))
pheatmap(newdata,scale = "row",fontsize = 7,
         cluster_rows = F,border_color = NA,cluster_cols = F,color = colorRampPalette(colors = c("blue","white","red"))(100))
pheatmap(newdata,scale = "row",fontsize = 7,show_rownames = F,
         cluster_rows = F,border_color = NA,cluster_cols = F,color = colorRampPalette(colors = c("blue","white","red"))(100))

pheatmap(newdata,fontsize = 7,
         cluster_rows = T,border_color = NA,cluster_cols = T,color = colorRampPalette(colors = c("white","red","darkred"))(100))

pheatmap(newdata,scale = "row",fontsize = 7,
         cluster_rows = T,border_color = NA,cluster_cols = T,color = colorRampPalette(colors = c("navyblue","white","red"))(100))
pheatmap(newdata,scale = "row",fontsize = 7,
         cluster_rows = F,border_color = NA,cluster_cols = F,color = colorRampPalette(colors = c("blue","white","red"))(100))


pheatmap(newdata,scale = "row",fontsize = 7,show_rownames = F,
         cluster_rows = F,border_color = NA,cluster_cols = F,color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))



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



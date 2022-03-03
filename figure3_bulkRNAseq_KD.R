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
setwd("D:/anxiiuli_data/ptp4a3/count")




#读入数据并合并
# ctrl_1<-read.table("Ctrl-1_FRAS220009649-1r.count",header=T)
# ctrl_2<-read.table("Ctrl-2_FRAS220009650-2r.count",header=T)
# ptp4a3_1<-read.table("PRL3-1-A_FRAS220009651-2r.count",header=T)
# ptp4a3_2<-read.table("PRL3-1-B_FRAS220009652-1r.count",header=T)

ctrl_1<-read.table("PRL3-1-A_FRAS220009651-2r.count",header=T)
ctrl_2<-read.table("Ctrl-2_FRAS220009650-2r.count",header=T)
ptp4a3_1<-read.table("Ctrl-1_FRAS220009649-1r.count",header=T)
ptp4a3_2<-read.table("PRL3-1-B_FRAS220009652-1r.count",header=T)


#合并数据
data<-merge(ctrl_1,ctrl_2,by=c("Geneid","Chr","Start","End","Strand","Length"))
data<-merge(data,ptp4a3_1,by=c("Geneid","Chr","Start","End","Strand","Length"))
data<-merge(data,ptp4a3_2,by=c("Geneid","Chr","Start","End","Strand","Length"))


WT_count<-data

head(WT_count)
dim(WT_count)
colnames(WT_count)


#删除ensmusg的版本号

#WT_count$Geneid<-gsub("\\.*","",WT_count$Geneid)
#WT_count$Geneid <- gsub("\\.[0-9]*$", "", WT_count$Geneid)
rownames(WT_count)<-WT_count$Geneid
#https://www.nhooo.com/note/qa02b6.html
#https://www.biostars.org/p/178726/
WT_count_all<-WT_count[,c(1,6:10)]
head(WT_count_all)
colnames(WT_count_all)<-c("Geneid","Length","ctrl_1","ctrl_2","ptp4a3_1","ptp4a3_2")


cts<-WT_count_all
head(cts)
dim(cts)

write.csv(cts,"ensembl_gene_expression_count_ptp4a3.csv")



#基因ID转换
library('biomaRt')
library("curl")
library(ensembldb)
library(dplyr)
library(AnnotationHub)

mart <- useDataset("mmapiens_gene_ensembl", useMart("ensembl"))
saveRDS(mart,"mart_human.rds")
mart<-readRDS("mart_human.rds")

gene<-read.csv("ensembl_gene_expression_count_anxiuli.csv")
gene<-as.matrix(gene$X)
head(gene)
colnames(gene)[1]<-"ensembl_gene_id"
#listAttributes(mart)
id_con<-getBM(attributes=c('ensembl_gene_id','external_gene_name',"description"),filters = 'ensembl_gene_id', values = gene, mart = mart)
head(id_con)
write.csv(id_con,"human_gene_ensembl_transition.csv")




library(stringr)
#cts$Geneid<-str_sub(cts$Geneid,1,str_locate(cts$Geneid,"\\.")[1]-1)
id_con<-read.csv("D:/anxiiuli_data/ptp4a3/count/mouse_gene_ensembl_transition.csv")
id_con<-id_con[,c(2,3)]
colnames(cts)[1]<-"ensembl_gene_id"
head(cts)
head(id_con)
cts<-merge(id_con,cts,by=c("ensembl_gene_id"))
dim(cts)
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
countdata <- cts[,3:6]
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
condition<-factor(c("ctrl","ctrl","ptp4a3","ptp4a3"),
                  levels = c("ctrl","ptp4a3"))

head(cts)
tmp<-cts[,3:6]
head(tmp)
colData <- data.frame(row.names=colnames(cts[,3:6]), condition)
head(colData,10)
dds_all<- DESeqDataSetFromMatrix(countData = cts[,3:6],colData = colData,design= ~condition)
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
cts1<-cts
cts<-cts[,3:6]
cor(cts,method = "pearson")
# ctrl_1    ctrl_2  ptp4a3_1  ptp4a3_2
# ctrl_1   1.0000000 0.9863111 0.8094934 0.5713300
# ctrl_2   0.9863111 1.0000000 0.7706288 0.5144732
# ptp4a3_1 0.8094934 0.7706288 1.0000000 0.9099334
# ptp4a3_2 0.5713300 0.5144732 0.9099334 1.0000000
heatmap(cor(cts))

pheatmap(cor(cts))

pheatmap(cor(cts,method = "pearson"),display_numbers = T)
#pheatmap(cor(cts,method = "spearman"),display_numbers = T)
cts<-cts[which(rowSums(cts) > 0),] 

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
#diff_gene_up<-read.csv(file = "diff_gene_up.csv",row.names = 1)
diff_gene_down <- subset(res_all, padj < 0.05 & (log2FoldChange < -1))
write.csv(diff_gene_down,file = "diff_gene_down.csv",row.names = T)
head(diff_gene_down)
dim(diff_gene)
dim(diff_gene_up)
dim(diff_gene_down)

head(diff_gene_up)

resdata <- merge(as.data.frame(res_all), as.data.frame(counts(dds_all, normalized=TRUE)),by="row.names",sort=FALSE)
write.csv(resdata,file = "DEG_all.csv",row.names = F)
resdata<-read.csv(file = "DEG_all.csv",row.names = 1)
head(resdata)
resdata = res_all[order(resdata$pvalue),]
summary(resdata)
table(resdata$padj<0.05)#number of true 小于0.05 的基因个数


diff_gene<-as.data.frame(diff_gene)
term<-rownames(diff_gene)
newdata<-cts[c(term),]
pheatmap(newdata,scale ="row",border_color = NA)
pheatmap(newdata,scale ="row",border_color = NA,show_rownames = F)


newdata<-cts[c(labelGene$gene),]
pheatmap(newdata,scale ="row",border_color = NA,color = colorRampPalette(colors = c("blue","white","red"))(100))

pheatmap(newdata,scale ="row",cluster_rows = T,border_color = NA,color = colorRampPalette(colors = c("blue","white","red"))(100))

pheatmap(newdata,cluster_rows =T,scale = "row",clustering_method = "average",fontsize=5,fontsize_row=5,fontsize_col=10,color=colorRampPalette(rev(c("red","white","blue")))(102))


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
diffData <- fread("DEG_all.csv")

colnames(diffData)[1] <- "gene"

diffData[is.na(padj), padj := 1][]
diffData[, p := -log10(padj)][]


diffData[, type := "ns"][]
diffData[log2FoldChange > 1 & padj < 0.05, type := "up"][log2FoldChange < -1 & padj < 0.05, type := "down"][]

labelGene1 <- diffData[order(p, decreasing = T)][type == "up"][1:10]
labelGene2 <- diffData[order(p, decreasing = T)][type == "down"][1:10]
labelGene <- rbind(labelGene1,labelGene2)
labelGene <- rbind(labelGene,diffData[1685,])#Ptp4a3

options(ggrepel.max.overlaps = Inf)
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
    axis.line = element_line())+labs(title= "ptp4a3 KD vocano plot")
ggsave("ptp4a3 KD_volcano2.pdf")

ggsave("ptp4a3 KD_volcano.pdf",width = 6,height = 6)


resname<-c("ptp4a3 KD")



#下调基因
library(clusterProfiler)
library("org.Mm.eg.db")
library(ggplot2)
b<-as.data.frame(rownames(diff_gene_down))
eg = bitr(b$`rownames(diff_gene_down)`, fromType="SYMBOL", toType="ENTREZID", OrgDb ="org.Mm.eg.db")
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
b<-as.data.frame(rownames(diff_gene_up))
eg = bitr(b$`rownames(diff_gene_up)`, fromType="SYMBOL", toType="ENTREZID", OrgDb ="org.Mm.eg.db")
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


kk2 = setReadable(ekegg,
                  OrgDb = "org.Mm.eg.db",
                  keyType = "ENTREZID")
head(kk2@result$geneID)

write.csv(kk2,file="KEGG_up_kk2.csv")



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




b<-as.data.frame(rownames(diff_gene))
eg = bitr(b$`rownames(diff_gene)`, fromType="SYMBOL", toType="ENTREZID", OrgDb ="org.Mm.eg.db")
gene <- eg[,2]
head(eg)



info<-as.data.frame(diff_gene)
library(dplyr)
diff_gene<- add_rownames(as.data.frame(diff_gene), "SYMBOL")
diff_gene<-as.data.frame(diff_gene)
head(diff_gene)

diff_gene<-merge(eg,diff_gene,by=c("SYMBOL"))



GSEA_input <- diff_gene$log2FoldChange
names(GSEA_input) = diff_gene$ENTREZID
GSEA_input = sort(GSEA_input, decreasing = TRUE)
GSEA_KEGG <- gseKEGG(GSEA_input, organism = 'mmu', pvalueCutoff = 0.05)#GSEA富集分析
#https://zhuanlan.zhihu.com/p/377356510


#GSEA富集图:
library(clusterProfiler)
library(DOSE)
library(enrichplot)
ridgeplot(GSEA_KEGG) 
gseaplot2(GSEA_KEGG,1)
gseaplot2(GSEA_KEGG,1:30)#30是根据ridgeplot中有30个富集通路得到的

gseaplot2(GSEA_KEGG,2,pvalue_table = TRUE)#输出第212个结果

##R语言准备gsea输入文件
library(Seurat)
library(tidyverse)
library(GSEA)
dir.create("GSEA")   
dir.create("GSEA/input")
dir.create("GSEA/output")

expr <- cts
expr <- data.frame(NAME=rownames(expr), Description=rep('na', nrow(expr)), expr, stringsAsFactors=F)
write('#1.2', "GSEA/input/expr.gct", ncolumns=1)
write(c(nrow(expr),(ncol(expr)-2)), "GSEA/input/expr.gct", ncolumns=2, append=T, sep='\t')
write.table(expr, "GSEA/input/expr.gct", row.names=F, sep='\t', append=T, quote=F)
line.1 <- c((ncol(expr)-2), 2, 1)
tmp <- c("ctrl","pta4a3")
line.2 <- c("#", c("ctrl","pta4a3"))
line.3 <- c("ctrl","ctrl","pta4a3","pta4a3")
write(line.1, 'GSEA/input/group.cls', ncolumns=length(line.1), append=T, sep='\t')
write(line.2, 'GSEA/input/group.cls', ncolumns=length(line.2), append=T, sep='\t')
write(line.3, 'GSEA/input/group.cls', ncolumns=length(line.3), append=T, sep='\t')

#gs.db参数对应的是基因集数据文件，这里只对kegg收录的基因集做了分析
GSEA::GSEA("GSEA/input/expr.gct", "GSEA/input/group.cls", 
           gs.db="D:/生信学习/GEO/MsigDB/symbols/c2.all.v6.2.symbols.gmt", output.directory="GSEA/output/")







# MESSAGE -----------------------------------------------------------------
#
# author: Yulin Lyu
# email: lvyulin@pku.edu.cn
#
# require: R whatever
#
# ---

# * 1. Load packages ------------------------------------------------------


# grammar
library(tidyverse)
library(magrittr)
library(glue)
library(data.table)

# analysis
library(DESeq2)
library(org.Mm.eg.db)
library(clusterProfiler)

# graphics
library(ggplot2)
library(ggsci)
library(latex2exp)
library(patchwork)

# * 2. Load data ----------------------------------------------------------

GOfile <- list.files(".", "go.csv")

# * 3. Plot ---------------------------------------------------------------

GOdata <- map(GOfile, ~ fread(.x)) %>% set_names(str_remove(GOfile, ".GO.*"))
plotData <- GOdata %>% imap(~ {.x[pvalue < 0.05 & Count >= 5][1:10][, type := .y]})
plotData %<>% map(~ {.x[, p := -log10(pvalue)]})
head(plotData$go.csv,20)
GOplot <- imap(plotData, ~ {
  ggplot(.x, aes(x = p, y = fct_reorder(Description, p,))) +
    geom_col(aes(alpha = p), fill = "black", width = .8, show.legend = F) +
    scale_alpha_continuous(range = c(.5, 1)) +
    scale_x_continuous(expand = expansion(c(0, 0.05))) +
    labs(x = TeX("$-log_{10}(\\textit{P}\\,value)$"), y = "", title = "") +
    theme(
      aspect.ratio = 0.75,
      panel.background = element_blank(),
      panel.grid = element_blank(),
      axis.line = element_line())
})

purrr::reduce(GOplot, `+`) + plot_layout(ncol = 1)

ggsave("graphics/GO.png", width = 10, height = 4 * length(GOfile))



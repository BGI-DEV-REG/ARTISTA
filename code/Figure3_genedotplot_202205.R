rm(list=ls())
gc()
library(data.table)
library(dplyr)
library(Seurat)
library(RColorBrewer)
library(ggthemes)
library(ggplot2)
library(plyr)
library(ggpubr)
library(viridis)
library(reshape2)
library(pheatmap)
library(ggthemes)
obj.integrated.Find <- readRDS("D:/02.project/06RYfish/02brain_dev/Injury_dev_inter0307_V2/Injury_dev_inter0307_V2_SCT_rpca.rds")
Idents(obj.integrated.Find) <- "Annotation_0306"
sum2=read.table("D:/02.project/06RYfish/02brain_dev/Injury_2_5_10_merge0307/Summary_Gene_Annotaion_0707.xls",header=T,fill=TRUE, na.strings = "")
#sum2$idmouse <- paste0(sum2$Axolotl_ID,":",sum2$mm_gene)
gene=as.data.frame(obj.integrated.Find@assays$SCT@data@Dimnames[[1]])
colnames(gene)="Axolotl_ID"
gene$Axolotl_ID=as.character(gene$Axolotl_ID)
gene_merge=join_all(list(gene,sum2),by="Axolotl_ID",type = "left")
obj.integrated.Find@assays$SCT@data@Dimnames[[1]]=gene_merge$Axolotl_tanaka_annotated_gene
gene=as.data.frame(obj.integrated.Find@assays$SCT@counts@Dimnames[[1]])
colnames(gene)="Axolotl_ID"
gene$Axolotl_ID=as.character(gene$Axolotl_ID)
gene_merge=join_all(list(gene,sum2),by="Axolotl_ID",type = "left")
obj.integrated.Find@assays$SCT@counts@Dimnames[[1]]=gene_merge$Axolotl_tanaka_annotated_gene
gene=as.data.frame(obj.integrated.Find@assays$RNA@counts@Dimnames[[1]])
colnames(gene)="Axolotl_ID"
gene$Axolotl_ID=as.character(gene$Axolotl_ID)
gene_merge=join_all(list(gene,sum2),by="Axolotl_ID",type = "left")
obj.integrated.Find@assays$RNA@counts@Dimnames[[1]]=gene_merge$Axolotl_tanaka_annotated_gene
gene=as.data.frame(obj.integrated.Find@assays$RNA@data@Dimnames[[1]])
colnames(gene)="Axolotl_ID"
gene$Axolotl_ID=as.character(gene$Axolotl_ID)
gene_merge=join_all(list(gene,sum2),by="Axolotl_ID",type = "left")
obj.integrated.Find@assays$RNA@data@Dimnames[[1]]=gene_merge$Axolotl_tanaka_annotated_gene

DefaultAssay(obj.integrated.Find) <- "SCT"
sum1=read.table("D:/02.project/06RYfish/02brain_dev/Injury_dev_inter0307_V2/dotplot_gene_20220307.txt",header=T,sep = "\t")
markers.to.plot <- sum1$Axolotl_tanaka_annotated_gene
Idents(obj.integrated.Find) <- "Annotation_0306"
DotPlot(obj.integrated.Find, features = markers.to.plot , cols = c("white", "red"), dot.scale = 8) + RotatedAxis()
gene_matrix <- DotPlot(obj.integrated.Find, features = markers.to.plot , cols = c("white", "red"), dot.scale = 8) + RotatedAxis()
matrix_gene_2 <- gene_matrix$data
matrix_gene_2$name <- rownames(matrix_gene_2)
vaccines <- acast(matrix_gene_2,markers.to.plot~id,value.var = "avg.exp")
result <- pheatmap(vaccines,cluster_rows = T,cluster_cols = F,scale = "row",color = viridis(30),show_colnames=T,border_color = NA,fontsize_row = 6)
row_oder=result$tree_row$order 
rn_new <- rownames(vaccines)[row_oder] 
DotPlot(obj.integrated.Find, features = rn_new , cols = c("white", "red"), dot.scale = 8) + RotatedAxis()
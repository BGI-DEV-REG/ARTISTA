library(Seurat)
library(monocle)
library(ggplot2)
library(dplyr)

setwd('../data/Figure4_data/Figure4F_Monocle2')
obj <- readRDS('Batch1_Injury_15DPI_C_rep1_FP200000266TR_E6.rds')
Idents(obj) = 'Annotation_0306'
obj = subset(obj,idents=c('WNTEGC','MPEX'))
root = 'MPEX'
co = read.csv('/mnt/e/share_wangshuai/fengweimin/Monocle2/00.data/Inj_24_Adult_Develop_color_0305.txt',sep = '\t')
celltype =as.data.frame(unique(obj$Annotation_0306))
colnames(celltype) = 'celltype'
co = merge(co,celltype,by.x = 'order',by.y = 'celltype',sort=F)
col = co$Color
names(col) = co$order

Mono_matrix<-GetAssayData(obj,slot = "counts")
feature_ann<-data.frame(gene_id=rownames(Mono_matrix),gene_short_name=rownames(Mono_matrix))
rownames(feature_ann)<-rownames(Mono_matrix)
Mono_fd<-new("AnnotatedDataFrame", data = feature_ann)
sample_ann<-obj@meta.data
Mono_pd<-new("AnnotatedDataFrame", data =sample_ann)
Mono.cds<-newCellDataSet(Mono_matrix,phenoData =Mono_pd,featureData =Mono_fd,expressionFamily=negbinomial.size())

rm(obj)
gc()
Mono.cds <- estimateSizeFactors(Mono.cds)
Mono.cds <- estimateDispersions(Mono.cds)
disp_table <- dispersionTable(Mono.cds)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.51)
Mono.cds <- setOrderingFilter(Mono.cds, unsup_clustering_genes$gene_id)
Mono.cds <- reduceDimension(
  Mono.cds,
  max_components = 5,
  method = 'DDRTree')
Mono.cds <- orderCells(Mono.cds)
head(pData(Mono.cds))
plot_cell_trajectory(Mono.cds,cell_size = 1)
plot_cell_trajectory(Mono.cds, color_by = "Annotation_0306",cell_size = 1)+ggplot2::scale_colour_manual(values = col)

Mono.cds <- orderCells(Mono.cds, root_state = c('1'))


plot_cell_trajectory(Mono.cds, color_by = "Pseudotime",cell_size = 1)+viridis::scale_color_viridis(option="viridis")
ggsave(paste0(root,'_Pseudotime.pdf'),width = 5,height=5)#width = 9.62,height = 4.81)

plot_cell_trajectory(Mono.cds, color_by = "Annotation_0306",cell_size = 1)+ggplot2::scale_colour_manual(values = col)
ggsave(paste0(root,'_Annotation_0306.pdf'),width = 5,height = 5)#width = 9.62,height = 4.81)




library(Seurat)
library(monocle)
library(ggplot2)
library(dplyr)

################

setwd('../data/FigureS16_data/FigureS16_Monocle2')
obj <- readRDS('Batch1_Stage54_telencephalon_rep2_DP8400015649BRD6_2.rds')
Idents(obj) = 'Annotation_0305'
obj = subset(obj,idents=c('DEGC','DNBL4','Immature NPTXEX'))
root = 'Immature_NPTXEX'
co = read.csv('Inj_24_Adult_Develop_color_0305.txt',sep = '\t')
celltype =as.data.frame(unique(obj$Annotation_0305))
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
unique(Mono.cds$Annotation_0305)
rm(obj)
gc()
Mono.cds <- estimateSizeFactors(Mono.cds)
Mono.cds <- estimateDispersions(Mono.cds)
disp_table <- dispersionTable(Mono.cds)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.9)
Mono.cds <- setOrderingFilter(Mono.cds, unsup_clustering_genes$gene_id)

Mono.cds <- reduceDimension(
  Mono.cds,
  max_components = 3,
  method = 'DDRTree')
Mono.cds <- orderCells(Mono.cds)
head(pData(Mono.cds))
plot_cell_trajectory(Mono.cds,cell_size = 1)
plot_cell_trajectory(Mono.cds, color_by = "Annotation_0305",cell_size = 1)+ggplot2::scale_colour_manual(values = col)

Mono.cds <- orderCells(Mono.cds, root_state = c('3'))

plot_cell_trajectory(Mono.cds, use_color_gradient = TRUE,cell_size=0.1)

plot_cell_trajectory(Mono.cds,cell_size = 1)
ggsave(paste0(root,'_State.png'))
plot_cell_trajectory(Mono.cds, color_by = "seurat_clusters",cell_size = 1)
ggsave(paste0(root,'_seurat_clusters.png'))

plot_cell_trajectory(Mono.cds, color_by = "Pseudotime",cell_size = 1)+viridis::scale_color_viridis(option="viridis")
ggsave(paste0(root,'_Pseudotime.pdf'),width = 5,height=5)#width = 9.62,height = 4.81)

plot_cell_trajectory(Mono.cds, color_by = "Annotation_0305",cell_size = 1)+ggplot2::scale_colour_manual(values = col)
ggsave(paste0(root,'_Annotation_0305.pdf'),width = 5,height = 5)#width = 9.62,height = 4.81)


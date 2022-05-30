library(Seurat)
library(monocle)
library(ggplot2)
library(dplyr)

sewd('../data/FigureS12_data/FigureS12J_Monocle2')
obj <- readRDS('Batch1_Injury_15DPI_rep4_FP200000266TR_E4.rds')
Idents(obj) = 'Annotation_0306'
obj = subset(obj,idents=c('REAEGC','RIPC2','DPEX','MPEX'))
root = 'DPEX'
co = read.csv('Inj_24_Adult_Develop_color_0305.txt',sep = '\t')
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
Idents(obj) = 'Annotation_0306'
deg_bulk=FindAllMarkers(obj,logfc.threshold = 0.5,only.pos =T)

marker = list()
id = c('EX')
for(i in id){
  data = deg_bulk
  for(j in unique(data$cluster)){
    data1 = data[which(data$cluster==j),]
    data1 = data1[order(-data1$avg_log2FC,data1$p_val),]
    if(length(data1[,1])>69){
      data2 = data1[c(1:69),]
      marker[[paste0(i,"_",j)]] = data2}
    else{marker[[paste0(i,"_",j)]] = data1}
  }
  m_list = do.call(rbind,marker)
  print(length(unique(m_list$gene)))
  m_list1 = as.data.frame(unique(m_list$gene))
  colnames(m_list1) ='top_gene'
  #write.csv(m_list1,paste0('Promoter_celltype_top50_gene_',i,'.csv'),quote=F)

}
rm(obj)
gc()
Mono.cds <- estimateSizeFactors(Mono.cds)
Mono.cds <- estimateDispersions(Mono.cds)
disp_table <- dispersionTable(Mono.cds)
Mono.cds <- setOrderingFilter(Mono.cds, m_list1$top_gene)

plot_ordering_genes(Mono.cds)

Mono.cds <- reduceDimension(
  Mono.cds,
  max_components = 3,
  method = 'DDRTree')
Mono.cds <- orderCells(Mono.cds)
head(pData(Mono.cds))
plot_cell_trajectory(Mono.cds,cell_size = 1)
plot_cell_trajectory(Mono.cds, color_by = "Annotation_0306",cell_size = 1)+ggplot2::scale_colour_manual(values = col)
Mono.cds <- orderCells(Mono.cds, root_state = c('1'))

plot_cell_trajectory(Mono.cds, use_color_gradient = TRUE,cell_size=0.1)


plot_cell_trajectory(Mono.cds,cell_size = 1)
ggsave(paste0(root,'_State.png'))
plot_cell_trajectory(Mono.cds, color_by = "seurat_clusters",cell_size = 1)
ggsave(paste0(root,'_seurat_clusters.png'))

plot_cell_trajectory(Mono.cds, color_by = "Pseudotime",cell_size = 1)+viridis::scale_color_viridis(option="viridis")
ggsave(paste0(root,'_Pseudotime.pdf'),width = 5,height=5)#width = 9.62,height = 4.81)

plot_cell_trajectory(Mono.cds, color_by = "Annotation_0306",cell_size = 1)+ggplot2::scale_colour_manual(values = col)
ggsave(paste0(root,'_Annotation_0306.pdf'),width = 5,height = 5)#width = 9.62,height = 4.81)


library(Seurat)
library(monocle)
library(ggplot2)
library(dplyr)

setwd('../data/FigureS11_data/FigureS11B_Monocle2')
obj <- readRDS('Batch2_Injury_2DPI_rep1_SS200000147BL_D5.rds')
Idents(obj) = 'Annotation_0306'
root = 'EGC'
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
unique(Mono.cds$Annotation_0306)

Idents(obj) = 'EGC'
deg_bulk=FindAllMarkers(obj,logfc.threshold = 0.5,only.pos =T)
marker = list()
id = c('EX')
for(i in id){
  data = deg_bulk
  for(j in unique(data$cluster)){
    data1 = data[which(data$cluster==j),]
    data1 = data1[order(-data1$avg_log2FC,data1$p_val),]
    if(length(data1[,1])>12){
      data2 = data1[c(1:12),]
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

Mono.cds <- reduceDimension(
  Mono.cds,
  max_components = 10,
  method = 'DDRTree')
Mono.cds <- orderCells(Mono.cds)
plot_cell_trajectory(Mono.cds, color_by = "EGC",cell_size = 1)


head(pData(Mono.cds))
plot_cell_trajectory(Mono.cds,cell_size = 1)
plot_cell_trajectory(Mono.cds, color_by = "EGC",cell_size = 1)
plot_cell_trajectory(Mono.cds, color_by = "Annotation_0306",cell_size = 1)+ggplot2::scale_colour_manual(values = col)

Mono.cds <- orderCells(Mono.cds, root_state = c('4'))

plot_cell_trajectory(Mono.cds, color_by = "Pseudotime",cell_size = 1)+viridis::scale_color_viridis(option="viridis")
ggsave(paste0(root,'_Pseudotime.pdf'),width = 5,height=5)#width = 9.62,height = 4.81)

plot_cell_trajectory(Mono.cds, color_by = "Annotation_0306",cell_size = 1)+ggplot2::scale_colour_manual(values = col)
ggsave(paste0(root,'_Annotation_0305.pdf'),width = 5,height = 5)#width = 9.62,height = 4.81)
plot_cell_trajectory(Mono.cds, color_by = "EGC",cell_size = 1)
ggsave(paste0(root,'EGC.pdf'),width = 5,height = 5)

####
### monocle 2
Time_diff <- differentialGeneTest(Mono.cds, cores = 1, 
                                  fullModelFormulaStr = "~sm.ns(Pseudotime)")
write.csv(Time_diff, "Time_diff_all.csv", row.names = F)

sig_gene_names <- subset(Time_diff, qval < 0.01)
p = plot_pseudotime_heatmap(Mono.cds[sig_gene_names$gene_short_name,],
                        num_clusters = 4,
                        cores = 1,
                        show_rownames = T,return_heatmap=T)


cluster = as.data.frame(cutree(p$tree_row,k=4))
colnames(cluster) = 'cluster'
cluster$id = row.names(cluster)
Time_diff_sig1 = merge(sig_gene_names,cluster,by.x = 'gene_short_name',by.y = 'id',sort =F)
or = c('3','1','4','2')
Time_diff_sig1$cluster = as.character(Time_diff_sig1$cluster)
Time_diff_sig1$cluster = factor(Time_diff_sig1$cluster,levels = or)
Time_diff_sig2 = Time_diff_sig1[order(Time_diff_sig1$cluster),]

p = plot_pseudotime_heatmap(Mono.cds[Time_diff_sig2$gene_short_name,],
                            num_clusters = 4,
                            cores = 1,cluster_rows=F,
                            show_rownames = F,return_heatmap=T)


p
ggsave(plot=p,filename = 'sig.heatmap.pdf')
write.csv(Time_diff_sig2,'Time_diff_sig2.csv',row.names=F)

###
lung_genes <- row.names(subset(fData(Mono.cds),
                               gene_short_name %in% c("AMEX60DD043470", "AMEX60DD046284","AMEX60DD029446")))

cds_subset <- Mono.cds[lung_genes,]
plot_genes_in_pseudotime(cds_subset,
                               color_by = "Annotation_0306",
                               ncol = 1)+ggplot2::scale_colour_manual(values = col)
ggsave('genes_in_pseudotime.pdf')

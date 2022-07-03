library(Seurat)
library(monocle3)
library(ggplot2)
library(viridis)
packageVersion('monocle3')

setwd('../data/FigureS15_data/FigureS15_Monocle3')
data = readRDS('Batch1_Stage57_telencephalon_rep2_DP8400015649BRD5_1.Annotation_0305.rds')

Idents(data) = 'Annotation_0305'
s15_c = c('VLMC','CP','DNBL5','Immature CMPN')
data = subset(data,idents= s15_c,invert = TRUE)

counts = data@assays$Spatial@counts
meta = data@meta.data
obj = CreateSeuratObject(counts= counts,meta.data =meta)
obj <- SCTransform(obj, verbose = FALSE)
obj <- RunPCA(obj, verbose = FALSE,npcs = 8)
obj <- RunUMAP(obj, dims = 1:6, verbose = FALSE,min.dist=0.4)

obj <- FindNeighbors(obj, dims = 1:30, verbose = FALSE)
obj <- FindClusters(obj, verbose = FALSE)
saveRDS(obj,paste0('Batch1_Stage57_telencephalon_rep2_DP8400015649BRD5_1.rds'))
co = read.csv('Inj_24_Adult_Develop_color_0305.txt',sep='\t')
celltype=as.data.frame(unique(data$Annotation_0305))
colnames(celltype) ='celltype'
co = merge(co,celltype,by.x='order',by.y='celltype',sort=F)
cor = co$Color
names(cor) = co$order

library(SeuratWrappers)
cds <- as.cell_data_set(obj)
cds <- cluster_cells(cds = cds, reduction_method = "UMAP")
cds <- learn_graph(cds, use_partition = TRUE)

get_earliest_principal_node <- function(cds, time_bin='DEGC'){
  # 首先找到出现在最早时间区间的细胞ID
  cell_ids <- which(colData(cds)[, "Annotation_0305"] == time_bin)
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
      (which.max(table(closest_vertex[cell_ids,]))))]

  root_pr_nodes

}

cds = order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds),reduction_method = "UMAP")

p1 = plot_cells(
  cds = cds,
  color_cells_by = "pseudotime",
  show_trajectory_graph = TRUE,
  group_cells_by = 'cluster',cell_size = 1,label_leaves=FALSE
 )+scale_color_viridis(option = "D")

ggsave(plot=p1,file=paste0('monocle3.pdf'),width=5,height=4)

p1 = plot_cells(cds,
           color_cells_by = "Annotation_0305",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,cell_size = 1, label_cell_groups=FALSE)+ ggplot2::scale_color_manual(values=cor)

ggsave(plot=p1,file=paste0('Celltype.pdf'),width=5.5,height=4)

saveRDS(cds,paste0('Batch1_Stage57_telencephalon_rep2_DP8400015649BRD5_1_monocle.rds'))
pseudotime <- as.data.frame(cds@principal_graph_aux$UMAP$pseudotime)
colnames(pseudotime) = 'pseudotime'
meta = as.data.frame(cds@colData)
meta$pseudotime = pseudotime$pseudotime
write.csv(meta,paste0('Batch1_Stage57_telencephalon_rep2_DP8400015649BRD5_1_monocle.anno.csv'),quote=F)


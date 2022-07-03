library(Seurat)
library(ggplot2)
library(SeuratWrappers)
library(monocle3)

setwd('../data/FigureS11_data/FigureS11B_Monocle3')
meta = read.csv('Batch2_Injury_2DPI_rep1_SS200000147BL_D5_Annotation_0306.meta.csv',row.names =1)
data$inj_M_L = meta$inj_M_L
Idents(data) = 'inj_uninj'
data = subset(data,idents = 'inj')
Idents(data) = 'Annotation_0306'

data$EGC = data$Annotation_0306
data$EGC[which(data$EGC=='REAEGC' & data$inj_M_L == 'inj_D_LP')] = 'REAEGC_LP'
data$EGC[which(data$EGC=='REAEGC')] = 'REAEGC_MP'

counts = data@assays$Spatial@counts
meta = data@meta.data
obj = CreateSeuratObject(counts= counts,meta.data =meta)
obj <- SCTransform(obj, verbose = FALSE)
obj <- RunPCA(obj, verbose = FALSE,npcs=3)
obj <- RunUMAP(obj, dims = 1:3, verbose = FALSE,min.dist=0.5)
obj <- FindNeighbors(obj, dims = 1:3, verbose = FALSE)
obj <- FindClusters(obj, verbose = FALSE)

saveRDS(obj,paste0('Batch2_Injury_2DPI_rep1_SS200000147BL_D5.rds'))

### Monocle3
cds <- as.cell_data_set(obj)
cds <- cluster_cells(cds = cds, reduction_method = "UMAP")
cds <- learn_graph(cds, use_partition = TRUE)

get_earliest_principal_node <- function(cds, time_bin='WNTEGC'){
  cell_ids <- which(colData(cds)[, "Annotation_0306"] == time_bin)
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

ggsave(plot=p1,file=paste0('monocle3.pdf'),width=4,height=3)

co = read.csv('Inj_24_Adult_Develop_color_0305.txt',sep='\t')
celltype=as.data.frame(unique(obj$Annotation_0306))
colnames(celltype) ='celltype'
co = merge(co,celltype,by.x='order',by.y='celltype',sort=F)
cor = co$Color
names(cor) = co$order


p1 = plot_cells(cds,
           color_cells_by = "Annotation_0306",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,cell_size = 1, label_cell_groups=FALSE)+ ggplot2::scale_color_manual(values=cor)

ggsave(plot=p1,file=paste0('Celltype.pdf'),width=4.5,height=3)

saveRDS(cds,paste0('Batch2_Injury_2DPI_rep1_SS200000147BL_D5_monocle.rds'))
pseudotime <- as.data.frame(cds@principal_graph_aux$UMAP$pseudotime)
colnames(pseudotime) = 'pseudotime'
meta = as.data.frame(cds@colData)
meta$pseudotime = pseudotime$pseudotime
write.csv(meta,paste0('Batch2_Injury_2DPI_rep1_SS200000147BL_D5_monocle.anno.csv'),quote=F)


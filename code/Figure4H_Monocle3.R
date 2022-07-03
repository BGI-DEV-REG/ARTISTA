library(Seurat)
library(Monocle3)
library(ggplot2)
library(viridis)

setwd('../data/Figure4_data/Figure4H')
slice=list.files(pattern="Batch")

rds=list()

for(i in seq(length(slice))){
 data=readRDS(slice[i])
 Idents(data) = 'Annotation_0306'
 ref = unique(data$Annotation_0306)
 s = c("REAEGC","IMN","RIPC1","RIPC2","RIPC3","RIPC4","NPTXEX")
 ss = intersect(ref,s)
 data = subset(data,idents= ss)
 Idents(data) = 'inj_uninj'
 data = subset(data,idents= "inj")
 Idents(data) = 'D_V'
 data = subset(data,idents= "D")
 counts = data@assays$Spatial@counts
 meta = data@meta.data
 data1 = CreateSeuratObject(counts= counts,meta.data =meta)
 name=gsub(".rds","",slice[i])
 data1$Batch=name
 rds[[i]]=data1
 print(nrow(data1@meta.data))
 print(paste0("Read ",name," Finished!!!"))
}

obj <- merge(rds[[1]], y=rds[2:length(slice)])
obj$batch <- unlist(lapply(strsplit(as.character(obj$Batch),"_"),"[",1))

ifnb.list <- SplitObject(object = obj, split.by = "batch")

ifnb.list <- lapply(X = ifnb.list, FUN = SCTransform, vars.to.regress = "nCount_RNA",assay = 'RNA')
features <- SelectIntegrationFeatures(object.list = ifnb.list, nfeatures = 3000)
ifnb.list <- PrepSCTIntegration(object.list = ifnb.list, anchor.features = features)
ifnb.list <- lapply(X = ifnb.list, FUN = RunPCA, features = features)

immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, normalization.method = "SCT",
    anchor.features = features, dims = 1:30, reduction = "rpca", k.anchor = 20)
immune.combined <- IntegrateData(anchorset = immune.anchors, normalization.method = "SCT", dims = 1:30)

immune.combined <- RunPCA(immune.combined, verbose = FALSE)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:30)

immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindClusters(immune.combined, resolution = 0.5)
saveRDS(immune.combined,paste0("result_in_D_reaegc_SCT_rpca.rds"))

### adjust
### load data
obj = readRDS('result_in_D_reaegc_SCT_rpca.rds')
Idents(obj) = obj$Annotation_0306

### removed RIPC2
obj = subset(obj,idents='RIPC2',invert = TRUE)
Idents(obj) = obj$Batch

### removed Control
obj = subset(obj,idents='Batch1_Injury_control_FP200000239BL_E3',invert = TRUE)
tag2=unlist(lapply(strsplit(as.character(obj$Batch),"_"),"[",3))
tag3=unlist(lapply(strsplit(as.character(obj$Batch),"_"),"[",4))
obj$abbBatch = paste(tag2,tag3,sep='_')

###
meta = as.data.frame(obj@meta.data)
Idents(obj) =obj$abbBatch

### removed regenaration early stage
sub = subset(obj,idents = c("2DPI_rep1","2DPI_rep2","2DPI_rep3","5DPI_rep1","5DPI_rep2","5DPI_rep3","10DPI_rep1","10DPI_rep2","10DPI_rep4"))
Idents(sub) =sub$Annotation_0306
sub = subset(sub,idents = c("NPTXEX"))
cells =rownames(sub@meta.data)
obj = subset(obj, cells = cells,invert = TRUE)

obj <- RunUMAP(obj, reduction = "pca", dims = 1:35,min.dist = 0.1)

### color
co = read.csv('Inj_24_Adult_Develop_color_0305.txt',sep = '\t')
ct = as.data.frame(unique(obj$Annotation_0306))
colnames(ct) = 'celltype'
co = merge(co,ct,by.x='order',by.y='celltype')
cor = co$Color
names(cor) = co$order

theme_black_black = theme(panel.border = element_blank(), axis.text = element_blank(), panel.grid = element_blank(), axis.ticks = element_blank(), plot.background = element_rect(fill = "black"), panel.background = element_rect(fill = "black"), axis.title = element_blank(), legend.text=element_text(color="black"), plot.title = element_text(color = "white"))

p <- DimPlot(obj, cols = cor,group.by='Annotation_0306',raster=FALSE)+theme_black_black
ggplot2::ggsave(plot = p ,file = paste0('result_in_D_reaegc.UMAP_adjust.pdf'),width = 5,height = 4)

### color
co = read.csv('RGC_color.txt',sep='\t')
cor = co$Color
names(cor) = co$order
p <- DimPlot(obj, cols = cor,group.by='abbBatch',raster=FALSE)+theme_black_black
ggplot2::ggsave(plot = p ,file = paste0('result_in_D_reaegc.UMAP_adjust_truetime.pdf'),width = 5,height = 4)

### Monocle 3
cds <- as.cell_data_set(obj)

cds <- cluster_cells(cds = cds, reduction_method = "UMAP")
cds <- learn_graph(cds, use_partition = TRUE)

get_earliest_principal_node <- function(cds, time_bin="REAEGC"){
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
saveRDS(cds,'REAEGC_0307.rds')

p1 = plot_cells(
  cds = cds,
  color_cells_by = "pseudotime",
  show_trajectory_graph = TRUE,
  group_cells_by = 'cluster'
 )+scale_color_viridis(option = "D")+theme_black_black

ggsave(plot=p1,file='monocle3.pdf',width=6,height=5)

pseudotime <- as.data.frame(cds@principal_graph_aux$UMAP$pseudotime)
colnames(pseudotime) = 'pseudotime'
meta = as.data.frame(cds@colData)
meta$pseudotime = pseudotime$pseudotime
write.csv(meta,'monocle.anno.csv',quote=F)

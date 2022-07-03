library(Seurat)
setwd("../data/Figure3_data/Figure3F")
stagedata=readRDS("../data/Figure3_data/Figure3F/stage_merge_merge.rds")
stagedata$Annotation_0306 <- stagedata$Annotation_0305
indata=readRDS("../data/Figure3_data/Figure3F/in_2_20_merge_2_20.rds")
obj <- merge(stagedata, indata)
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
immune.combined <- FindClusters(immune.combined, resolution = 1)
saveRDS(immune.combined,"Injury_dev_inter0307_V2_SCT_rpca.rds")
rds=RunUMAP(immune.combined,dims = 1:20,min.dist = 0.1)
p <- DimPlot(rds,group.by = "Annotation_0306") + scale_color_manual(values = c("#faa307","#8983BF","#A30059","#00868B"))

sum2=read.table("../data/Figure3_data/Figure3F/Summary_Gene_Annotaion_0707.xls",header=T,fill=TRUE, na.strings = "")
#sum2$idmouse <- paste0(sum2$Axolotl_ID,":",sum2$mm_gene)
obj.integrated.Find <- rds
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
sum1=read.table("../data/Figure3_data/Figure3F/dotplot_gene_20220307.txt",header=T,sep = "\t")
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


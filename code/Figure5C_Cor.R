library(ggplot2)
library(Seurat)
library(ArchR)

rds1=readRDS("15DPI_4.rds")
rds2=readRDS("Control(Juv).rds")
rds3=readRDS("Stage44.rds")
rds4=readRDS("Stage54.rds")
rds5=readRDS("Stage57.rds")

rds1$Batch="15DPI_4"
rds2$Batch="Juv"
rds3$Batch="Stage44"
rds4$Batch="Stage54"
rds5$Batch="Stage57"

meta1=rds1@meta.data
meta2=rds2@meta.data
meta3=rds3@meta.data
meta4=rds4@meta.data
meta5=rds5@meta.data

rds1$Annotation_Final=rds1$Annotation_0306
rds2$Annotation_Final=rds2$Annotation_0306
rds3$Annotation_Final=rds3$Annotation_0305
rds4$Annotation_Final=rds4$Annotation_0305
rds5$Annotation_Final=rds5$Annotation_0305

meta1_sub=subset(meta1,meta1$inj_uninj=="inj" & meta1$D_V=="D")
meta2_sub=subset(meta2,meta2$inj_uninj=="inj" & meta2$D_V=="D")
meta3_sub=subset(meta3,meta3$inj_uninj=="inj" & meta3$D_V=="D")
meta4_sub=subset(meta4,meta4$inj_uninj=="inj" & meta4$D_V=="D")
meta5_sub=subset(meta5,meta5$inj_uninj=="inj" & meta5$D_V=="D")

rds1_sub=subset(rds1,cells = rownames(meta1_sub))
rds2_sub=subset(rds2,cells = rownames(meta2_sub))
rds3_sub=subset(rds3,cells = rownames(meta3_sub))
rds4_sub=subset(rds4,cells = rownames(meta4_sub))
rds5_sub=subset(rds5,cells = rownames(meta5_sub))

rds_all=merge(rds1_sub,c(rds2_sub,rds3_sub,rds4_sub,rds5_sub))
celltype=c("dEGC","reaEGC","ribEGC","wntEGC","sfrpEGC")
Idents(rds_all)=rds_all$Annotation_Final
rds_all_celltype=subset(rds_all,idents = celltype)
rds_all_celltype1=CreateSeuratObject(counts = rds_all_celltype@assays$Spatial@counts,project = "artista",min.cells = 2)
rds_all_celltype1=AddMetaData(rds_all_celltype1,metadata = rds_all_celltype@meta.data)

###SCT_RPCA
ifnb.list <- SplitObject(rds_all_celltype, split.by = "Batch")
ifnb.list <- lapply(X = ifnb.list, FUN = SCTransform, method = "glmGamPoi",assay = "Spatial")
features <- SelectIntegrationFeatures(object.list = ifnb.list, nfeatures = 2000)
ifnb.list <- PrepSCTIntegration(object.list = ifnb.list, anchor.features = features)
ifnb.list <- lapply(X = ifnb.list, FUN = RunPCA, features = features)
immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, normalization.method = "SCT",
                                         anchor.features = features, dims = 1:30, reduction = "rpca", k.anchor = 20)
immune.combined.sct <- IntegrateData(anchorset = immune.anchors, normalization.method = "SCT", dims = 1:30,k.weight = 50)
immune.combined.sct <- RunPCA(immune.combined.sct, verbose = FALSE)
immune.combined.sct <- RunUMAP(immune.combined.sct, reduction = "pca", dims = 1:50,min.dist = 0.1)

###Cor
immune.combined.sct$Group=paste0(immune.combined.sct$Annotation_Final,"__",immune.combined.sct$Batch)
Idents(immune.combined.sct)=immune.combined.sct$Group
ave=AverageExpression(immune.combined.sct)
ave_RNA=as.data.frame(ave$integrated)

pheatmap::pheatmap(cor(ave_RNA[,-4]),cellwidth = 10,cellheight = 10,treeheight_row = 15,treeheight_col = 15,
                   color = c(paletteContinuous(set = "solarExtra", n = 30),"#A31D1D","#A31D1D","#A31D1D","#A31D1D","#A31D1D"))

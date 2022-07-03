#####Figure S7##########
###### mut  section inter####
co <- c('#2092F5','#1F96F6','#1E99F7','#1C9CF8','#1B9FF9','#1AA3FA','#16ADFD','#14B0FE','#16B3FE','#1FB5FD','#29B7FC','#32BAFA','#57C2F5','#61C5F4','#6AC7F3','#74C9F1','#7DCBF0','#95CFEA','#A7D1E4','#BAD4DE','#BED4DC','#C2D4D9','#C5D4D3','#CFD4C2','#D3D4BD','#D6D3B7','#D9D3B2','#DDD3AC','#E0D3A7','#E3D3A1','#EBD08B','#EDCD81','#F3C359','#F8BB3A','#F9B830','#FBB626','#FCB31C','#FBAA1A','#F99F1C','#F7941D','#F5891E','#F37E20','#F17321','#EF6822','#ED5D24','#EB5225','#E94726','#E73B27','#E53029','#E53028','#E53027','#E53024','#E53023','#E53022','#E53021','#E53019','#E53018','#E53017','#E53016','#E53015','#E53013','#E53012','#E53009','#E53007','#E22905','#DC2828','#D72727','#D22626','#CD2525','#C72424','#C22323','#BD2222','#B82121','#B22020','#AD1F1F','#A81E1E','#A31D1D')
colorlist<- c( "#604E97","#F6A600", "#B3446C","#0cf56d","#E25822","#DCD300","#882D17", "#8DB600","#654522")
###################2DPI####

obj = readRDS('../data/FigureS8_data/Injury_2DPI/Injury_2DPInfeatures3000_SCT_rpca.rds')
obj$batch <- unlist(lapply(strsplit(as.character(obj$Batch),"_"),"[",4))
Idents(obj) <- "Annotation_0306"
color <- read.table("D:/02.project/06RYfish/02brain_dev/00.modul/Inj_24_Adult_Develop_color_0307.txt",sep = "\t",comment.char = "*",header = T)
celltype_ry <- as.data.frame(unique(obj$Annotation_0306))
colnames(celltype_ry) <- "name"
cellcolor <- merge(celltype_ry,color,by.x = "name",by.y = "order")
obj$Annotation_0306_2 <- factor(obj$Annotation_0306,levels = cellcolor$name)
pdf("../data/FigureS8_data/Injury_2DPI/2DPInfeatures3000_Celltype_umap_0507.pdf",width = 9.33,height = 7.22)
p2 <- DimPlot(obj,group.by = "Annotation_0306_2",cols = cellcolor$Color,raster=FALSE) + ggtitle("2DPInfeatures3000")
print(p2)
dev.off()
pdf("../data/FigureS8_data/Injury_2DPI/2DPInfeatures3000_batch_umap_0507.pdf",width = 7.45,height = 6.77)
p2 <- DimPlot(obj, cols = colorlist,group.by='batch',raster=FALSE) + ggtitle("2DPInfeatures3000")
print(p2)
dev.off()
obj$group <- paste(obj$batch,obj$Annotation_0306,sep = '_') 
Idents(obj) <- obj$group
cluster.averages <- AverageExpression(obj)
inte <- as.data.frame(cluster.averages$integrated)
cluster_cor <- cor(inte,inte)
p1 <- pheatmap(cluster_cor)
gn=rownames(cluster_cor)[p1$tree_row[["order"]]]
sn=colnames(cluster_cor)[p1$tree_col[["order"]]]
new_test=cluster_cor[gn,sn]
p2 <- pheatmap(new_test,color =co,cluster_rows = F,cluster_cols = F,border_color = NA)
pdf("../data/FigureS8_data/Injury_2DPI/2DPInfeatures3000_Celltype_cor_0507.pdf",width = 9.4,height = 9) 
print(p2)
dev.off()

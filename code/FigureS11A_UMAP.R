library(Seurat)
library(ggplot2)

setwd('../data/FigureS11_data/FigureS11A')
data = readRDS('Batch2_Injury_2DPI_rep1_SS200000147BL_D5.rds')

data <- RunUMAP(data, dims = 1:20,verbose = FALSE,min.dist=0.5,seed.use = 31)
p = DimPlot(data, label = FALSE,group.by='Annotation_0306',cols=cor)
ggsave(plot=p,file=paste0('Batch2_Injury_2DPI_rep1_SS200000147BL_D5.rds_Annotation_0306.pdf'),width=6,height = 5)

co = read.csv('Inj_24_Adult_Develop_color_0305.txt',sep='\t')
celltype=as.data.frame(unique(data$Annotation_0306))
colnames(celltype) ='celltype'
co = merge(co,celltype,by.x='order',by.y='celltype',sort=F)
cor = co$Color
names(cor) = co$order
p = DimPlot(data, label = TRUE,group.by='Annotation_0306',cols=cor)
ggsave(plot=p,file=paste0('Batch2_Injury_2DPI_rep1_SS200000147BL_D5.rds_Annotation_0306.pdf'),width=6,height = 5)


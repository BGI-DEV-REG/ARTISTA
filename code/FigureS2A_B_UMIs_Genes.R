library(Seurat)
library(ggplot2)
library(reshape2)
library(dplyr)

#### UMI ### Gene ###UMI

setwd('../data/FigureS2_data/FigureS2A_B')
data <-  read.table('Batch1_Adult_telencephalon_rep2_DP8400015234BLA3_1.SCT_Removed_0305.csv',sep = ',',header = T,row.names = 1)
#unique(da$Batch)
#data <- da[which(da$Batch=='Adult_telencephalon_rep2_DP8400015234BLA3_1'),]
colnames(data)
cor = read.csv('Inj_24_Adult_Develop_color_0305.txt',header=T,sep='\t')
co = cor$Color
names(co) = cor$order

####修改
p = ggplot(data,aes(x=Annotation_0305,y=nCount_Spatial,fill=Annotation_0305))+
  geom_boxplot(outlier.size=0.01)+
  scale_fill_manual(values=co)+
  theme_bw()+theme(panel.border = element_blank(),panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   axis.line = element_line(colour = "black"),strip.background = element_rect(
                     color = "white", fill = "white"),axis.text.x =  element_text(hjust = 0.5, angle = 90))+ylab('UMI count')+ylim(c(0,8000))

ggsave(plot= p, file ='Adult_UMIs.pdf')
## Gene 
p = ggplot(data,aes(x=Annotation_0305,y=nFeature_Spatial,fill=Annotation_0305))+
  geom_boxplot(outlier.size=0.01)+
  scale_fill_manual(values=co)+
  theme_bw()+theme(panel.border = element_blank(),panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   axis.line = element_line(colour = "black"),strip.background = element_rect(
                     color = "white", fill = "white"),axis.text.x =  element_text(hjust = 0.5, angle = 90))+ylab('Gene count')+ylim(c(0,4000))

ggsave(plot= p, file ='Adult_Genes.pdf')

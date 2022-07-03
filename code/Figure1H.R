##
library(Seurat)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(DOSE)
library(org.Hs.eg.db)
library(topGO)
library(clusterProfiler)
library(pathview)
library(org.Mm.eg.db)

deg=read.table("/../Figure1_RGC_MK.xls",header=T)
deg2=subset(deg,deg$avg_log2FC>0.7 & deg$p_val<0.01)
deg2=deg2[order(deg2$avg_log2FC,decreasing = TRUE),]
deg22=deg2

deg7=list()
for(i in unique(deg22$cluster)){
  deg7[[i]] <- deg22 %>%
    filter(cluster==i) %>%
    pull(mm_gene) %>%
    unique()
}

compGO <- compareCluster(geneCluster   = deg7,
                         fun           = "enrichGO",
                         pvalueCutoff  = 0.1,
                         pAdjustMethod = "BH",
                         keyType = 'SYMBOL',
                         OrgDb = org.Mm.eg.db,
                         ont = "BP",
                         readable = FALSE)

p=dotplot(compGO,showCategory = 15)+
  theme(axis.text.x = element_text(size = 10, vjust = 0.5, hjust = 0.5, angle = 90,color = "black"))

out=compGO@compareClusterResult
select=as.data.frame(out)

##
select=select[-27,]
RIBRGC=select[select$Cluster=="RIBRGC",]

###
RIBRGC_gene=RIBRGC_gene[!is.na(RIBRGC)]
RIBRGC_gene=as.data.frame(RIBRGC_gene)
colnames(RIBRGC_gene)="mm_gene"
#RIBRGC_gene$Mouse_EntreID=as.integer(RIBRGC_gene$Mouse_EntreID)
RIBRGC_gene2=left_join(RIBRGC_gene,deg22,by="mm_gene")
#


SFRPRGC_gene2=SFRPRGC_gene2[SFRPRGC_gene2$cluster=="SFRPRGC",]
RIBRGC_gene2=RIBRGC_gene2[RIBRGC_gene2$cluster=="RIBRGC",]
SFRPRGC_gene2=SFRPRGC_gene2[SFRPRGC_gene2$cluster=="SFRPRGC",]

all_gene=rbind.data.frame(WBTRGC_gene2,SFRPRGC_gene2,RIBRGC_gene2)

##

rds=readRDS("/../Figure1_Adult_RGC.rds")
p=DotPlot(object = rds, features =all_gene$gene,group.by = "Anno")+
  theme(axis.text.x=element_text(angle=90, hjust=1))+
  scale_color_viridis_c(option = "B",direction = 1,end = 0.8)+
  scale_x_discrete("New Type", labels = all_gene$mm_gene)
###GO





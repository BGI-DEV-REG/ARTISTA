library(Seurat)
library(ggplot2)

rds=readRDS("Batch1_Adult_telencephalon_rep2_DP8400015234BLA3_1.SCT.rds")
rds1=subset(rds,idents = c(20,22,15))
rds1$Anno="NA"
rds1$Anno[rds1$seurat_clusters==20]="ribEGC"
rds1$Anno[rds1$seurat_clusters==22]="wntEGC"
rds1$Anno[rds1$seurat_clusters==15]="sfrpEGC"
rds1$Anno=factor(rds1$Anno,levels = c("wntEGC","sfrpEGC","ribEGC"))
Idents(rds1)=rds1$Anno
col=c("#71c33a","#e14aec","#A52A2A")
VlnPlot(rds1,features = "AMEX60DD001640",pt.size = 0)+ggtitle("Sox2")+scale_fill_manual(values = col)
VlnPlot(rds1,features = "AMEX60DD022108",pt.size = 0)+ggtitle("Vim")+scale_fill_manual(values = col)
VlnPlot(rds1,features = "AMEX60DD041912",pt.size = 0)+ggtitle("Slc1a3")+scale_fill_manual(values = col)
VlnPlot(rds1,features = "AMEX60DD009861",pt.size = 0)+ggtitle("Gfap")+scale_fill_manual(values = col)

###
mk_rgc=FindAllMarkers(rds1)
anno=read.table("Summary_Gene_Annotaion_0707.xls",header = T,sep = "\t")
colnames(anno)[1]="gene"
merge_gene=merge(mk_rgc,anno,by="gene")
write.table(merge_gene,"Figure1_RGC_DEG/Figure1_RGC_MK.xls",sep = "\t",quote = FALSE,row.names = FALSE)

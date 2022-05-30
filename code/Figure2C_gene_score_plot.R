rm(list=ls())
gc()
library(data.table)
library(dplyr)
library(Seurat)
library(RColorBrewer)
library(ggthemes)
library(ggplot2)
library(plyr)
library(ggpubr)
library(viridis)
library(reshape2)
library(pheatmap)
library(ggthemes)
mytheme <- theme(legend.position = 'right',
                 plot.title = element_text(size = rel(1.2), hjust = 0),
                 axis.text = element_text(angle = 45,size = rel(1.0)),
                 axis.title = element_text(size = rel(1.2)),
                 legend.title = element_text(size = rel(1.2),hjust = 0),
                 panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))

###  RY data_compare
setwd("D:/02.project/06RYfish/02brain_dev/00.modul/")
obj.integrated.Find<- readRDS(paste("D:/02.project/06RYfish/02brain_dev/00.python_spifig/data/all_dev_recluster.rds",sep=''))
ref =read.table("D:/02.project/06RYfish/02brain_dev/00.modul/Summary_Gene_Annotaion_0707.xls",header=T,fill=TRUE, na.strings = "")
geneTRan <- read.table("D:/02.project/06RYfish/02brain_dev/00.modul/data_Cell_Cycle.uniq.txt")
a <- geneTRan$V1
DefaultAssay(obj.integrated.Find)<-"RNA"
obj <- AddModuleScore(
  object = obj.integrated.Find,
  features = list(unique(a)),
  ctrl = 100,
  nbin = 24,
  name = "Cell_Cycle"
)
lymeta <- as.data.frame(obj@meta.data)  
lymeta$batch3 <- factor(lymeta$batch,levels = c('Stage44','Stage54','Stage57','Injury','Adult','Meta'))
library(ggplot2)
p <- ggplot(data=lymeta, aes(x=batch3, y=Cell_Cycle1, fill=batch3)) + geom_violin(alpha=0.8,width=1) +geom_boxplot(position = position_dodge(width = 0.1), outlier.size=0, width = 0.1, lwd=0.01 ,show.legend = FALSE) +
  scale_fill_manual(values=c(viridis(n=8,option = 'B',direction = -1)[0:9]))
p <- p + mytheme                  
pdf("F2C_Cell_Cycle_0818RNA.pdf",width = 5.7,height = 3.1)
print(p)
dev.off()

geneTRan <- read.table("D:/02.project/06RYfish/02brain_dev/00.modul/data_Translation.uniq.txt")
a <- geneTRan$V1
DefaultAssay(obj.integrated.Find)<-"RNA"
obj <- AddModuleScore(
  object = obj.integrated.Find,
  features = list(unique(a)),
  ctrl = 100,
  nbin = 24,
  name = "Translation"
)
lymeta <- as.data.frame(obj@meta.data)  
lymeta$batch3 <- factor(lymeta$batch,levels = c('Stage44','Stage54','Stage57','Injury','Adult','Meta'))
library(ggplot2)
p <- ggplot(data=lymeta, aes(x=batch3, y=Translation1, fill=batch3)) + geom_violin(alpha=0.8,width=1) +geom_boxplot(position = position_dodge(width = 0.1), outlier.size=0, width = 0.1, lwd=0.01 ,show.legend = FALSE) +
  scale_fill_manual(values=c(viridis(n=8,option = 'B',direction = -1)[0:9]))
p <- p + mytheme                  
pdf("F2C_Translation_0818RNA.pdf",width = 5.7,height = 3.1)
print(p)
dev.off()

geneTRan <- read.table("D:/02.project/06RYfish/02brain_dev/00.modul/data_NSC_gene.uniq.txt")
a <- geneTRan$V2
DefaultAssay(obj.integrated.Find)<-"RNA"
obj <- AddModuleScore(
  object = obj.integrated.Find,
  features = list(unique(a)),
  ctrl = 100,
  nbin = 24,
  name = "NSCgene"
)
lymeta <- as.data.frame(obj@meta.data)  
lymeta$batch3 <- factor(lymeta$batch,levels = c('Stage44','Stage54','Stage57','Injury','Adult','Meta'))
library(ggplot2)
p <- ggplot(data=lymeta, aes(x=batch3, y=NSCgene1, fill=batch3)) + geom_violin(alpha=0.8,width=1) +geom_boxplot(position = position_dodge(width = 0.1), outlier.size=0, width = 0.1, lwd=0.01 ,show.legend = FALSE) +
  scale_fill_manual(values=c(viridis(n=8,option = 'B',direction = -1)[0:9]))
p <- p + mytheme                  
pdf("F2C_NSCgene_0818RNA.pdf",width = 5.7,height = 3.1)
print(p)
dev.off()

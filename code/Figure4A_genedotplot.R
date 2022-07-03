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

obj.integrated.Find <- readRDS("../data/Figure4_data/Figure4A/RY_RNA_inter_0827.rds")
Idents(obj.integrated.Find) <- "Batch"
sub_set_obj.integrated.Find  <- subset(obj.integrated.Find,idents = c( "Injury_15DPI_C_rep1_FP200000266TR_E6.ssDNA.0722.rds","Injury_15DPI_rep2_FP200000266TR_E2.ssDNA.0722.rds","Injury_15DPI_rep3_FP200000266TR_E3.ssDNA.0722.rds","Injury_15DPI_rep4_FP200000266TR_E4.ssDNA.0722.rds"))
Idents(sub_set_obj.integrated.Find) <- "Celltype"
sub_sub_set_obj.integrated.Find  <- subset(sub_set_obj.integrated.Find,idents = c( "RIPC1","REARGC","Immature NPTXEX"))
sub_sub_set_obj.integrated.Find$CelltypeV1 <- factor(sub_sub_set_obj.integrated.Find$Celltype, levels = c("REARGC","RIPC1","Immature NPTXEX"))
DotPlot(sub_sub_set_obj.integrated.Find, features = c("AMEX60DD052557",
"AMEX60DD008688",
"AMEX60DD004311",
"AMEX60DD002984",
"AMEX60DD003720",
"AMEX60DD039970",
"AMEX60DD032735",
"AMEX60DD047001",
"AMEX60DD046284",
"AMEX60DD016851",
"AMEX60DD022108") , cols = c("white", "red"),group.by = "CelltypeV1", dot.scale = 8) + RotatedAxis() + coord_flip()

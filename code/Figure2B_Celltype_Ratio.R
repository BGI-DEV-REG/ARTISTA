library(Seurat)
library(ggplot2)
library(reshape2)
library(dplyr)
library(pheatmap)

### Celltype
setwd('../data/Figure2_data/Figure2B')

### load meta data
Adult = read.csv('Adult_telencephalon_rep2_DP8400015234BLA3_1_SCT_Removed.meta_0305.csv')[c('Batch','Annotation_0305')]
Meta = read.csv('Batch1_Meta_telencephalon_rep1_DP8400015234BLB2_1_Annotation_0305.meta.csv')[c('Batch','Annotation_0305')]
Control  = read.csv('Batch1_Injury_control_FP200000239BL_E3_Annotation_0306.csv')[c('Batch','Annotation_0306')]
Stage44 = read.csv('Batch1_Stage44_telencephalon_rep2_FP200000239BL_E4_Annotation_0305.meta.csv')[c('Batch','Annotation_0305')]
Stage54 = read.csv('Batch1_Stage54_telencephalon_rep2_DP8400015649BRD6_2_Annotation_0305.meta.csv')[c('Batch','Annotation_0305')]
Stage57 = read.csv('Batch1_Stage57_telencephalon_rep2_DP8400015649BRD5_1_Annotation_0305.meta.csv')[c('Batch','Annotation_0305')]
colnames(Control) = c('Batch','Annotation_0305')

Adult$batch[which(Adult$Batch=='Adult_telencephalon_rep2_DP8400015234BLA3_1')] <- 'Adult'
Meta$batch[which(Meta$Batch=='Meta_telencephalon_rep1_DP8400015234BLB2_1')] <- 'Meta'
Stage44$batch[which(Stage44$Batch=='Stage44_telencephalon_rep2_FP200000239BL_E4')] <- 'Stage44'
Stage54$batch[which(Stage54$Batch=='Stage54_telencephalon_rep2_DP8400015649BRD6_2')] <- 'Stage54'
Stage57$batch[which(Stage57$Batch=='Stage57_telencephalon_rep2_DP8400015649BRD5_1')] <- 'Stage57'
Control$batch[which(Control$Batch=='Injury_control_FP200000239BL_E3')] <- 'Control'

### rbind
da = rbind(Adult,Meta)
da = rbind(da,Control)
da =rbind(da,Stage44)
da =rbind(da,Stage54)
da =rbind(da,Stage57)
da =rbind(da,Control)

data <- da[which(da$batch %in% c('Adult','Meta','Control','Stage44','Stage54','Stage57')),]

### Count data
data_in <- data
data_in_count2 <- as.data.frame(data_in %>% group_by(batch,Annotation_0305) %>% summarise(count=n()))
meta_in_count4 <- as.data.frame(dcast(data_in_count2, batch~Annotation_0305 ,value.var = 'count'))
meta_in_count4[is.na(meta_in_count4)] <- 0
meta_section_in <- as.data.frame(data_in %>% group_by(batch) %>% summarise(count=n()))
meta_in_count4$Spot <- meta_section_in$count
colnames(meta_in_count4)

meta_in <- melt(meta_in_count4,id.vars=c('batch','Spot'))
meta_in$Ratio <- meta_in$value/meta_in$Spot

### Color
co1 = read.csv('Inj_Adult_Develop_color_0305.txt',sep='\t',header=T)
d = as.data.frame(unique(unique(meta_in$variable)))
colnames(d) = 'order'
cor = merge(co1,d,by ='order')
co = cor$Color
names(co) = cor$order

### Order stage
ord <- rev(c('Stage44','Stage54','Stage57','Control','Adult','Meta'))
meta_in$batch <- factor(meta_in$batch,levels = ord)
da <- meta_in[order(meta_in$batch),]

### Removed UnKnown, CP, VLMC cell types 
da<-da[which(da$variable !=c('UnKnown')),]
da<-da[which(da$variable !=c('CP')),]
da<-da[which(da$variable !=c('VLMC')),]

### Order cell types
l <- c(
  "Immature MSN",
  "Immature CMPN",
  "Immature NPTXEX",
  "Immature MPEX",
  "Immature CCKIN",
  "Immature DPEX",
  "Immature DMIN",
  "DNBL1",
  "DNBL2",
  "DNBL3",
  "DNBL4",
  "DNBL5",
  "OBNBL",
  "TLNBL1" ,"TLNBL2","TLNBL3","TLNBL4",
  "TLNBL" ,
  "DEGC",
  "WNTEGC",
  "SFRPEGC",
  "RIBEGC",
  "NPTXEX",
  "MSN",
  "MPEX",
  "CMPN",
  "DPEX",
  "SCGIN",
  "CCKIN","SCGNIN","NTNG1EX","SSTIN","MPIN", "NPYIN","OLIGO","TAC2IN","MCG"
  )

da$variable = factor(da$variable,levels =l)

### Plot
p = ggplot(da,aes(variable,batch,colour=variable,fill=variable))+geom_point(aes(size=Ratio),alpha=0.8)+theme_set(theme_bw())+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+scale_colour_manual(values = co)+scale_size_continuous(range=c(0,12))+
  theme(axis.text.x = element_text(vjust = 0.5, hjust = 0.5,angle = 90))

### Save
ggsave(plot = p, file = 'Dev_Ratio_Polt.pdf')

rm(list = ls())
library(ggplot2)
library(dplyr)
library(reshape2)
library(tidyr)
library(tibble)
library(ggradar)
setwd('D:/02.project/06RYfish/02brain_dev/00.python_spifig/03.RY_replot_2022/rataplot/00.data')
id = c('Batch1_Injury_15DPI_rep2_FP200000266TR_E2','Batch1_Injury_15DPI_rep3_FP200000266TR_E3',
       'Batch1_Injury_15DPI_rep4_FP200000266TR_E4','Batch1_Injury_30DPI_rep2_FP200000264BL_A6','Batch1_Injury_60DPI_rep3_FP200000264BL_A6',
       'Batch1_Injury_control_FP200000239BL_E3','Batch2_Injury_10DPI_rep1_SS200000147BL_B5','Batch2_Injury_10DPI_rep2_SS200000147BL_B2',
       'Batch2_Injury_10DPI_rep4_SS200000147BL_B3','Batch2_Injury_20DPI_rep1_SS200000147BL_B4','Batch2_Injury_20DPI_rep2_SS200000147BL_B4',
       'Batch2_Injury_20DPI_rep3_SS200000147BL_B5','Batch2_Injury_2DPI_rep1_SS200000147BL_D5','Batch2_Injury_2DPI_rep2_SS200000147BL_D5',
       'Batch2_Injury_2DPI_rep3_SS200000147BL_D4',
       'Batch2_Injury_5DPI_rep1_SS200000147BL_D2','Batch2_Injury_5DPI_rep2_SS200000147BL_D2','Batch2_Injury_5DPI_rep3_SS200000147BL_D3')
l = list()
for (i in id ) {
  data = read.csv(paste0(i,'_meta.csv'))
  data = data[,c('CellID','Batch','inj_uninj','D_V','Annotation_0306')]
  l[[i]] = data
}
data = do.call(rbind,l)
BIA <- as.data.frame(data %>% group_by(Batch, inj_uninj,Annotation_0306) %>% dplyr::summarise(BIA=n()))
BI = as.data.frame(data %>% group_by(Batch, inj_uninj) %>% dplyr::summarize(BI=n()))
dBIA = dcast(BIA, Batch ~ inj_uninj+Annotation_0306, value.var = 'BIA')
dBIA[is.na(dBIA)] <- 0
dBIA[dBIA <= 5] <- 5
dBI = dcast(BI, Batch ~ inj_uninj,value.var = 'BI')
inj = dBIA[,c(1:29)]
colnames(inj)
rownames(inj) = inj[,1]
inj = inj[,-1]
all_inj = as.data.frame(dBI$inj)
inj_r = inj/all_inj$`dBI$inj`
uninj = dBIA[,c(30:56)]
rownames(uninj) = dBIA$Batch
colnames(uninj)
all_uninj = as.data.frame(dBI$uninj)
uninj_r = uninj/all_uninj$`dBI$uninj`
colnames(inj_r)
colnames(uninj_r)
inj_r$b1 <- unlist(lapply(strsplit(as.character(rownames(inj_r)),"_"),"[",2))
inj_r$b2 <- unlist(lapply(strsplit(as.character(rownames(uninj_r)),"_"),"[",3))
inj_r$Batch = paste(inj_r$b1,inj_r$b2,sep = '_')
unique(inj_r$Batch)
uninj_r$b1 <- unlist(lapply(strsplit(as.character(rownames(uninj_r)),"_"),"[",2))
uninj_r$b2 <- unlist(lapply(strsplit(as.character(rownames(uninj_r)),"_"),"[",3))
uninj_r$Batch = paste(uninj_r$b1,uninj_r$b2,sep = '_')
l = sort(unique(data$Annotation_0306))
color_cell <- read.table("D:/02.project/06RYfish/02brain_dev/00.python_spifig/03.RY_replot_2022/rataplot/00.data/Inj_24_Adult_Develop_color_0307.txt",comment.char = "*",header = T,sep = "\t")
rownames(color_cell) <- color_cell$order
final_cell <- color_cell[sort(c("CMPN","DPEX","MPEX","MSN","SCGNIN","SFRPEGC","SSTIN","VLMC","IMN", "MCG","NPTXEX","REAEGC","WNTEGC","WSN")), ]
X <- as.character(final_cell$order)
Y <- as.character(final_cell$Color)
final_cell$order1 <- factor(final_cell$order,levels = X)
final_cell$Color1 <- factor(final_cell$Color,levels = Y)
df_final <- data.frame()
for (i in 1:length(final_cell$order1)) {
  dinj = as.data.frame(inj_r[,paste0('inj_',X[i])])
  colnames(dinj) = X[i]
  dinj$type = 'Inj'
  dinj$Batch = inj_r$Batch
  duninj =as.data.frame(uninj_r[,paste0('uninj_',X[i])])
  colnames(duninj) = X[i]
  duninj$type = 'uninj'
  duninj$Batch = uninj_r$Batch
  cd = rbind(dinj,duninj)
  cd_wide_1 <- spread(cd,type,X[i])
  cd_wide_1$fd <- log(cd_wide_1$Inj +1)/log(cd_wide_1$uninj + 1)
  rownames(cd_wide_1) <- cd_wide_1$Batch
  cd_wide_1 <- cd_wide_1[,-1]
  cd_wide_1_T <- t(cd_wide_1)
  cd_wide_1_T <- as.data.frame(cd_wide_1_T)
  data <- cd_wide_1_T[c("control_FP200000239BL","2DPI_rep1",
                        "2DPI_rep2","2DPI_rep3",
                        "5DPI_rep1","5DPI_rep2",
                        "5DPI_rep3","10DPI_rep1",
                        "10DPI_rep2","10DPI_rep4",
                        "15DPI_rep2","15DPI_rep3",
                        "15DPI_rep4",
                        "20DPI_rep1","20DPI_rep2",
                        "20DPI_rep3","30DPI_rep2",
                        "60DPI_rep3")]
  colnames(data) <- c("control_FP200000239BL","2DPI_rep1",
                      "2DPI_rep2","2DPI_rep3",
                      "5DPI_rep1","5DPI_rep2",
                      "5DPI_rep3","10DPI_rep1",
                      "10DPI_rep2","10DPI_rep4",
                      "15DPI_rep2","15DPI_rep3",
                      "15DPI_rep4",
                      "20DPI_rep1","20DPI_rep2",
                      "20DPI_rep3","30DPI_rep2",
                      "60DPI_rep3")
  final <- data[c("fd"), ]
  rownames(final) <- X[i]
  df_merge <- final
  df_merge[df_merge==Inf] <- 4
  df_merge[df_merge >= 4] <- 4
  df_merge2 <- unstack(within(stack(df_merge), values[is.nan(values)] <- 0))
  df_merge3 <- as.data.frame(t(df_merge2))
  #df_merge2[is.na(df_merge)] <- 0
  colnames(df_merge3) <- c("control_FP200000239BL","2DPI_rep1",
                           "2DPI_rep2","2DPI_rep3",
                           "5DPI_rep1","5DPI_rep2",
                           "5DPI_rep3","10DPI_rep1",
                           "10DPI_rep2","10DPI_rep4",
                           "15DPI_rep2","15DPI_rep3",
                           "15DPI_rep4",
                           "20DPI_rep1","20DPI_rep2",
                           "20DPI_rep3","30DPI_rep2",
                           "60DPI_rep3")
  rownames(df_merge3) <- X[i]
  df_merge3 <- df_merge3 %>% rownames_to_column("group")
  p1 <- ggradar(
    df_merge3, 
    values.radar = c("0", "2", "4"),
    grid.min = 0, grid.mid = 2, grid.max = 4,
    # Polygons
    group.line.width = 0.8, 
    group.point.size = 1,
    group.colours = c("#BA0900"),
    # Background and grid lines
    background.circle.colour = "white",
    gridline.mid.colour = "grey",
    legend.position = "bottom"
  )
  
  pdf(paste0(X[i],"_celltype_final_4.pdf"),height = 10,width = 10)
  print(p1)
  dev.off()
  df_final <- rbind(df_final, df_merge3)
  write.table(df_merge3,paste0(X[i],"_celltype_inj_final_10.txt"), sep="\t", row.names=T)
} 

#2DPI
df_final$I2DPI_mea <- apply(df_final[,3:5],1,mean)
#5DPI
df_final$I5DPI_mea <- apply(df_final[,6:8],1,mean)
#10DPI
df_final$I10DPI_mea <- apply(df_final[,9:11],1,mean)
#15DPI
df_final$I15DPI_mea <- apply(df_final[,12:14],1,mean)
#20DPI
df_final$I20DPI_mea <- apply(df_final[,15:17],1,mean)
data_final <- df_final[c("group","control_FP200000239BL","I2DPI_mea",
                         "I5DPI_mea","I10DPI_mea",
                         "I15DPI_mea","I20DPI_mea",
                         "30DPI_rep2",
                         "60DPI_rep3")]
gd1_long1<-melt(data_final,
                id.vars = c('group'),
                measure.vars = c("control_FP200000239BL","I2DPI_mea","I5DPI_mea","I10DPI_mea","I15DPI_mea","I20DPI_mea","30DPI_rep2","60DPI_rep3"),
                variable.name='time',
                value.name='ratio')

gd1_long1$time1 <- as.integer(gd1_long1$time)
p <- ggplot(gd1_long1,aes(x = time1 , y = ratio, colour = group)) + 
  geom_line(size=0.8)+  ##更改线型 线粗
  geom_point(size=2)+   ##加入点 并调整点的大小
  theme_light()+
  theme(axis.text.x=element_text(angle=60, hjust=1))+theme_bw() +
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank())
p + scale_color_manual(values= c("#FFE4E1", "#1CE6FF", "#8983BF", "#808000", "#FF4A46", "#927d85", "#A30059", "#BA0900", "#8FB0FF","#e14aec", "#63FFAC", "#FFFF00", "#71c33a", "#00868B")) +
  theme(axis.text.x=element_text(angle=60, hjust=1))


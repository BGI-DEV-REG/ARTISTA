rm(list = ls())
library(tidyr)
library(reshape2)

setwd("D:/02.project/06RYfish/RYpaper_Table_0309/meta_data/")
m = c("Batch2_Injury_2DPI_rep1_SS200000147BL_D5",
"Batch2_Injury_2DPI_rep2_SS200000147BL_D5",
"Batch2_Injury_2DPI_rep3_SS200000147BL_D4",
"Batch2_Injury_5DPI_rep1_SS200000147BL_D2",
"Batch2_Injury_5DPI_rep2_SS200000147BL_D2",
"Batch2_Injury_5DPI_rep3_SS200000147BL_D3",
"Batch2_Injury_10DPI_rep1_SS200000147BL_B5",
"Batch2_Injury_10DPI_rep2_SS200000147BL_B2",
"Batch2_Injury_10DPI_rep4_SS200000147BL_B3",
"Batch1_Injury_15DPI_rep2_FP200000266TR_E2",
"Batch1_Injury_15DPI_rep3_FP200000266TR_E3",
"Batch1_Injury_15DPI_rep4_FP200000266TR_E4",
"Batch1_Injury_15DPI_C_rep1_FP200000266TR_E6",
"Batch2_Injury_20DPI_rep1_SS200000147BL_B4",
"Batch2_Injury_20DPI_rep2_SS200000147BL_B4",
"Batch2_Injury_20DPI_rep3_SS200000147BL_B5",
"Batch1_Injury_30DPI_rep2_FP200000264BL_A6",
"Batch1_Injury_60DPI_rep3_FP200000264BL_A6",
"Batch1_Injury_control_FP200000239BL_E3")
l = list()
for (i in m) {
  data1 = read.csv(paste0(i,'_newmeta0307.txt'),sep='\t')[c('nCount_Spatial','nFeature_Spatial','Annotation_0306')]
  data1$Batch = i
  l[[i]] =data1
}

data1 = do.call(rbind,l)
m = c('Batch1_Stage44_telencephalon_rep2_FP200000239BL_E4','Batch1_Stage54_telencephalon_rep2_DP8400015649BRD6_2','Batch1_Stage57_telencephalon_rep2_DP8400015649BRD5_1','Batch1_Adult_telencephalon_rep2_DP8400015234BLA3_1','Batch1_Meta_telencephalon_rep1_DP8400015234BLB2_1')
setwd("D:/02.project/06RYfish/RYpaper_Table_0309/meta_data/")
l = list()
for (i in m) {
  data2 = read.csv(paste0(i,'.meta.csv'))[c('nCount_Spatial','nFeature_Spatial','Annotation_0305')]
  data2$Batch = i
  l[[i]] =data2
}
data2 = do.call(rbind,l)
colnames(data2) <- c("nCount_Spatial","nFeature_Spatial","Annotation_0306","Batch")
data <- rbind(data1,data2)
data$batch_id <- paste0(data$Batch,'___',data$Annotation_0306)

cellnumber_fiffantprop<- as.data.frame(prop.table(table( data$Annotation_0306,data$Batch), margin = 2))


Batch_groups<-group_by(data,batch_id)
data.fc_by_Batch<-dplyr::summarise(Batch_groups, nCount_Spatial_mean=mean(nCount_Spatial), nCount_Spatial_median=median(nCount_Spatial),nFeature_Spatial_mean=mean(nFeature_Spatial), nFeature_Spatial_median=median(nFeature_Spatial), n=n(),ra=)
data.fc_by_Batch = as.data.frame(data.fc_by_Batch)



#data.fc_by_Batch$Batch = factor(data.fc_by_Batch$Batch,levels = m)
#dat = as.data.frame(data.fc_by_Batch[order(data.fc_by_Batch$Batch),])
#write.csv(dat,'/hwfssz5/ST_SUPERCELLS/P21Z10200N0171/PROJECT/Axolotl_Brain_Regeneration/80.Spatial_add_data/25.count/03.result/Inj_19_count.csv',quote=F,row.names=F)







setwd("D:/02.project/06RYfish/RYpaper_Table_0309/DNBcount/")
m = c("Batch2_Injury_2DPI_rep1_SS200000147BL_D5",
      "Batch2_Injury_2DPI_rep2_SS200000147BL_D5",
      "Batch2_Injury_2DPI_rep3_SS200000147BL_D4",
      "Batch2_Injury_5DPI_rep1_SS200000147BL_D2",
      "Batch2_Injury_5DPI_rep2_SS200000147BL_D2",
      "Batch2_Injury_5DPI_rep3_SS200000147BL_D3",
      "Batch2_Injury_10DPI_rep1_SS200000147BL_B5",
      "Batch2_Injury_10DPI_rep2_SS200000147BL_B2",
      "Batch2_Injury_10DPI_rep4_SS200000147BL_B3",
      "Batch1_Injury_15DPI_rep2_FP200000266TR_E2",
      "Batch1_Injury_15DPI_rep3_FP200000266TR_E3",
      "Batch1_Injury_15DPI_rep4_FP200000266TR_E4",
      "Batch1_Injury_15DPI_C_rep1_FP200000266TR_E6",
      "Batch2_Injury_20DPI_rep1_SS200000147BL_B4",
      "Batch2_Injury_20DPI_rep2_SS200000147BL_B4",
      "Batch2_Injury_20DPI_rep3_SS200000147BL_B5",
      "Batch1_Injury_30DPI_rep2_FP200000264BL_A6",
      "Batch1_Injury_60DPI_rep3_FP200000264BL_A6",
      "Batch1_Injury_control_FP200000239BL_E3")
l = list()
for (i in m) {
  data1 = read.csv(paste0(i,'_DNB.count.csv'))[c('Annotation_0306','DNB_count')]
  data1$Batch = i
  l[[i]] =data1
}

data1 = do.call(rbind,l)
m = c('Batch1_Stage44_telencephalon_rep2_FP200000239BL_E4','Batch1_Stage54_telencephalon_rep2_DP8400015649BRD6_2','Batch1_Stage57_telencephalon_rep2_DP8400015649BRD5_1','Batch1_Adult_telencephalon_rep2_DP8400015234BLA3_1','Batch1_Meta_telencephalon_rep1_DP8400015234BLB2_1')
setwd("D:/02.project/06RYfish/RYpaper_Table_0309/DNBcount/")
l = list()
for (i in m) {
  data2 = read.csv(paste0(i,'_DNB.count.csv'))[c('Annotation_0305','DNB_count')]
  data2$Batch = i
  l[[i]] =data2
}
data2 = do.call(rbind,l)

colnames(data2) <- c("Annotation_0306","DNB_count","Batch")
data <- rbind(data1,data2)
data$batch_id <- paste0(data$Batch,'___',data$Annotation_0306)
Batch_groups<-group_by(data,batch_id)
DNB.fc_by_Batch<-dplyr::summarise(Batch_groups, DNB_count_mean=mean(DNB_count), DNB_count_median=median(DNB_count), n=n())
DNB.fc_by_Batch = as.data.frame(DNB.fc_by_Batch)




###celltype
m = c("Batch2_Injury_2DPI_rep1_SS200000147BL_D5",
      "Batch2_Injury_2DPI_rep2_SS200000147BL_D5",
      "Batch2_Injury_2DPI_rep3_SS200000147BL_D4",
      "Batch2_Injury_5DPI_rep1_SS200000147BL_D2",
      "Batch2_Injury_5DPI_rep2_SS200000147BL_D2",
      "Batch2_Injury_5DPI_rep3_SS200000147BL_D3",
      "Batch2_Injury_10DPI_rep1_SS200000147BL_B5",
      "Batch2_Injury_10DPI_rep2_SS200000147BL_B2",
      "Batch2_Injury_10DPI_rep4_SS200000147BL_B3",
      "Batch1_Injury_15DPI_rep2_FP200000266TR_E2",
      "Batch1_Injury_15DPI_rep3_FP200000266TR_E3",
      "Batch1_Injury_15DPI_rep4_FP200000266TR_E4",
      "Batch1_Injury_15DPI_C_rep1_FP200000266TR_E6",
      "Batch2_Injury_20DPI_rep1_SS200000147BL_B4",
      "Batch2_Injury_20DPI_rep2_SS200000147BL_B4",
      "Batch2_Injury_20DPI_rep3_SS200000147BL_B5",
      "Batch1_Injury_30DPI_rep2_FP200000264BL_A6",
      "Batch1_Injury_60DPI_rep3_FP200000264BL_A6",
      "Batch1_Injury_control_FP200000239BL_E3",'Batch1_Stage44_telencephalon_rep2_FP200000239BL_E4','Batch1_Stage54_telencephalon_rep2_DP8400015649BRD6_2','Batch1_Stage57_telencephalon_rep2_DP8400015649BRD5_1','Batch1_Adult_telencephalon_rep2_DP8400015234BLA3_1','Batch1_Meta_telencephalon_rep1_DP8400015234BLB2_1')
l = list()
setwd("D:/02.project/06RYfish/RYpaper_Table_0309/marker/")
for (i in m) {
  data = read.csv(paste0(i,'_top20_markers.xls'),sep='\t')[c('cluster','Axolotl_tanaka_annotated_gene')]
  data$Batch = i
  l[[i]] =data
}
data = do.call(rbind,l)
data$batch_id <- paste0(data$Batch,'___',data$cluster)
#Batch_groups<-group_by(data,batch_id)
ll = list()
for (i in unique(data$batch_id)) {
  me = data[which(data$batch_id==i),]
  me$seq =seq(length(me[,1]))
  ll[[i]] = me
}
da = do.call(rbind,ll)
marker.fc_by_Batch<- spread(da,seq,Axolotl_tanaka_annotated_gene)


dat1 = as.data.frame(data.fc_by_Batch[order(data.fc_by_Batch$batch_id),])
dat2 = as.data.frame(DNB.fc_by_Batch[order(DNB.fc_by_Batch$batch_id),])
dat3 = as.data.frame(marker.fc_by_Batch[order(marker.fc_by_Batch$batch_id),])

merge_data0308 <- cbind(dat1,dat2,dat3)

write.table(merge_data0308,"D:/02.project/06RYfish/RYpaper_Table_0309/meta_data/merge_data0308.xls",sep = "\t")





###celltype
m = c("Batch2_Injury_2DPI_rep1_SS200000147BL_D5",
      "Batch2_Injury_2DPI_rep2_SS200000147BL_D5",
      "Batch2_Injury_2DPI_rep3_SS200000147BL_D4",
      "Batch2_Injury_5DPI_rep1_SS200000147BL_D2",
      "Batch2_Injury_5DPI_rep2_SS200000147BL_D2",
      "Batch2_Injury_5DPI_rep3_SS200000147BL_D3",
      "Batch2_Injury_10DPI_rep1_SS200000147BL_B5",
      "Batch2_Injury_10DPI_rep2_SS200000147BL_B2",
      "Batch2_Injury_10DPI_rep4_SS200000147BL_B3",
      "Batch1_Injury_15DPI_rep2_FP200000266TR_E2",
      "Batch1_Injury_15DPI_rep3_FP200000266TR_E3",
      "Batch1_Injury_15DPI_rep4_FP200000266TR_E4",
      "Batch1_Injury_15DPI_C_rep1_FP200000266TR_E6",
      "Batch2_Injury_20DPI_rep1_SS200000147BL_B4",
      "Batch2_Injury_20DPI_rep2_SS200000147BL_B4",
      "Batch2_Injury_20DPI_rep3_SS200000147BL_B5",
      "Batch1_Injury_30DPI_rep2_FP200000264BL_A6",
      "Batch1_Injury_60DPI_rep3_FP200000264BL_A6",
      "Batch1_Injury_control_FP200000239BL_E3",'Batch1_Stage44_telencephalon_rep2_FP200000239BL_E4','Batch1_Stage54_telencephalon_rep2_DP8400015649BRD6_2','Batch1_Stage57_telencephalon_rep2_DP8400015649BRD5_1','Batch1_Adult_telencephalon_rep2_DP8400015234BLA3_1','Batch1_Meta_telencephalon_rep1_DP8400015234BLB2_1')
l = list()
setwd("D:/02.project/06RYfish/RYpaper_Table_0309/marker/")
for (i in m) {
  data = read.csv(paste0(i,'_top20_markers.xls'),sep='\t')[c('cluster','Axolotl_tanaka_annotated_gene')]
  data$Batch = i
  l[[i]] =data
}
data = do.call(rbind,l)

data$batch_id <- paste0(data$Batch,'___',data$Annotation_0306)
Batch_groups<-group_by(data,batch_id)
data.fc_by_Batch<-dplyr::summarise(Batch_groups, nCount_Spatial_mean=mean(nCount_Spatial), nCount_Spatial_median=median(nCount_Spatial),nFeature_Spatial_mean=mean(nFeature_Spatial), nFeature_Spatial_median=median(nFeature_Spatial), n=n())
data.fc_by_Batch = as.data.frame(data.fc_by_Batch)







unique(data$cluster)
data.fc_by_Batch<-dplyr::summarise(Batch_groups, nCount_Spatial_mean=mean(nCount_Spatial), nCount_Spatial_median=median(nCount_Spatial),nFeature_Spatial_mean=mean(nFeature_Spatial), nFeature_Spatial_median=median(nFeature_Spatial), n=n())
data.fc_by_Batch = as.data.frame(data.fc_by_Batch)

library(dplyr)
a <- data%>%dplyr::arrange(desc(batch_id))%>%dplyr::mutate(no=rownames(data))
ll = list()
for (i in unique(data$batch_id)) {
  me = data[which(data$batch_id==i),]
  if(length(me)>=20){
    me1 = me[c(1:20),]
    me1$seq = seq(20)
  }
  else{
    me1 =me
    me1$seq = seq(length(me1[,1]))
  }
  ll[[i]] = me1
}







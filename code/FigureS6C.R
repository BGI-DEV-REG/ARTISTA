####
library(Seurat)
library(ggplot2)
library(dplyr)

rds=readRDS("/../.rds")
ref=read.table("/../meta.csv",sep=",",header=T)

meta=rds@meta.data
meta$CellID=rownames(meta)

meta2=left_join(meta,ref,by="CellID")
rds$spatial_leiden=meta2$spatial_leiden

rds$Region=rds$spatial_leiden

rds$Region=gsub("0","VLMC",rds$Region)
rds$Region=gsub("1","Striatum",rds$Region)
rds$Region=gsub("2","Septum",rds$Region)
rds$Region=gsub("3","VZ",rds$Region)
rds$Region=gsub("4","MP",rds$Region)
rds$Region=gsub("5","LP",rds$Region)
rds$Region=gsub("6","VZ",rds$Region)
rds$Region=gsub("7","DP",rds$Region)
rds$Region=gsub("8","Septum",rds$Region)
rds$Region=gsub("9","LP",rds$Region)

data=as.data.frame(table(rds$Region))
data$all=length(rds$Region)

colnames(data)[1]="Region"
data$ratio=data$Freq/data$all
data$Batch <- factor(all$Batch,levels=c("Stage44","Stage54","Stage57","Injury_control","Adult","Meta"))

p=ggplot(data,aes(x=Batch,y=ratio,fill=Region))+
  geom_bar(stat="identity",color="black",width = .7, position = 'fill')+
  geom_text(aes(label=ratio),position=position_stack(vjust=0.5))+
  scale_fill_manual(values=c("#D51F26","#272E6A","#208A42","#89288F","#F47D2B", "#FEE500","#8A9FD1","#C06CAB","#D8A767", "#90D5E4","#89C75F","#F37B7D","#9983BD"))+
  theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour="black"))


p=ggplot(data,aes(x=Batch,y=ratio,fill=Region))+
  geom_bar(stat="identity",color="black",width = 1, position = 'fill')+
  scale_fill_manual(values=c("#EF833A","#FFDE00","#40E0D0","#238ffc","#AEC5EB", "#CCFFCC","#E71D36"))+
  theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour="black"))



#######

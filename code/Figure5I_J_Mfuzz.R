library(ggplot2)
library(Seurat)
library(ArchR)
library(Mfuzz)

rds1=readRDS("15DPI_4.rds")
rds2=readRDS("Stage57.rds")

DefaultAssay(rds1)="Spatial"
DefaultAssay(rds2)="Spatial"
rds1=NormalizeData(rds1)
rds2=NormalizeData(rds2)
rds1_sub=subset(rds1,subset = "reaEGC","rIPC1","IMN","nptxEX")
rds2_sub=subset(rd21,subset = "dEGC","dNBL4","Immature nptxEX","nptxEX")
rds1_sub$Annotation_0306=factor(rds1_sub$Annotation_0306,levels = "reaEGC","rIPC1","IMN","nptxEX")
rds2_sub$Annotation_0305=factor(rds2_sub$Annotation_0305,levels = "dEGC","dNBL4","Immature nptxEX","nptxEX")
Idents(rds1_sub)=rds1_sub$Annotation_0306
Idents(rds2_sub)=rds2_subAnnotation_0305

###Mfuzz-15DPI
ave=AverageExpression(rds1_sub)
ave_RNA=as.data.frame(ave$Spatial)
ave_sub=ave_RNA[,c(1:4)]
###
nc=4
ave_sub$Mean=apply(ave_sub[,1:nc],1,mean)
ave_sub$Sd=apply(ave_sub[,1:nc],1,sd)
ave_sub$CV=ave_sub$Sd/ave_sub$Mean*100
ave_sub$Max=apply(ave_sub[,1:nc],1,max)
ave_sub_1=subset(ave_sub,ave_sub$Max>0.5)

Set <- ExpressionSet(assayData=as.matrix(ave_sub_1[,c(1:nc)]))
Set$index=seq(0,(nc-1))
Set$label=as.factor(colnames(ave_sub_1[,c(1:nc)]))
Set$time=seq(0,(nc-1))
Set.s <- standardise(Set)
set.seed(13)
cl <- mfuzz(Set.s,c=2,m=mestimate(Set.s))
mfuzz.plot(Set.s,cl=cl,mfrow=c(2,1),time.labels=colnames(ave_sub_1[,c(1:nc)]),new.window = F)
membership=as.data.frame(cl$membership)
membership$max=apply(membership,1,max)
membership_sub=subset(membership,membership$max>0.7)
n=ncol(membership_sub)+1
a=ncol(membership_sub)
for(i in 1:nrow(membership_sub)){
  membership_sub[i,n]=paste0("Cluster",colnames(membership_sub)[which(membership_sub[i,]==membership_sub[i,a])[1]])
}
colnames(membership_sub)[n]="Cluster"
membership_sub$Gene=rownames(membership_sub)
ref=read.table("Summary_Gene_Annotaion_0707.xls",header = T,sep = "\t")
colnames(ref)[1]="Gene"
out=merge(membership_sub,ref,by="Gene")
write.table(out,"List_Mfuzz_15DPI.xls",sep = "\t",quote = FALSE,row.names = FALSE)

center <- as.data.frame(cl$centers)
center$cluster <- rownames(center)
center <- center %>%
  melt(id.var = "cluster",variable.name = "stage", value.name = "expression")
center$cluster <- paste0("Cluster",center$cluster)

Set.s1=melt(as.data.frame(Set.s))
colnames(Set.s1)=c("stage","Gene","expression")
Set.s1_sub=Set.s1

center$stage=factor(center$stage,levels = colnames(ave_sub_1)[1:nc])
Set.s1_sub$stage=factor(Set.s1_sub$stage,levels = colnames(ave_sub_1)[1:nc])
Set.s1_sub_plot=merge(Set.s1_sub,membership_sub,by="Gene")
colnames(center)[1]="Cluster"

plot=Set.s1_sub_plot[,-c(4,5,6)]
plot$Gene=as.character(plot$Gene)

plot$Cluster=factor(plot$Cluster,levels = paste0("Cluster",seq(12)))
plot$stage=factor(plot$stage,levels = c("reaEGC","rIPC1","IMN","nptxEX"))

center$Cluster=factor(center$Cluster,levels = paste0("Cluster",seq(12)))
center$stage=factor(center$stage,levels = c("reaEGC","rIPC1","IMN","nptxEX"))

p=ggplot(data=plot, aes(x=stage, y=expression, group=Gene))+
  geom_line(size = 0.05,color="#3361A5",alpha=0.5) +
  facet_wrap(~Cluster,ncol = 1)+
  geom_line(data = center,aes(stage,expression,group = 1),size=1,color = "#F7941D") + 
  geom_point(data = center,aes(stage,expression,group = 1),color = "#F7941D") +
  ylab("Expression changes") +
  xlab("") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_classic()+
  theme(axis.text.x = element_text(size = 10, vjust = 0.5, hjust = 0.5, angle = 90,color = "black"))

###Mfuzz-stage57
ave=AverageExpression(rds2_sub)
ave_RNA=as.data.frame(ave$Spatial)
ave_sub=ave_RNA[,c(1:4)]
###
nc=4
ave_sub$Mean=apply(ave_sub[,1:nc],1,mean)
ave_sub$Sd=apply(ave_sub[,1:nc],1,sd)
ave_sub$CV=ave_sub$Sd/ave_sub$Mean*100
ave_sub$Max=apply(ave_sub[,1:nc],1,max)
ave_sub_1=subset(ave_sub,ave_sub$Max>0.5)

Set <- ExpressionSet(assayData=as.matrix(ave_sub_1[,c(1:nc)]))
Set$index=seq(0,(nc-1))
Set$label=as.factor(colnames(ave_sub_1[,c(1:nc)]))
Set$time=seq(0,(nc-1))
Set.s <- standardise(Set)
set.seed(13)
cl <- mfuzz(Set.s,c=2,m=mestimate(Set.s))
mfuzz.plot(Set.s,cl=cl,mfrow=c(2,1),time.labels=colnames(ave_sub_1[,c(1:nc)]),new.window = F)
membership=as.data.frame(cl$membership)
membership$max=apply(membership,1,max)
membership_sub=subset(membership,membership$max>0.7)
n=ncol(membership_sub)+1
a=ncol(membership_sub)
for(i in 1:nrow(membership_sub)){
  membership_sub[i,n]=paste0("Cluster",colnames(membership_sub)[which(membership_sub[i,]==membership_sub[i,a])[1]])
}
colnames(membership_sub)[n]="Cluster"
membership_sub$Gene=rownames(membership_sub)
ref=read.table("Summary_Gene_Annotaion_0707.xls",header = T,sep = "\t")
colnames(ref)[1]="Gene"
out=merge(membership_sub,ref,by="Gene")
write.table(out,"List_Mfuzz_Stage57.xls",sep = "\t",quote = FALSE,row.names = FALSE)

center <- as.data.frame(cl$centers)
center$cluster <- rownames(center)
center <- center %>%
  melt(id.var = "cluster",variable.name = "stage", value.name = "expression")
center$cluster <- paste0("Cluster",center$cluster)

Set.s1=melt(as.data.frame(Set.s))
colnames(Set.s1)=c("stage","Gene","expression")
Set.s1_sub=Set.s1

center$stage=factor(center$stage,levels = colnames(ave_sub_1)[1:nc])
Set.s1_sub$stage=factor(Set.s1_sub$stage,levels = colnames(ave_sub_1)[1:nc])
Set.s1_sub_plot=merge(Set.s1_sub,membership_sub,by="Gene")
colnames(center)[1]="Cluster"

plot=Set.s1_sub_plot[,-c(4,5,6)]
plot$Gene=as.character(plot$Gene)

plot$Cluster=factor(plot$Cluster,levels = paste0("Cluster",seq(12)))
plot$stage=factor(plot$stage,levels = c("dEGC","dNBL4","Immature nptxEX","nptxEX"))

center$Cluster=factor(center$Cluster,levels = paste0("Cluster",seq(12)))
center$stage=factor(center$stage,levels = c("dEGC","dNBL4","Immature nptxEX","nptxEX"))

p=ggplot(data=plot, aes(x=stage, y=expression, group=Gene))+
  geom_line(size = 0.05,color="#3361A5",alpha=0.5) +
  facet_wrap(~Cluster,ncol = 1)+
  geom_line(data = center,aes(stage,expression,group = 1),size=1,color = "#F7941D") + 
  geom_point(data = center,aes(stage,expression,group = 1),color = "#F7941D") +
  ylab("Expression changes") +
  xlab("") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_classic()+
  theme(axis.text.x = element_text(size = 10, vjust = 0.5, hjust = 0.5, angle = 90,color = "black"))


#########Overlap Plot#########
list_15DPI=read.table("List_Mfuzz_15DPI.xls",header = T,sep = "\t")[,c(1,5)]
list_Stage57=read.table("List_Mfuzz_Stage57.xls",header = T,sep = "\t")[,c(1,5)]

list_15DPI$Group_15DPI[list_15DPI$Cluster=="Cluster1"]="Reg_UP"
list_15DPI$Group_15DPI[list_15DPI$Cluster=="Cluster2"]="Reg_Down"
list_Stage57$Group_Stage57[list_Stage57$Cluster=="Cluster1"]="Dev_UP"
list_Stage57$Group_Stage57[list_Stage57$Cluster=="Cluster2"]="Dev_Down"
mm=merge(list_15DPI[,c(1,3)],list_Stage57[,c(1,3)],by="Gene",all=T)
mm$Group_15DPI[is.na(mm$Group_15DPI)]="Reg_None"
mm$Group_Stage57[is.na(mm$Group_Stage57)]="Dev_None"
mm$Group=paste0(mm$Group_Stage57,"__",mm$Group_15DPI)

rownames(mm)=mm$Gene
mm=mm[,-1]

ave_plot=ave_SCT[match(rownames(mm),rownames(ave_SCT)),]

library(pheatmap)

pp=ave_plot[,c(5:8,1:4)]
mm$Group_Stage57=factor(mm$Group_Stage57,levels = c("Dev_UP","Dev_Down","Dev_None"))
mm$Group=factor(mm$Group,levels = c("Dev_UP__Reg_UP","Dev_UP__Reg_Down","Dev_UP__Reg_None","Dev_Down__Reg_Down","Dev_Down__Reg_UP","Dev_Down__Reg_None","Dev_None__Reg_UP","Dev_None__Reg_Down"))
mm_o=mm[order(mm$Group_Stage57,mm$Group),]
#mm_o=mm[order(mm$Group),]
pp_1=pp[match(rownames(mm_o),rownames(pp)),]

library(RColorBrewer)
getPalette = colorRampPalette(brewer.pal(12, "Paired"))
cols=getPalette(8)
names(cols)=sort(unique(mm_o$Group))

pheatmap(pp_1[,c(1:4)],scale = "row",show_rownames = FALSE,cluster_cols = FALSE,annotation_row = mm[3],
         cluster_rows = F,color = paletteContinuous(set = "solarExtra", n = 100),annotation_colors = list(Group=cols))
pheatmap(pp_1[,c(5:8)],scale = "row",show_rownames = FALSE,cluster_cols = FALSE,annotation_row = mm[3],
         cluster_rows = F,color = paletteContinuous(set = "solarExtra", n = 100),annotation_colors = list(Group=cols))

ref=read.table("Summary_Gene_Annotaion_0707.xls",header = T,sep = "\t")
colnames(ref)[1]="Gene"

mm$Gene=rownames(mm)
mm_1=merge(mm,ref,by="Gene")

list=mm_1[,c(4,7,8)]
list$mm_gene[list$mm_gene=="-"]=list$Axolotl_tanaka_annotated_gene[list$mm_gene=="-"]
list1=list[grep("[|]",list$mm_gene),]
list1$mm_gene1=unlist(lapply(strsplit(list1$mm_gene,"[|]"),"[",2))
list1$mm_gene1=unlist(lapply(strsplit(list1$mm_gene1,"[[]"),"[",1))
list1$mm_gene=list1$mm_gene1
list1=list1[,-4]
list_out=rbind(list,list1)
list_out_1=list_out[-grep("[|]",list_out$mm_gene),]

g=unique(list_out_1$Group)

for(i in seq(length(g))){
  sub=subset(list_out_1,list_out_1$Group==g[i])
  write.table(sub,paste0("CP1021_",g[i],".xls"),sep = "\t",quote = FALSE,row.names = FALSE)
}


###Line plot Stage57###
nc=4
pp_stage57=pp_1[,c(1:4)]
Set <- ExpressionSet(assayData=as.matrix(pp_stage57))
Set$index=seq(0,(nc-1))
Set$label=as.factor(colnames(pp_stage57[,c(1:nc)]))
Set$time=seq(0,(nc-1))
Set.s <- standardise(Set)
pp_stage57_out1=melt((as.data.frame(Set.s)))
pp_stage57_out1$value[is.na(pp_stage57_out1$value)]=0
colnames(pp_stage57_out1)=c("Celltype","Gene","Expression")
pp_stage57_melt=merge(pp_stage57_out1,mm,by="Gene")
pp_stage57_melt$cc=paste0(pp_stage57_melt$Celltype,":",pp_stage57_melt$Group)
pp_stage57_melt$Celltype=factor(pp_stage57_melt$Celltype,levels = c("dEGC","dNBL4","Immature nptxEX","nptxEX"))

center=as.data.frame(tapply(pp_stage57_melt$Expression,pp_stage57_melt$cc, mean))
colnames(center)[1]="Expression"
center$Celltype=unlist(lapply(strsplit(as.character(rownames(center)),":"),"[",1))
center$Group=unlist(lapply(strsplit(as.character(rownames(center)),":"),"[",2))
center$Celltype=factor(center$Celltype,levels = c("dEGC","dNBL4","Immature nptxEX","nptxEX"))

ggplot(data=pp_stage57_melt, aes(x=Celltype, y=Expression, group=Gene))+
  geom_line(size = 0.1,color="#7fbf7b",alpha=0.5) +
  facet_wrap(~Group,ncol = 2) +
  geom_line(data = center,aes(Celltype,Expression,group = 1),size=1,color = "#F7941D") + 
  geom_point(data = center,aes(Celltype,Expression,group = 1),color = "#F7941D") +
  ylab("Expression changes") +
  xlab("") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_classic()+
  theme(axis.text.x = element_text(size = 10, vjust = 0.5, hjust = 0.5, angle = 90,color = "black"))

###Line plot 15DPI###
nc=4
pp_15DPI=pp_1[,c(5:8)]
Set <- ExpressionSet(assayData=as.matrix(pp_15DPI))
Set$index=seq(0,(nc-1))
Set$label=as.factor(colnames(pp_15DPI[,c(1:nc)]))
Set$time=seq(0,(nc-1))
Set.s <- standardise(Set)
pp_15DPI_out1=melt((as.data.frame(Set.s)))
pp_15DPI_out1$value[is.na(pp_15DPI_out1$value)]=0
colnames(pp_15DPI_out1)=c("Celltype","Gene","Expression")
pp_15DPI_melt=merge(pp_15DPI_out1,mm,by="Gene")
pp_15DPI_melt$cc=paste0(pp_15DPI_melt$Celltype,":",pp_15DPI_melt$Group)
pp_15DPI_melt$Celltype=factor(pp_15DPI_melt$Celltype,levels = c("reaEGC","rIPC1","IMN","nptxEX"))

center=as.data.frame(tapply(pp_15DPI_melt$Expression,pp_15DPI_melt$cc, mean))
colnames(center)[1]="Expression"
center$Celltype=unlist(lapply(strsplit(as.character(rownames(center)),":"),"[",1))
center$Group=unlist(lapply(strsplit(as.character(rownames(center)),":"),"[",2))
center$Celltype=factor(center$Celltype,levels = c("reaEGC","rIPC1","IMN","nptxEX"))

ggplot(data=pp_15DPI_melt, aes(x=Celltype, y=Expression, group=Gene))+
  geom_line(size = 0.1,color="#7fbf7b",alpha=0.5) +
  facet_wrap(~Group,ncol = 2) +
  geom_line(data = center,aes(Celltype,Expression,group = 1),size=1,color = "#F7941D") + 
  geom_point(data = center,aes(Celltype,Expression,group = 1),color = "#F7941D") +
  ylab("Expression changes") +
  xlab("") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_classic()+
  theme(axis.text.x = element_text(size = 10, vjust = 0.5, hjust = 0.5, angle = 90,color = "black"))








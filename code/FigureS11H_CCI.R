library(Seurat)
setwd("../Anno_RDS_V2")
sample=list.files(pattern = "Anno_Batch2_")
sample1=gsub("Anno_Batch2_|.rds","",sample)

for (i in seq(length(sample))) {
  
  path_rds=paste0("../Anno_RDS_V2/",sample[i])
  path_meta=paste0("../Meta_Sptaial/",sample1[i],".meta.txt")
  rds=readRDS(path_rds)
  meta2=read.table(path_meta)
  meta2=meta2[match(colnames(rds),rownames(meta2)),]
  meta1=rds@meta.data
  name=gsub("Anno_|.rds","",sample[i])
  
  rds=subset(rds,cells = )
  
  rds$Spatial_Cluster=paste0("S_",meta2$spatial_leiden_e30_s8)
  Idents(rds)=rds$Spatial_Cluster
  ave=AverageExpression(rds)
  
  ave_sct=as.data.frame(ave$SCT)
  lr_ref=read.table("human_lr_pair.txt",header = T)
  lr_ref1=lr_ref[,c(2,8)]
  lr_ref2=lr_ref[,c(3,9)]
  colnames(lr_ref1)[2]="hs_ID"
  colnames(lr_ref2)[2]="hs_ID"
  colnames(lr_ref1)[1]="hs_Gene"
  colnames(lr_ref2)[1]="hs_Gene"
  lr_ref1$Num_Ligand=seq(nrow(lr_ref1))
  lr_ref2$Num_Receptor=seq(nrow(lr_ref2))
  
  gene=read.table("/Summary_Gene_Annotaion_0707V2.xls",header = T)
  ave_sct$Axolotl_ID=rownames(ave_sct)
  ave_sct_sub=merge(ave_sct,gene,by="Axolotl_ID")
  
  ave_sct_sub_ligand=merge(ave_sct_sub,lr_ref1,by="hs_ID")
  ave_sct_sub_receptor=merge(ave_sct_sub,lr_ref2,by="hs_ID")
  
  all=sort(intersect(ave_sct_sub_ligand$Num_Ligand,ave_sct_sub_receptor$Num_Receptor))
  
  ave_sct_sub_ligand_u=ave_sct_sub_ligand[match(all,ave_sct_sub_ligand$Num_Ligand),]
  ave_sct_sub_receptor_u=ave_sct_sub_receptor[match(all,ave_sct_sub_receptor$Num_Receptor),]
  
  n=length(unique(rds$Spatial_Cluster))
  
  mx_ligand=ave_sct_sub_ligand_u[,c(3:(2+n))]
  mx_receptor=ave_sct_sub_receptor_u[,c(3:(2+n))]
  sum_lr=mx_ligand*mx_receptor
  
  ratio=as.data.frame(t(apply(sum_lr,1,function(x) x/sum(x))))
  ratio[ratio==0]=0.000001
  entropy <- function(a){
    b<- -sum(a*log2(a))
    return(b)
  }
  entro=as.data.frame(apply(ratio,1,entropy))
  
  ratio$entropy=entro[,1]
  ratio$ID_Ligand=ave_sct_sub_ligand_u$Axolotl_ID
  ratio$ID_Receptor=ave_sct_sub_receptor_u$Axolotl_ID
  
  ratio$LR=paste0(ave_sct_sub_ligand_u$hs_Gene,"_",ave_sct_sub_receptor_u$hs_Gene)
  ratio=na.omit(ratio)
  ratio$CL="NA"
  for (i in seq(nrow(ratio))) {
    ratio[i,ncol(ratio)]=colnames(ratio)[which(ratio[i,1:n]==max(ratio[i,1:n]))]
  }
  ratio_sub=subset(ratio,ratio$entropy<2)
  rownames(ratio_sub)=ratio_sub$LR
  
  out1=paste0("CCI_Ratio_",name,".txt")
  write.table(ratio_sub,out1,sep = "\t",quote = FALSE,row.names = FALSE)
  
  out3=paste0("CCI_Ratio_ALL_",name,".txt")
  write.table(ratio,out3,sep = "\t",quote = FALSE,row.names = FALSE)
  
  
  library(pheatmap)
  out2=paste0("CCI_Heatmap_",name,".pdf")
  p=pheatmap(ratio_sub[,1:n],fontsize = 5,treeheight_row = 0,treeheight_col = 0)
  pdf(out2,width = 8,height = 15)
  print(p)
  dev.off()
}



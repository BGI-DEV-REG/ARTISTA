library(motifmatchr)

all=read.csv("Time_diff_all.csv")
all_sub=subset(all,all$qval<0.01)

###
tf=readxl::read_xlsx("Human_TF_MotifList_v_1.01 (2).xlsx")
tf_gene=unique(tf$`HGNC symbol`)
ref=read.table("Summary_Gene_Annotaion_0707.xls",header = T,sep = "\t")

tab4=all_sub[5]
colnames(tab4)[1]="Axolotl_ID"
tab4_ref=merge(tab4,ref,by="Axolotl_ID")

out3=list()
for (i in seq(length(tf_gene))) {
  g1=grep(tf_gene[i],tab4_ref$Axolotl_tanaka_annotated_gene,ignore.case = T)
  g2=grep(tf_gene[i],tab4_ref$mm_gene,ignore.case = T)
  g3=grep(tf_gene[i],tab4_ref$hs_gene,ignore.case = T)
  g4=grep(tf_gene[i],tab4_ref$Human_Annotation,ignore.case = T)
  g5=grep(tf_gene[i],tab4_ref$Mouse_Annotation,ignore.case = T)
  
  g=unique(c(g1,g2,g3,g4,g5))
  
  if(length(g)>0){
    print(length(g))
    print(g)
    print(i)
    print(tf_gene[i])
    o=tab4_ref[g,]
    o$TF_Gene=tf_gene[i]
    out3[[i]]=o
  }
}
out3_data=do.call(rbind,out3)
tf_selected=c("BHLHE22","DLX6","EMX1","HES5","NEUROD2","NEUROD4","NEUROD6","NFIB","NR2E1","SOX2","SOX8","STAT5B","ZBTB18","ZBTB20","ZEB2","ZIC1","ZNF484")
###---------motif match-----------###
gene=read.table("tss_geneid.txt",header = T)
seq=read.table("gene_tss_2k_seq.bed")

gene$Tss_start_2k=gene$Tss_start_2k+1
gene$ID=paste(gene$Chr,gene$Tss_start_2k,gene$Tss_end_2k,sep = "_")
gene=gene[,c(4,7)]
seq$V1=gsub(">","",seq$V1)
seq$V1=gsub(":.","",seq$V1)
seq$V1=gsub("-","_",seq$V1)
colnames(seq)[1]="ID"
mm=merge(seq,gene,by="ID")

list_sub=all_sub
colnames(list_sub)[5]="Gene"
mm1=merge(list_sub,mm,by="Gene",all.x = T)
mm2=merge(list_sub,mm,by="Gene")

tf_used=tf[,c(2,7,8)]
tf_selected=as.data.frame(tf_selected)
colnames(tf_used)[1]="Gene"
colnames(tf_selected)[1]="Gene"

tt=merge(tf_used,tf_selected,by="Gene")
tt_sub=subset(tt,tt$`Best Motif(s)? (Figure 2A)`=="TRUE")[-2,]
#tt_sub[5,2]="M03671_1.94d"

motif_wxy=list()
for(i in seq(nrow(tt_sub))){
  tf_test=read.table(paste0("../PWMs/",tt_sub[i,2],".txt"),header = T)
  tf_test1=t(tf_test[,-1])
  tf_test2=tf_test1*25
  pfm <- PFMatrix(ID=tt_sub[i,2], name=tt_sub[i,1],strand="+", bg=c(A=0.25, C=0.25, G=0.25, T=0.25), profileMatrix=tf_test2)
  motif_wxy[[i]]=pfm
}
names(motif_wxy)=tt_sub$Gene

pfmList <- do.call(PFMatrixList, motif_wxy)

sum=list()
for (i in seq(nrow(mm2))) {
  motif_ix_1 <- matchMotifs(pfmList, as.character(mm2[i,9]))
  r1=t(as.matrix(motifMatches(motif_ix_1)))
  out=as.data.frame(r1)
  out$TF=rownames(out)
  out$Gene=mm2[i,1]
  sum[[i]]=out
  print(paste0(i,"  finished!!!"))
}
sum_all=do.call(rbind,sum)
sum_all_sub=subset(sum_all,sum_all$V1=="TRUE")
colnames(sum_all_sub)[3]="Axolotl_ID"
tab4_used_plot=merge(sum_all_sub,tab4_ref,by="Axolotl_ID")
write.table(tab4_used_plot,"Summary_reaEGC_TF_network.xls",sep = "\t",quote = FALSE,row.names = FALSE)



















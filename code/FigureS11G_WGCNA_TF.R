library(WGCNA)
library(Seurat)
library(TFBSTools)
library(motifmatchr)

rds=readRDS("Batch2_Injury_2DPI_rep1_SS200000147BL_D5.rds")
DefaultAssay(rds)
Idents(rds)=rds$Annotation_0306

rds_sub=subset(rds,idents = c("reaEGC","ribEGC","sfrpEGC","wntEGC"))
mk=FindAllMarkers(rds_sub)

mk_1=FindMarkers(rds,ident.1 = "reaEGC")
mk_2=FindMarkers(rds,ident.1 = "sfrpEGC")
mk_3=FindMarkers(rds,ident.1 = "wntEGC")
mk_4=FindMarkers(rds,ident.1 = "ribEGC")

mk_1_sub=subset(mk_1,mk_1$avg_log2FC>0.3 & mk_1$p_val_adj<0.05)
mk_2_sub=subset(mk_2,mk_2$avg_log2FC>0.3 & mk_2$p_val_adj<0.05)
mk_3_sub=subset(mk_3,mk_3$avg_log2FC>0.3 & mk_3$p_val_adj<0.05)
mk_4_sub=subset(mk_4,mk_4$avg_log2FC>0.3 & mk_4$p_val_adj<0.05)
mk_sub=rbind(mk_1_sub,mk_2_sub,mk_3_sub,mk_4_sub)
mk_sub$gene=rownames(mk_sub)

exp=as.data.frame(rds_sub@assays$SCT@data)
exp_sub=exp[match(unique(mk_sub$gene),rownames(exp)),]

meta=rds_sub@meta.data
meta_sub=meta[8]
colnames(meta_sub)="subtype"
meta_sub$gsm=rownames(meta_sub)
meta_sub=meta_sub[,c(2,1)]

###merge cells
pseudocell.size = 10
meta_sub_1=subset(meta_sub,meta_sub$subtype=="reaEGC")
meta_sub_2=subset(meta_sub,meta_sub$subtype=="sfrpEGC")
meta_sub_3=subset(meta_sub,meta_sub$subtype=="ribEGC")
meta_sub_4=subset(meta_sub,meta_sub$subtype=="wntEGC")

meta_sub_1$num=seq(nrow(meta_sub_1))
meta_sub_2$num=seq(nrow(meta_sub_2))
meta_sub_3$num=seq(nrow(meta_sub_3))
meta_sub_4$num=seq(nrow(meta_sub_4))

meta_sub_1$ID=paste0("reaEGC__",meta_sub_1$num%/%pseudocell.size)
meta_sub_2$ID=paste0("sfrpEGC__",meta_sub_2$num%/%pseudocell.size)
meta_sub_3$ID=paste0("ribEGC__",meta_sub_3$num%/%pseudocell.size)
meta_sub_4$ID=paste0("wntEGC__",meta_sub_4$num%/%pseudocell.size)

meta_sub_mm=rbind(meta_sub_1,meta_sub_2,meta_sub_3,meta_sub_4)
exp_sub_t=as.data.frame(t(exp_sub))
exp_sub_t1=cbind(meta_sub_mm[4],exp_sub_t)
exp_sub_t_mean=aggregate(exp_sub_t1,list(exp_sub_t1$ID),mean)
rownames(exp_sub_t_mean)=exp_sub_t_mean$Group.1
exp_sub_t_mean=exp_sub_t_mean[,-c(1,2)]

###WGCNA 

datExpr=exp_sub_t_mean
datTraits=unique(meta_sub_mm[,c(4,2)])
colnames(datTraits)[1]="gsm"

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# Plot the results:
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="green")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
#choose beta score
best_beta=sft$powerEstimate
best_beta
#blockwiseConsensusModules
multiExpr <- list()
multiExpr[['REAEGC']] <- list(data=datExpr)
#softPower=best_beta
softPower=best_beta
net=blockwiseConsensusModules(multiExpr, blocks = NULL,
                              maxBlockSize = 2000, ## This should be set to a smaller size if the user has limited RAM
                              randomSeed = 12345,
                              corType = "pearson",
                              power = softPower,
                              consensusQuantile = 0.3,
                              networkType = "signed",
                              TOMType = "unsigned",
                              TOMDenom = "min",
                              scaleTOMs = TRUE, scaleQuantile = 0.8,
                              sampleForScaling = TRUE, sampleForScalingFactor = 1000,
                              useDiskCache = TRUE, chunkSize = NULL,
                              deepSplit = 4,
                              pamStage=FALSE,
                              detectCutHeight = 0.995, minModuleSize = 20,
                              mergeCutHeight = 0.25)

table(net$colors)
# Convert labels to colors for plotting
#mergedColors = labels2colors(net$colors)
mergedColors = net$colors

table(mergedColors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
#dev.off()
saveRDS(net,"NET_WGCNA_2DPI.RDS")
#########
nSamples = nrow(datExpr)
design=model.matrix(~0+ datTraits$subtype)
datTraits$subtype=factor(datTraits$subtype)
colnames(design)=levels(datTraits$subtype)
moduleColors <- net$colors
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0); 
moduleTraitCor = cor(MEs, design , use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
#pdf("Nor_20000_labeledHeatmap4.pdf",width = 10,height = 8)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(design),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
#dev.off()

list=as.data.frame(net$colors)
list$Gene=rownames(list)
lits_sub=subset(list,list$`net$colors`=="brown")
lits_sub$Gene=rownames(lits_sub)
write.table(lits_sub,"WGCNA_2DPI_REAEGC(brown).xls",sep = "\t",quote = FALSE)

###-----------TF Analysis------------###
tf=readxl::read_xlsx("Human_TF_MotifList_v_1.01 (2).xlsx")
tf_gene=unique(tf$`HGNC symbol`)
ref=read.table("Summary_Gene_Annotaion_0707.xls",header = T,sep = "\t")

tab4=lits_sub[2]
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

tf_selected=c("CEBPD","EGR1","JUN","JUNB","RUNX1","ZNF585A")

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
mm1=merge(lits_sub,mm,by="Gene",all.x = T)
mm2=merge(lits_sub,mm,by="Gene")

tf_used=tf[,c(2,7,8)]
tf_selected=as.data.frame(tf_selected)
colnames(tf_used)[1]="Gene"
colnames(tf_selected)[1]="Gene"

tt=merge(tf_used,tf_selected,by="Gene")
tt_sub=subset(tt,tt$`Best Motif(s)? (Figure 2A)`=="TRUE")[-2,]
tt_sub[5,2]="M03671_1.94d"

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
  motif_ix_1 <- matchMotifs(pfmList, as.character(mm2[i,4]))
  r1=t(as.matrix(motifMatches(motif_ix_1)))
  out=as.data.frame(r1)
  out$TF=rownames(out)
  out$Gene=mm2[i,1]
  sum[[i]]=out
}
sum_all=do.call(rbind,sum)
sum_all_sub=subset(sum_all,sum_all$V1=="TRUE")
colnames(sum_all_sub)[3]="Axolotl_ID"
tab4_used_plot=merge(sum_all_sub,tab4_ref,by="Axolotl_ID")
write.table(tab4_used_plot,"Summary_WGCNA_REAEGC_brown.xls",sep = "\t",quote = FALSE,row.names = FALSE)
tab4_used_plot1=tab4_used_plot[-grep("AMEX",tab4_used_plot$Axolotl_tanaka_annotated_gene),]
write.table(tab4_used_plot1[,c(3,6)],"Plot_WGCNA_REAEGC_brown.txt",sep = "\t",quote = FALSE,row.names = FALSE)
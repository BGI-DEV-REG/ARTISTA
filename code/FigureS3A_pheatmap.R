library(Seurat)
library(ggplot2)
library(reshape2)
library(dplyr)
library(pheatmap)
library(RColorBrewer)
library(viridis)

setwd('../FigureS3_data/FigureS3A')
### Adult gene annotation
ref = read.csv('Adult_annotation.txt',header = F,sep = '\t')
### AverageExpression
sct <- read.csv('Batch1_Adult_telencephalon_rep2_DP8400015234BLA3_1.SCT_Removed_AverageExpression_0305.csv',header = T,sep = ',')
### genes order
gene = read.csv('Adult.0826.mat_cluster.csv')
d = as.data.frame(colnames(sct))
colnames(d) = 'order'

l =c('SSTIN','TLNBL','MSN','CCKIN','MPIN','DPEX','MPEX','NPTXEX','CMPN','SCGNIN','NTNG1EX','CP','VLMC','RIBEGC','SFRPEGC','WNTEGC')
setdiff(l,d$order)

d$order = factor(d$order,levels = rev(l))
d = as.data.frame(d[order(d$order),])
colnames(d) = 'order'

da = merge(ref,sct,by.x = 'V1',by.y = 'X',sort=F) 
gene = as.data.frame(gene$X)
colnames(gene) = 'gene'
da = merge(gene,da,by.x = 'gene',by.y = 'V2',sort = F)
colnames(da)
rownames(da) = da$gene

d = as.data.frame(d[c(1:16),])
colnames(d) = 'order'

da_order_1=da[,match(d$order,colnames(da))]
co <- c('#3361A5','#3164AB','#3067B1','#2F6AB7','#2E6EBE','#2C71C4','#2B74CA','#2A78D1','#297BD7','#287EDD','#2682E4','#2585EA','#2488F0','#238CF3','#218FF4','#2092F5','#1F96F6','#1E99F7','#1C9CF8','#1B9FF9','#1AA3FA','#18A6FB','#17A9FC','#16ADFD','#14B0FE','#16B3FE','#1FB5FD','#29B7FC','#32BAFA','#3BBCF9','#45BEF8','#4EC0F6','#57C2F5','#61C5F4','#6AC7F3','#74C9F1','#7DCBF0','#86CDEF','#8CCEED','#90CFEC','#95CFEA','#99D0E9','#9ED0E7','#A3D1E5','#A7D1E4','#ACD2E2','#B0D3E1','#B5D3DF','#BAD4DE','#BED4DC','#C2D4D9','#C5D4D3','#C9D4CE','#CCD4C8','#CFD4C2','#D3D4BD','#D6D3B7','#D9D3B2','#DDD3AC','#E0D3A7','#E3D3A1','#E7D39B','#EAD295','#EBD08B','#EDCD81','#EECA77','#F0C86D','#F1C563','#F3C359','#F4C04F','#F6BD44','#F8BB3A','#F9B830','#FBB626','#FCB31C','#FBAA1A','#F99F1C','#F7941D','#F5891E','#F37E20','#F17321','#EF6822','#ED5D24','#EB5225','#E94726','#E73B27','#E53029','#E22929','#DC2828','#D72727','#D22626','#CD2525','#C72424','#C22323','#BD2222','#B82121','#B22020','#AD1F1F','#A81E1E','#A31D1D',"#800026","#800026","#800026","#800026","#800026","#800026","#800026")

p = pheatmap(da_order_1,cluster_cols = F,cluster_rows = F,scale = "row",color =co)
ggsave(plot=p,file = 'Adult_pheatmap.pdf')

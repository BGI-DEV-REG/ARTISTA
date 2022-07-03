library(Seurat)
rds=readRDS("/../.rds")
meta=read.table("/../meta_Region.csv",sep=",",header=T)
rds$Region=meta$Region
Idents(rds)=rds$Region
DefaultAssay(rds)="SCT"
marker=FindAllMarkers(rds,group.by="Region")#,logfc.threshold = 0.25,only.pos=T)

cluster.averages <- AverageExpression(rds)
data=as.data.frame(cluster.averages$SCT)

for(i in c(1:length(colnames(data)))){
  colnames(data)[i]=paste0("Adult_",colnames(data)[i])
}

gene_exist=intersect(gene_exist,e)

data=data[match(unique(gene_exist),rownames(data)),]

data_all=cbind.data.frame(data_Stage44,data_Stage54,data_Stage57,data_Injury_control,data_Meta)


matrix <- cor(data_all,data_all,method="spearman",use="complete.obs")

co <- c('#3361A5','#3361A5','#3361A5','#3361A5','#3361A5','#3361A5','#3361A5','#3361A5','#3361A5','#3361A5','#3361A5','#3361A5','#3361A5',
        '#3361A5','#3361A5','#3361A5','#3361A5','#3361A5','#3361A5','#3361A5','#3164AB','#3067B1','#2F6AB7','#2E6EBE','#2C71C4','#2B74CA','#2A78D1',
        '#297BD7','#287EDD','#2682E4','#2585EA','#2488F0','#238CF3','#218FF4','#2092F5','#1F96F6','#1E99F7','#1C9CF8','#1B9FF9','#1AA3FA','#18A6FB',
        '#17A9FC','#16ADFD','#14B0FE','#16B3FE','#1FB5FD','#29B7FC','#32BAFA','#3BBCF9','#45BEF8','#4EC0F6','#57C2F5','#61C5F4','#6AC7F3','#74C9F1',
        '#7DCBF0','#86CDEF','#8CCEED','#90CFEC','#95CFEA','#99D0E9','#9ED0E7','#A3D1E5','#A7D1E4','#ACD2E2','#B0D3E1','#B5D3DF','#BAD4DE','#BED4DC',
        '#C2D4D9','#C5D4D3','#C9D4CE','#CCD4C8','#CFD4C2','#D3D4BD','#D6D3B7','#D9D3B2','#DDD3AC','#E0D3A7','#E3D3A1','#E7D39B','#EAD295','#EBD08B',
        '#EDCD81','#EECA77','#F0C86D','#F1C563','#F3C359','#F4C04F','#F6BD44','#F8BB3A','#F9B830','#FBB626','#FCB31C','#FBAA1A','#F99F1C','#F7941D',
        '#F5891E','#F37E20','#F17321','#EF6822','#ED5D24','#EB5225','#E94726','#E73B27','#E53029','#E22929','#DC2828','#D72727','#D22626','#CD2525',
        '#C72424','#C22323','#BD2222','#B82121','#B22020','#AD1F1F','#A81E1E','#A31D1D','#A31D1D','#A31D1D','#A31D1D','#A31D1D','#A31D1D','#A31D1D',
        '#A31D1D','#A31D1D','#A31D1D','#A31D1D','#A31D1D','#A31D1D','#A31D1D','#A31D1D','#A31D1D','#A31D1D','#A31D1D','#A31D1D')
##黑色'#000000','#000000','#000000','#000000','#000000','#000000',
####
mat <- dist(matrix)
hclust_mat <- hclust(mat)
hclust_mat$order
hclust_mat$labels

index=c(5,9,15,27,41,34,12,16,23,39,33,3,8,18,22,40,32,6,10,17,26,38,31,4,11,19,25,37,29,1,13,20,24,36,30,2,7,14,21,35,28)
hclust_mat$order <- index

y=c("MP.DP_Stage44","DP_Stage54","DP_Stage57","DP_Injury_control","DP_Adult","DP_Meta",
                    "MP_Stage54","MP_Stage57","MP_Injury_control","MP_Adult","MP_Meta",
                    "LP_Stage44","LP_Stage54","LP_Stage57","LP_Injury_control","LP_Adult","LP_Meta",
                    "VZ_Stage44","VZ_Stage54","VZ_Stage57","VZ_Injury_control","VZ_Adult","VZ_Meta",
                    "Striatum_Stage44","Striatum_Stage54","Striatum_Stage57","Striatum_Injury_control","Striatum_Adult","Striatum_Meta",
                    "Septum_Stage44","Septum_Stage54","Septum_Stage57","Septum_Injury_control","Septum_Adult","Septum_Meta",
                   "VLMC_Stage44","VLMC_Stage54","VLMC_Stage57","VLMC_Injury_control","VLMC_Adult","VLMC_Meta")
hclust_mat$labels=y

gn=rownames(matrix)[p$tree_row[["order"]]]
sn=colnames(matrix)[p$tree_col[["order"]]]
new_matrix=matrix[gn,sn]



pheatmap(new_matrix,color =co,cluster_rows = F,cluster_cols = F)#, cluster_rows = hclust_mat,cluster_cols =hclust_mat)#,cluster_rows = F,cluster_cols = F)

####

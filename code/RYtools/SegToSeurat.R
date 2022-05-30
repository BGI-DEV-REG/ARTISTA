# segmented data to Seurat SpatailObject
# Author: Feng Weimin
parser = argparse::ArgumentParser(description="Segmented data to Seurat SpatailObject, saving as .rds and .h5ad")
parser$add_argument('-I','--input', help='input directory')
parser$add_argument('-D','--id', help='tissue ID')
parser$add_argument('-O','--output', help='out directory')
parser$add_argument('-N','--nCount_Spatial', help='nCount_RNA')
args = parser$parse_args()

library(dplyr)
library(data.table)
library(Matrix)
library(Seurat)
library(ggplot2)
library(SeuratObject)
library(SeuratDisk)

source('LoadBGI_Spatial.R')

id = args$id
wd = paste0(args$input)

setwd(wd)

#### Creat obj
obj = LoadBGI_Spatial(paste0(id,'_scgem.csv.gz'),
                             outdir = getwd(),
                             bin_data = F,
                             bin_size = 50,
                             cell_mask = T,
                             area_mask = F,
                             save_as = "rds",
                             pro_name = id,
                             UMI_GreyScale_Image = F,
                             assay = "Spatial",
                             slice = id,
                             delete_bg = T,csv=F)



#### QC and selecting cells for further analysis
obj <- obj[,obj$nCount_Spatial > as.numeric(args$nCount_Spatial)] # nCount_RNA can be adjusted

saveRDS(obj,paste0(args$output,'/',args$id,'.rds'))


#### Convert as h5ad file
library(reticulate)

path_to_python = "" # replace with your python path
use_python(path_to_python)

genes <- as.data.frame(rownames(obj), row.names = rownames(obj))
names(genes) <- "Gene"

cells <- as.data.frame(colnames(obj), row.names = colnames(obj))
names(cells) <- "CellID"

row <- obj@images[[1]]@coordinates$row
col <- obj@images[[1]]@coordinates$col
coordinates <- list(matrix(c(row, col), ncol = 2))
names(coordinates) <- "spatial"

ann <- import("anndata")
ad <- ann$AnnData(X = obj@assays$Spatial@data, obs = genes, var = cells, varm = coordinates,layers = list(counts = obj@assays$Spatial@counts))
ad <- ad$T

ad$write_h5ad(file.path(args$output, paste0(args$id, ".h5ad")))



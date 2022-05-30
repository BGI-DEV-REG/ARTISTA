########
# modified: median --> mean

LoadBGI_Spatial <- function(filename,
                            outdir = getwd(),
                            bin_data = T,
                            bin_size = 50,
                            cell_mask = F,
                            area_mask = F,
                            save_as = "rds",
                            pro_name = "Spatial",
                            UMI_GreyScale_Image = T,
                            assay = "Spatial",
                            slice = "slice1",
			    delete_bg = T,csv=F) {
    if (csv)
	{dat<-read.csv(file=filename)
  print('mark')}
    else
    {dat <- data.table::fread(file = filename)}
    colnames(dat)[1:4] <- c("geneID", "x", "y", "UMICount")
    if (bin_data) {
        rowname_style <- "BIN"

        dat$x <- trunc((dat$x - min(dat$x)) / bin_size + 1)
        dat$y <- trunc((dat$y - min(dat$y)) / bin_size + 1)

        if (area_mask) {
            dat$bin_ID <- max(dat$x) * (dat$y - 1) + dat$x
            Area_ <- dat$AreaID[!duplicated(dat$bin_ID)]
        }

        if ("MIDCounts" %in% colnames(dat)) {
            dat <- dat[, sum(MIDCounts), by = .(geneID, x, y)]
        } else {
            dat <- dat[, sum(UMICount), by = .(geneID, x, y)]
        }
        dat$bin_ID <- max(dat$x) * (dat$y - 1) + dat$x
        bin.coor <- dat[, sum(V1), by = .(x, y)]
        if (UMI_GreyScale_Image) {
            scale_grey <- bin.coor$V1 / quantile(bin.coor$V1, probs = 0.95)
            scale_grey[scale_grey > 1] <- 1
            tissue_lowres_image <- as.matrix(Matrix::sparseMatrix(
            i = bin.coor$y,
            j = bin.coor$x,
            x = scale_grey
            ))
            tissue_lowres_image_r <- raster::raster(tissue_lowres_image)
            tissue_lowres_image_r <- raster::writeRaster(tissue_lowres_image_r, file.path(outdir, paste0(pro_name, "_UMIGreyScaleFakeIMG.tiff")), overwrite=T)
        } else tissue_lowres_image <- matrix(0, max(bin.coor$y), max(bin.coor$x))
    }

    if (cell_mask) {
        
        rowname_style <- "CELL"
        colnames(dat)[1:5] <- c("geneID", "x", "y", "UMICount", "label")
        if (delete_bg){
	dat<-subset(dat,label>0)}
        dat$x <- dat$x - min(dat$x) + 1
        dat$y <- dat$y - min(dat$y) + 1

        dat_bkp <- dat
        dat_bkp <- dat_bkp[, sum(UMICount), by = .(geneID, x, y)]
        bin.coor.bkp <- dat_bkp[, sum(V1), by = .(x, y)]

        if (UMI_GreyScale_Image) {
            scale_grey <- bin.coor.bkp$V1 / quantile(bin.coor.bkp$V1, probs = 0.95)
            scale_grey[scale_grey > 1] <- 1
            tissue_lowres_image <- as.matrix(Matrix::sparseMatrix(
                i = bin.coor.bkp$y,
                j = bin.coor.bkp$x,
                x = scale_grey
            ))
            tissue_lowres_image_r <- raster::raster(tissue_lowres_image)
            tissue_lowres_image_r <- raster::writeRaster(tissue_lowres_image_r, file.path(outdir, paste0(pro_name, "_UMIGreyScaleFakeIMG.tiff")), overwrite=T)
        } else tissue_lowres_image <- matrix(1, max(bin.coor.bkp$y), max(bin.coor.bkp$x))

        dat <- dat[dat$label != 0,]
        dat.x <- as.data.frame(dat[, ceiling(mean(x)), by = .(label)])
        hash.x <- data.frame(row.names = dat.x$label, values = dat.x$V1)
        dat.y <- as.data.frame(dat[, ceiling(mean(y)), by = .(label)])
        hash.y <- data.frame(row.names = dat.y$label, values = dat.y$V1)
        dat$x <- hash.x[sprintf("%d", dat$label), "values"]
        dat$y <- hash.y[sprintf("%d", dat$label), "values"]

        te <- data.frame(unique(paste(dat$x, dat$y, dat$label, sep = "_"))) # slow
        colnames(te) <- "xyb"
        split_b <- stringr::str_split(te$xyb, "_")
        te$x <- as.vector(sapply(split_b, "[", 1))
        te$y <- as.vector(sapply(split_b, "[", 2))
        te$b <- as.vector(sapply(split_b, "[", 3))
        te$xy <- paste(te$x, te$y, sep = "_")
        dp <- which(duplicated(te$xy))
        dplen <- length(dp)

        for (i in dplen) {
            tmp <- dp[i]
            xy <- te[tmp, ]$xy
            cxy <- which(te$xy %in% xy)
            cxylen <- length(cxy)
            keep <- te[cxy[1], ]$b
            for (m in 2:cxylen) {
                c <- which(dat$label %in% te[cxy[m], ]$b)
                dat$label[c] <- as.numeric(keep)
            }
        }
        dat$label <- as.numeric(dat$label)
        dat <- dat[, sum(UMICount), by = .(geneID, x, y,label)]
        dat$bin_ID <- dat$label
        bin.coor <- dat[, sum(V1), by = .(x, y)]
    }

    geneID <- seq(length(unique(dat$geneID)))
    hash.G <- data.frame(row.names = unique(dat$geneID), values = geneID)
    gen <- hash.G[dat$geneID, "values"]
    bin_ID <- unique(dat$bin_ID)
    hash.B <- data.frame(row.names = sprintf("%d", bin_ID), values = bin_ID)
    bin <- hash.B[sprintf("%d", dat$bin_ID), "values"]
    cnt <- dat$V1
    mat <- Matrix::sparseMatrix(i = gen, j = bin, x = cnt)
    rownames(mat) <- rownames(hash.G)
    colnames(mat) <- paste(rowname_style, sprintf("%d", seq(max(hash.B[, "values"]))), sep = ".")
    if (area_mask) {
        hash.A <- data.frame(row.names = paste(rowname_style, rownames(hash.B), sep = '.'), area = Area_)
        object <- Seurat::CreateSeuratObject(mat, project = assay, assay = assay, meta.data = hash.A)
    } else object <- Seurat::CreateSeuratObject(mat, project = assay, assay = assay)
    rm(dat)
    gc()

    generate_spatialObj <- function(image, scale.factors, tissue.positions, filter.matrix = TRUE) {
        if (filter.matrix) {
            tissue.positions <- tissue.positions[which(tissue.positions$tissue == 1), , drop = FALSE]
        }
        unnormalized.radius <- scale.factors$fiducial_diameter_fullres * scale.factors$tissue_lowres_scalef
        spot.radius <- unnormalized.radius / max(dim(image))
        return(new(Class = "VisiumV1",
                    image = image,
                    scale.factors = Seurat::scalefactors(
                        spot = scale.factors$tissue_hires_scalef,
                        fiducial = scale.factors$fiducial_diameter_fullres,
                        hires = scale.factors$tissue_hires_scalef,
                        lowres = scale.factors$tissue_lowres_scalef
                    ),
                    coordinates = tissue.positions,
                    spot.radius = spot.radius
                ))
    }

    tissue_positions_list <- data.frame(
        row.names = paste(rowname_style, rownames(hash.B), sep = "."),
        tissue = 1,
        row = bin.coor$y, col = bin.coor$x,
        imagerow = bin.coor$y, imagecol = bin.coor$x
    )
    scalefactors_json <- rjson::toJSON(list(
        fiducial_diameter_fullres = 1,
        tissue_hires_scalef = 1,
        tissue_lowres_scalef = 1
    ))
    spatialObj <- generate_spatialObj(
        image = tissue_lowres_image,
        scale.factors = rjson::fromJSON(scalefactors_json),
        tissue.positions = tissue_positions_list
    )

    spatialObj <- spatialObj[Seurat::Cells(object)]
    Seurat::DefaultAssay(spatialObj) <- assay
    object[[slice]] <- spatialObj

    object <- subset(object, subset = nCount_Spatial > 0)

    if (save_as == "rds") {
        saveRDS(object, file = file.path(outdir, paste0(pro_name, "_spatialObj.rds")))
    }
    if (save_as == "h5ad") {
        genes <- as.data.frame(rownames(object), row.names = rownames(object))
        names(genes) <- "Gene"

        cells <- as.data.frame(colnames(object), row.names = colnames(object))
        names(cells) <- "CellID"

        row <- object@images[[slice]]@coordinates$row
        col <- object@images[[slice]]@coordinates$col
        coordinates <- list(matrix(c(row, col), ncol = 2))
        names(coordinates) <- "spatial"

        ad <- anndata::AnnData(X = object@assays$Spatial@counts, obs = genes, var = cells, varm = coordinates)
        ad <- ad$T

        ad$write_h5ad(file.path(outdir, paste0(pro_name, "_SquidpySpatialObj.h5ad")))
    }
    
    supported_format <- c("rds", "h5ad")
    if (!save_as %in% supported_format) {
        write("Warning: not saving to file or file format (specified in save_as) not supported, supported format: rds, h5ad", stderr())
    }

    return(object)

    ##bin data
    # out <- as.data.frame(dat)
    # out <- cbind(dat$y,dat$x,out)
    # colnames(out)[1:2] <- c(paste0(rowname_style,bs,".y"),paste0(rowname_style,bs,".x"))
    # fwrite(out,paste0(pro_name,"_bin",bs,"_information.txt"),col.names=T,row.names=F,sep="\t",quote=F)
    # out <- as.data.frame(cbind(paste0(rowname_style,unique(dat$bin_ID)),bin.coor$y,bin.coor$x))
    # colnames(out) <- c(paste0(rowname_style,bs),paste0(rowname_style,bs,".y"),paste0(rowname_style,bs,".x"))
    # rownames(out) <- out[,1]
    # fwrite(out,paste0(pro_name,"_bin",bs,"_position.txt"),col.names=T,row.names=F,sep="\t",quote=F)
    # bkpos <- out

    ##cell mask
    # out <- as.data.frame(dat)
    # out <- cbind(dat$y,dat$x,out)
    # colnames(out)[1:2] <- c(paste0(rowname_style,".y"),paste0(rowname_style,".x"))
    # fwrite(out,paste0(pro_name,"_cell","_information.txt"),col.names=T,row.names=F,sep="\t",quote=F)
    # out <- as.data.frame(cbind(paste0(rowname_style,unique(dat$bin_ID)),bin.coor$y,bin.coor$x))
    # colnames(out) <- c(paste0(rowname_style,bs),paste0(rowname_style,bs,".y"),paste0(rowname_style,bs,".x"))
    # rownames(out) <- out[,1]
    # fwrite(out,paste0(pro_name,"_bin",bs,"_position.txt"),col.names=T,row.names=F,sep="\t",quote=F)
    # bkpos <- out
}

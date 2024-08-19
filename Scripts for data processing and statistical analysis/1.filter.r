library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

## Create Seurat object
data_name = samples
object.list <- list()
for (i in seq(data_name)) {
	mat <- Read10X(data.dir = data_name[i], gene.column = 1)
	object.list[[i]] <- CreateSeuratObject(counts = mat, project = data_name[i], assay = assay)
}

## merge Seurat object
object <- merge(x = object.list[[1]], y = unlist(object.list[-1]),add.cell.ids = data_name)


# filter
object <- subset(object, subset = Percent_mito < 20)
object <- subset(object, subset = nFeature_RNA > 400)

# 过滤红细胞基因
rbc_genes = c("ENSGALG00000047152", "ENSGALG00000017347", "ENSGALG00000031597", "ENSGALG00000043234", "ENSGALG00000035309", "ENSGALG00000028273", "ENSGALG00000023740")
object <- PercentageFeatureSet(object, features = rbc_genes, col.name = "percent.hb")
object <- subset(object, percent.hb <= 0.01*100)

# 筛选表达 DAZL 的细胞
cell.use = Cells(obj)[obj[['RNA']]@data["ENSGALG00000011243", ] > 0]
object <- object[,cell.use]

save(object, file="obj.filter.Rda")


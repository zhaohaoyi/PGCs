library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

object <- Read10X("obj.filter.Rda")

## integration
object <- NormalizeData(object, normalization.method = "LogNormalize",scale.factor = 10000)
object <- FindVariableFeatures(object, selection.method = "vst",nfeatures = 2000)
object <- ScaleData(object = object, vars.to.regress = NULL,features = rownames(object))

object.list <- SplitObject(object, split.by = "orig.ident")
object.list <- SplitObject.Image(object.list)

object <- FindIntegrationAnchors(object.list,
								 dims = 1:50, normalization.method = "LogNormalize",
								 anchor.features = 3000, k.filter = 200)

object <- IntegrateData(anchorset = object,dims = 1:50, normalization.method = "LogNormalize")
object <- ScaleData(object, verbose = FALSE)
object <- RunPCA(object,assay = "integrated",npcs = 50,features = VariableFeatures(object), verbose = FALSE)
object <- RunUMAP(object, dims = 1:50, umap.method = "uwot")
object <- RunTSNE(object, dims = seq(50))

## Find clusters
object <- FindNeighbors(object, reduction = "pca", dims = seq(50),force.recalc = TRUE)
object <- FindClusters(object, resolution = 0.5, temp.file.location = getwd())
# set color of clusters
object@misc$color.cluster <- rainbow(nlevels(object@meta.data$seurat_clusters))
names(object@misc$color.cluster) <- levels(object@meta.data$seurat_clusters)

## Draw t-SNE UMAP plot
p1 <- DimPlot(object, reduction = 'tsne', group.by = "orig.ident", cols=object@misc$color.cluster, label = FALSE)
p2 <- DimPlot(object, reduction = 'tsne', group.by = "seurat_clusters", cols=object@misc$color.cluster, label = TRUE)
ggsave(p1+p2, file = "TSNE.pdf")
p1 <- DimPlot(object, reduction = 'umap', group.by = "orig.ident", cols=object@misc$color.cluster, label = FALSE)
p2 <- DimPlot(object, reduction = 'umap', group.by = "seurat_clusters", cols=object@misc$color.cluster, label = TRUE)
ggsave(p1+p2, file = "UMAP.pdf")

## Save data object 
DefaultAssay(object) <- "RNA"
save(object, file="obj.Rda")


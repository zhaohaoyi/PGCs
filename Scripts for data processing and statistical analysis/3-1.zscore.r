library(Seurat)
library(dplyr)
library(ggplot2)

obj = Load("obj.Rda")
genesets <- readLines("gene.list")

base_res = .GetMetaData(obj, c(Cluster = "seurat_clusters", Samples = "orig.ident", Groups = "Groups"))

## z-score
zscore_res = base_res
for (i in genesets) {
  zscore_res[, i] = as.vector(scale(Matrix::colMeans(obj[['RNA']]@data[gene_list[[i]], ])))
}
saveRDS(zscore_res, file = "zscore_results.Rds")
write.table(zscore_res, file = "zscore_results.xls", quote = FALSE, sep = "\t", row.names = FALSE)

# plot violin
for (j in genesets) {
	j2 = gsub(" |\\/", "_", j)
	df = df[, c("Cluster", j)]
	p = ggviolin(df, x = "", y = "Score", 
				 fill = "Cluster", color = "Cluster", 
				 palette = obj@misc$color.cluster, alpha = 1, width = 0.75) +
	    ggtitle(j) +
	ggsave(p, filename = paste0(j2, ".violin.pdf"), height = 6, width=8)
}

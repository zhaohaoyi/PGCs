library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

obj     <- Load("obj.Rda")
cc.gene <- "marker.yaml"
cc.gene <- .FindCCGene(obj, cc.gene)
database2table(cc.gene, "CellCycle.genes.xls", obj)

name <- "Cell.Cycle"
null <- "non-cycling"
ident.name = "Phase"
obj.cc <- AddModuleScore(obj = obj, features = cc.gene, name = name,
							ctrl = min(vapply(X = cc.gene, FUN = length, FUN.VALUE = numeric(length = 1))))
cc.columns <- grep(pattern = name, x = colnames(x = obj.cc[[]]), value = TRUE)
cc.scores <- obj.cc[[cc.columns]]
assign <- names(cc.gene)
assignments <- .CC.assignments(cc.scores, thres = 0, assign = assign, null = null)
cc.scores <- merge(x = cc.scores, y = data.frame(assignments), by = 0)
colnames(x = cc.scores) <- c("rownames", assign, ident.name)
rownames(x = cc.scores) <- cc.scores$rownames
cc.scores <- cc.scores[, c(assign, ident.name)]
obj[[colnames(x = cc.scores)]] <- cc.scores
obj@meta.data[[ident.name]] <- factor(obj@meta.data[[ident.name]], levels = c(assign, null))

vars <- c("G1", "S", "G2", "M")
cc_annot <- obj@meta.data %>%
	tibble::rownames_to_column("Cells") %>%
	mutate(CellCycle.Score = apply(.[,vars], 1, max)) %>%
	select(Cells, Samples = orig.ident, Groups = Groups, Clusters = seurat_clusters,
		!!!enquos(vars), CellCycle.Score, Phase)
write.table(cc_annot, file = "CellCycle.annot.xls"quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)


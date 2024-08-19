library(Seurat)
library(monocle)
library(cowplot)
library(igraph)
library(dplyr)

#1.Trajectory 
##obj created by Seurat
load("obj.Rda")
obj@meta.data$Samples  <- obj@meta.data$orig.ident
obj@meta.data$Clusters <- obj@meta.data$seurat_clusters

data <- as(as.matrix(obj@assays$RNA@counts), "sparseMatrix")
pd <- new("AnnotatedDataFrame", data = obj@meta.data)
fData <- data.frame(gene_short_name = row.names(data),row.names = row.names(data))
fd <- new("AnnotatedDataFrame", data = fData)
mnc_obj <- newCellDataSet(data, phenoData = pd, featureData = fd,lowerDetectionLimit = 0, expressionFamily = negbinomial.size())
mnc_obj <- setOrderingFilter(mnc_obj, ordering_genes = obj@assays[["RNA"]]@var.features)

mnc_obj <- reduceDimension(mnc_obj, method = 'DDRTree', verbose = F, scaling = T, max_components = 2, norm_method = 'none' , maxIter = 10, ncenter = NULL, tol = 0.001, sigma = 1, lambda = NULL, param.gamma = 10 )
mnc_obj <- orderCells(mnc_obj)
pData(mnc_obj) <- droplevels(pData(mnc_obj))

##draw
p1 <- plot_cell_trajectory(mnc_obj, color_by = "Samples",show_branch_points = F)
p2 <- p1 + facet_wrap(~Samples,nrow = 1)
p3 <- plot_cell_trajectory(mnc_obj, color_by = "Clusters",show_branch_points = F)
p4 <- p3 + facet_wrap(~Clusters,nrow = 1)
p5 <- plot_cell_trajectory(mnc_obj, color_by = "State",show_branch_points = F)
p6 <- p5 + facet_wrap(~State,nrow = 1)

##state diff
diff_state_res <- differentialGeneTest(mnc_obj, fullModelFormulaStr = "~State", cores = 4)
saveRDS( diff_state_res, file = "diff_state.rds" )
sig_gene_names <- rownames(subset( diff_state_res[order(diff_state_res$qval),], qval < 1e-7))
p1 <- plot_genes_jitter(mnc_obj[head(sig_gene_names, 10),], grouping = "State", color_by = "State", ncol = 5)
ggsave( "Diff.state.pdf", p1, width = 12, height = 5, limitsize = FALSE )


##pseudotime diff
diff_Pseudotime_res  <- differentialGeneTest(mnc_obj, fullModelFormulaStr = "~sm.ns(Pseudotime)", cores = cores)
saveRDS( diff_Pseudotime_res, file = "diff_Pseudotime.rds" )
sig_gene_names <- rownames(subset( diff_Pseudotime_res[order(diff_Pseudotime_res$qval),], qval < 1e-7))
p1 <- plot_genes_in_pseudotime(mnc_obj[head(sig_gene_names,10),], ncol = 5, color_by = "Samples")
p2 <- plot_pseudotime_heatmap(mnc_obj[sig_gene_names,], num_clusters = nlevels(pData(mnc_obj)$Clusters), cores = 4, show_rownames = T, return_heatmap = T)
ggsave("Diff.genes_in_pseudotime.pdf", p1, width = 12, height = 5, limitsize = FALSE )
ggsave("Diff.pseudotime_heatmap.pdf", p2, width = 10, height = 30, limitsize = FALSE )

##branch diff
BEAM_res <- BEAM(mnc_obj, branch_point = 1, cores = 4)
saveRDS( BEAM_res, file = "BEAM.1.rds")
sig_gene_names <- rownames(subset( BEAM_res[order(BEAM_res$qval),], qval < 1e-7))
p1 <- plot_genes_branched_pseudotime(mnc_obj[head(sig_gene_names,10),], branch_point = 1, color_by = "Samples", ncol = 5)
p2 <- plot_genes_branched_heatmap(mnc_obj[sig_gene_names,], branch_point = 1, num_clusters = 5, cores = 4, show_rownames = T, return_heatmap=T)
ggsave("Branch.1.genes_pseudotime.pdf", p1, width = 12, height = 5, limitsize = FALSE )
ggsave("Branch.1.genes_heatmap.pdf", p2, width = 10, height = 30, limitsize = FALSE )


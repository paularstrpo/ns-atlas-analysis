library(tidyverse)
library(Matrix)
library(Seurat)
library(SeuratObject)
setwd("D:/Dropbox/_projects/NormalSkinAtlas")

# library(future)
# future::plan("multisession", workers = 4, gc = TRUE)
# future::plan("sequential")
# options(future.globals.maxSize = 16e+09, future.rng.onMisuse = "ignore", future.seed=TRUE)


obj <- readRDS('data/scrna/seurat_objects/normal_skin.scRNA.merged.seurat_object.RDS')

cells_per_sample <- deframe(data.frame(table(obj@meta.data$donor_id)))
samples_to_keep <- names(cells_per_sample[cells_per_sample > 100])
obj <- subset(obj, donor_id %in% samples_to_keep)
ncol(obj)
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$donor_id)
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj, nfeatures = 2000)
obj <- ScaleData(obj, vars.to.regress = c("nCount_RNA", "pct.mito"))
obj <- RunPCA(obj, npcs = 30)
saveRDS(obj, file = 'data/scrna/seurat_objects/normal_skin.scRNA.unintegrated.seurat_object.RDS')


# obj <- readRDS('data/scrna/seurat_objects/normal_skin.scRNA.unintegrated.seurat_object.RDS')
gc()
obj <- IntegrateLayers(
  object = obj, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = TRUE
)

obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:20, graph.name = paste0('harmony.', c('nn', 'snn')))
# future::plan("sequential")
obj <- FindClusters(obj, resolution = c(0.2, 0.4, 0.5, 0.6, 0.8, 1, 1.2, 1.5, 2, 2.5), graph.name = 'harmony.snn')
obj <- RunUMAP(obj, dims = 1:20, reduction = "harmony", reduction.name = "harmony_umap",
               n.neighbors = 50,  n.epochs = 500, min.dist = 0.1, negative.sample.rate = 10)
obj[['RNA']] <- JoinLayers(obj[['RNA']])
saveRDS(obj, file = 'data/scrna/seurat_objects/normal_skin.scRNA.harmony.integrated.seurat_object.RDS')

cluster_prefix <- "harmony.snn"
res.cols <- colnames(obj@meta.data)[grepl(cluster_prefix, colnames(obj@meta.data))]

DimPlot(obj, reduction = "harmony_umap", group.by = res.cols, ncol = 4, label = TRUE) & NoLegend()

library(tidyverse)
library(presto)
cluster.marker.list <- list()
for (res in res.cols){
  print(paste0("Getting Markers for resolution ", res))

  obj@meta.data[, res] <- factor(obj@meta.data[, res], levels = as.character(0:max(as.numeric(as.character(na.omit(obj@meta.data[, res]))))))
  Idents(obj) <- res

  print("Cells per cluster:")
  print(table(obj@meta.data[, res]))

  markerdf <- wilcoxauc(obj, res, assay = 'data')
  markerdf$resolution <- res
  markerdf$diff.pct <- markerdf$pct_in - markerdf$pct_out
  markerdf$is.signif <- (markerdf$padj < 0.05)
  markerdf$is.positive <- markerdf$logFC > 0

  print("Markers per cluster:")
  print(table(markerdf$is.signif & markerdf$is.positive, markerdf$group))

  cluster.marker.list[[res]] <- markerdf
  markerdf
  write.csv(markerdf, file = paste0('data/scrna/cluster_markers/normal_skin.scRNA.harmony.cluster_markers.', res, '.csv'))
}

saveRDS(cluster.marker.list, file =  'data/scrna/cluster_markers/normal_skin.scRNA.harmony.harmony.cluster_marker_list.RDS')

qs::qsave(obj, file = 'data/scrna/seurat_objects/normal_skin.scRNA.harmony.integrated.seurat_object.QS')
saveRDS(obj, file = 'data/scrna/seurat_objects/normal_skin.scRNA.harmony.integrated.seurat_object.RDS')


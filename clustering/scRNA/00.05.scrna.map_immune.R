library(Seurat)
library(SeuratObject)
library(Matrix)
library(tidyverse)
library(harmony)
library(presto)
library(tidyverse)
library(BPCells)
root_path <- "D:/Dropbox/_projects/NormalSkinAtlas"
setwd(root_path)

process_reclustering <- function(obj, idents_to_subset,
                                 k.param = 50,
                                 ndims = 20,
                                 resolutions = c(0.2, 0.4, 0.5, 0.6, 0.8, 1, 1.5, 2, 2.5)
) {
  ## process the data
  reclust.obj <- NormalizeData(obj)
  reclust.obj <- FindVariableFeatures(reclust.obj, nfeatures = 2000)

  reclust.obj<- ScaleData(reclust.obj, features = VariableFeatures(reclust.obj), vars.to.regress = c("nCount_RNA", "pct.mito"))
  reclust.obj<- RunPCA(reclust.obj, features = VariableFeatures(reclust.obj), reduction.name = "reclust.pca")

  ## skip unintegrated clustering, go straight to harmony
  reclust.obj <- IntegrateLayers(reclust.obj, method = HarmonyIntegration, orig.reduction = 'reclust.pca', new.reduction = 'reclust.harmony')


  reclust.obj <- FindNeighbors(reclust.obj, reduction = "reclust.harmony", dims = 1:ndims, k.param = k.param)
  reclust.obj <- RunUMAP(reclust.obj, reduction = "reclust.harmony", dims = 1:ndims, reduction.name = 'reclust.umap',
                         n.epochs = 500, min.dist = 0.1, negative.sample.rate = 10, n.neighbors = 50)


  for (res in resolutions){
    res_col <- paste0("reclust.louvain_", res)
    print(res_col)

    reclust.obj <- FindClusters(reclust.obj, resolution = res, cluster.name = res_col)
    lvls <- as.character(0:max(as.numeric(as.character(na.omit(reclust.obj@meta.data[, res_col])))))
    reclust.obj@meta.data[, res_col] <- factor(as.character(reclust.obj@meta.data[, res_col]), levels = lvls)

    print(table(reclust.obj@meta.data[, res_col]))
  }

  return(reclust.obj)
}

## Immune
scrna <- qs::qread("data/scrna/seurat_objects/normal_skin.scRNA.harmony.integrated.seurat_object.annotated.filtered.QS")
scrna <- subset(scrna, cell_category == 'Immune')
gc()
scrna[['RNA']] <- JoinLayers(scrna[['RNA']])

samples <- table(scrna$donor_id)
samples <- names(samples[samples>30])
Idents(scrna) <- "donor_id"
scrna <- subset(scrna, idents=samples)

scrna[['RNA']] <- split(scrna[['RNA']], f = scrna$donor_id)
scrna <- process_reclustering(scrna)
qs::qsave(scrna, file = "data/scrna/reclustering/Immune.reclustered.seurat_object.QS")

merfish <- qs::qread("data/merfish/reclustering/Immune.reclustered.seurat_object.annotated.QS")
merfish <- subset(merfish, cell_type.detailed != "Doublet")
merfish.meta <- merfish@meta.data

overlap.genes <- intersect(rownames(merfish), rownames(scrna))
scrna[['RNA']] <- JoinLayers(scrna[['RNA']])

transfer.anchors <- FindTransferAnchors(reference = merfish, query = scrna, dims = 1:20,  features = overlap.genes)
saveRDS(transfer.anchors, file = "data/scrna/reclustering/Immune.merfish_cca_predictions.transfer_anchors.RDS")

predictions <- TransferData(anchorset = transfer.anchors, refdata = merfish.meta[colnames(merfish), "cell_type.detailed"], dims = 1:20)
saveRDS(predictions, file = "data/scrna/reclustering/Immune.merfish_cca_predictions.scrna.cell_type.detailed.RDS")

scrna <- AddMetaData(scrna, metadata = predictions$predicted.id, col.name = "cell_type.detailed.reclustered.cca.predicted_id")
scrna <- AddMetaData(scrna, metadata = predictions$prediction.score.max, col.name = "cell_type.detailed.reclustered.cca.prediction_score")
scrna <- AddMetaData(scrna, metadata = predictions)
qs::qsave(scrna, file = "data/scrna/reclustering/Immune.reclustered.seurat_object.QS")

DimPlot(scrna, reduction = 'reclust.umap',
        group.by = c("cell_type.detailed.reclustered.cca.predicted_id", "reclust.louvain_0.8", "reclust.louvain_1.5"),
        label=TRUE, alpha = 0.25, ncol=3) & NoLegend() & NoAxes() & coord_fixed()


library(presto)
library(tidyverse)
res.cols <- paste0("reclust.louvain_", c(0.2, 0.4, 0.5, 0.6, 0.8, 1, 1.5, 2, 2.5))
cluster.marker.list <- list()
for (res in res.cols){
  print(res)
  Idents(scrna) <- res
  markerdf <- presto::wilcoxauc(scrna, res, assay = 'data', seurat_assay = 'RNA')
  markerdf$resolution <- res
  markerdf$diff.pct <- markerdf$pct_in - markerdf$pct_out
  markerdf$is.signif <- (markerdf$padj < 0.05)
  markerdf$is.positive <- markerdf$logFC > 0
  print(table(markerdf$is.signif & markerdf$is.positive, markerdf$group))

  cluster.marker.list[[res]] <- markerdf

  markerdf <- markerdf %>% dplyr::filter(is.positive & is.signif)
  write.csv(markerdf, file = paste0('data/scrna/cluster_markers/Immune/Immune.cluster_markers.', res, '.csv'))
}
saveRDS(cluster.marker.list, file =  'data/scrna/cluster_markers/Immune/Immune.reclustering.cluster_marker_list.RDS')



scrna <- qs::qread("data/scrna/reclustering/Immune.reclustered.seurat_object.QS")
## now add the harmonized annotations
chosen_res = 'reclust.louvain_1'
celltype_annots <- readxl::read_xlsx("data/scrna/scrna.cell_type.annotations.xlsx", sheet = "Immune")%>% data.frame()
celltype_annots[, chosen_res] <- as.character(celltype_annots[, chosen_res])

meta <- merge(scrna@meta.data[, c("cell_barcode", chosen_res)], celltype_annots, by = chosen_res, all.x = TRUE, sort=FALSE)
rownames(meta) <- meta$cell_barcode
scrna <- AddMetaData(scrna, metadata = meta, col.name = 'cell_type.reclustered')
scrna <- UpdateSeuratObject(scrna)

scrna@meta.data$cell_type.reclustered <- factor(scrna@meta.data$cell_type.reclustered)
scrna@meta.data$cell_type.detailed <- scrna@meta.data$cell_type.reclustered


## save final objects
qs::qsave(scrna, file = "data/scrna/reclustering/Immune.reclustered.seurat_object.annotated.QS")
meta <- scrna@meta.data[, c("cell_barcode", "cell_type.detailed")]
saveRDS(meta, file =  "data/scrna/reclustering/Immune.reclustered.metadata.RDS")



# scrna <- qs::qread("data/scrna/reclustering/Immune.reclustered.seurat_object.annotated.QS")
#
# lym <- subset(scrna, cell_type.reclustered %in% c("T Lym", "B Lym"))
# lym <- JoinLayers(lym)
# lym[['RNA']] <- JoinLayers(lym[['RNA']])
#
# samples <- table(lym$donor_id)
# samples <- names(samples[samples>30])
#
# Idents(lym) <- "donor_id"
# lym <- subset(lym, idents=samples)
#
# lym[['RNA']] <- split(lym[['RNA']], f = lym$donor_id)
# lym <- process_reclustering(lym)
# lym <- JoinLayers(lym)
# lym[['RNA']] <- JoinLayers(lym[['RNA']])
# qs::qsave(lym, file = "data/scrna/reclustering/Immune.Lym.reclustered.seurat_object.QS")
#
# res.cols <- paste0("reclust.louvain_", c(0.2, 0.4, 0.5, 0.6, 0.8, 1, 1.5, 2, 2.5))
# library(tidyverse)
# library(presto)
# cluster.marker.list <- list()
# for (res in res.cols){
#   print(res)
#   Idents(lym) <- res
#   markerdf <- presto::wilcoxauc(lym, res, assay = 'data', seurat_assay = 'RNA')
#   markerdf$resolution <- res
#   markerdf$diff.pct <- markerdf$pct_in - markerdf$pct_out
#   markerdf$is.signif <- (markerdf$padj < 0.05)
#   markerdf$is.positive <- markerdf$logFC > 0
#   print(table(markerdf$is.signif & markerdf$is.positive, markerdf$group))
#
#   cluster.marker.list[[res]] <- markerdf
#
#   markerdf <- markerdf %>% dplyr::filter(is.positive & is.signif)
#   write.csv(markerdf, file = paste0('data/scrna/cluster_markers/Immune/Immune.Lym.cluster_markers.', res, '.csv'))
# }
#
# saveRDS(cluster.marker.list, file =  'data/scrna/cluster_markers/Immune/Immune.Lym.reclustering.cluster_marker_list.RDS')
#
# # lym <- qs::qread("data/scrna/reclustering/Immune.Lym.reclustered.seurat_object.QS")
#
# myl <- subset(scrna, cell_type.reclustered %in% c("Mac", "Mono", "DC", "LC"))
# myl <- JoinLayers(myl)
# myl[['RNA']] <- JoinLayers(myl[['RNA']])
#
# samples <- table(myl$donor_id)
# samples <- names(samples[samples>30])
#
# Idents(myl) <- "donor_id"
# myl <- subset(myl, idents=samples)
#
# myl[['RNA']] <- split(myl[['RNA']], f = myl$donor_id)
# myl <- process_reclustering(myl)
# myl <- JoinLayers(myl)
# myl[['RNA']] <- JoinLayers(myl[["RNA"]])
# qs::qsave(myl, file = "data/scrna/reclustering/Immune.Myeloid.reclustered.seurat_object.QS")
#
# res.cols <- paste0("reclust.louvain_", c(0.2, 0.4, 0.5, 0.6, 0.8, 1, 1.5, 2, 2.5))
#
# library(tidyverse)
# library(presto)
# cluster.marker.list <- list()
# for (res in res.cols){
#   print(res)
#   Idents(myl) <- res
#   markerdf <- presto::wilcoxauc(myl, res, assay = 'data', seurat_assay = 'RNA')
#   markerdf$resolution <- res
#   markerdf$diff.pct <- markerdf$pct_in - markerdf$pct_out
#   markerdf$is.signif <- (markerdf$padj < 0.05)
#   markerdf$is.positive <- markerdf$logFC > 0
#   print(table(markerdf$is.signif & markerdf$is.positive, markerdf$group))
#
#   cluster.marker.list[[res]] <- markerdf
#
#   markerdf <- markerdf %>% dplyr::filter(is.positive & is.signif)
#   write.csv(markerdf, file = paste0('data/scrna/cluster_markers/Immune/Immune.Myeloid.cluster_markers.', res, '.csv'))
# }
# saveRDS(cluster.marker.list, file =  'data/scrna/cluster_markers/Immune/Immune.Myeloid.reclustering.cluster_marker_list.RDS')



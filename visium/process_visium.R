library(tidyverse)
library(Seurat)
library(Matrix)
setwd("D:/Dropbox/_projects/NormalSkinAtlas")

study_ids <- c("ThraneJID2023", "JiCell2020",  "BergenstrahleNatBiotech2022", "GanierPNAS2024", "MitamuraAllergy2023", "MaNatComms2024", 'YuImmunity2024')
seurat_list <- lapply(study_ids, function(x) qs::qread(paste0("data/visium/seurat_objects/", x, ".visium_data.raw.seurat_object.QS")))
names(seurat_list) <- study_ids

obj <- merge(seurat_list[[1]], seurat_list[2:length(seurat_list)])
obj@meta.data$cell_barcode <- rownames(obj@meta.data)
obj <- JoinLayers(obj)
rm(seurat_list); gc()
obj <- split(obj, f = obj$sample_id)

## export spatial coords as a dimmreduc
spatial_coords <- lapply(obj@images, function(x){
  coords <- x$centroids@coords
  rownames(coords) <- x$centroids@cells
  colnames(coords) <- c("coords_1", "coords_2")
  return(coords)
}) %>% do.call(what=rbind)

obj[['coords']] <- CreateDimReducObject(embeddings = spatial_coords, assay="Spatial")
obj$coords_1 <- spatial_coords[,1]
obj$coords_2 <- spatial_coords[,2]
obj@meta.data$spatial_1 <- spatial_coords[,1]
obj@meta.data$spatial_2 <- spatial_coords[,2]
rm(spatial_coords)

obj <- UpdateSeuratObject(obj)
qs::qsave(obj, file = 'data/visium/seurat_objects/pan-skin_merged.visium_data.raw.seurat_object.QS')

VlnPlot(obj, features = c("nCount_Spatial", "nFeature_Spatial"), pt.size = 0, group.by = "study_id", log = TRUE) &
  geom_boxplot(aes(fill=NULL), width = 0.2) &
  #geom_hline(yintercept = 10, color = 'black', linetype = "dashed") +
  #labs(x = NULL, y = "nCount_Spatial", title = NULL) +
  theme_bw() & NoLegend() & coord_flip()


obj <- subset(obj, nCount_Spatial > 10)
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj, vars.to.regress = "nFeature_Spatial")
obj <- RunPCA(obj)

obj <- IntegrateLayers(
  object = obj, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = TRUE
)

obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:20, graph.name = paste0('harmony.', c('nn', 'snn')))
obj <- FindClusters(obj, resolution = c(0.2, 0.4, 0.5, 0.6, 0.8, 1, 1.2, 1.5, 2, 2.5), graph.name = 'harmony.snn')
obj <- RunUMAP(obj, dims = 1:20, reduction = "harmony", reduction.name = "harmony_umap",
               n.neighbors = 50,  n.epochs = 500, min.dist = 0.1, negative.sample.rate = 10)
obj <- JoinLayers(obj)



qs::qsave(obj, file = 'data/visium/seurat_objects/pan-skin_merged.visium_data.harmony_integrated.seurat_object.QS')
saveRDS(obj, file =  'data/visium/seurat_objects/pan-skin_merged.visium_data.harmony_integrated.seurat_object.RDS')
p1 <- DimPlot(obj, reduction = "harmony_umap", group.by = c( "sample_id","donor_id"), label=FALSE, raster.dpi = c(600,600), alpha=0.5, ncol=1) & NoLegend() & NoAxes() & coord_fixed()
p2 <- DimPlot(obj, reduction = "harmony_umap", group.by = c( "study_id", "tissue_type"), label=TRUE, raster.dpi = c(600,600), alpha=0.5, ncol=1) & NoAxes() & coord_fixed()
p3 <- DimPlot(obj, reduction = "harmony_umap", group.by = c("harmony.snn_res.0.6", "harmony.snn_res.0.8"), label=TRUE, raster.dpi = c(600,600), alpha=0.5, ncol=1) & NoLegend() & NoAxes() & coord_fixed()

featureqc_panel <- FeaturePlot(obj, reduction = 'harmony_umap', features = c("nCount_Spatial", "nFeature_Spatial", "pct.mito"), raster.dpi = c(600,600), alpha=0.5, ncol=3) & NoLegend() & NoAxes() & coord_fixed()


((p3| p1 | p2) / featureqc_panel) + patchwork::plot_layout(heights = c(2, 1))

## get spatial cluster markers!
data_path <- "data/visium/"
run_prefix <- "pan-skin_merged.visium_data"
resolutions <- c(0.2, 0.4, 0.5, 0.6, 0.8, 1, 1.2, 1.5, 2, 2.5)
res.cols <- paste0("harmony.snn_res.", resolutions)
library(presto)
library(tidyverse)
cluster.marker.list <- list()
for (res in res.cols){
  print(res)
  Idents(obj) <- res
  markerdf <- wilcoxauc(obj, res, assay = 'data', seurat_assay = 'Spatial')
  markerdf$resolution <- res
  markerdf$diff.pct <- markerdf$pct_in - markerdf$pct_out
  markerdf$is.signif <- (markerdf$padj < 0.05)
  markerdf$is.positive <- markerdf$logFC > 0
  print(table(markerdf$is.signif & markerdf$is.positive, markerdf$group))

  cluster.marker.list[[res]] <- markerdf

  markerdf <- markerdf %>% dplyr::filter(is.positive & is.signif)
  write.csv(markerdf, file = paste0(data_path, 'cluster_markers/', run_prefix, '.cluster_markers.', res, '.csv'))
}

saveRDS(cluster.marker.list, file =  paste0(data_path, 'cluster_markers/', run_prefix, '.full.harmony.cluster_marker_list.RDS'))
saveRDS(cluster.marker.list, file =  paste0(data_path, 'cluster_markers/', run_prefix, '.harmony.cluster_marker_list.RDS'))

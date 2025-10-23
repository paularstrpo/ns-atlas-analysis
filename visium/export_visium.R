library(tidyverse)
library(Seurat)
library(Matrix)
library(rhdf5)
library(reticulate)
library(anndata)
setwd("D:/Dropbox/_projects/NormalSkinAtlas")


obj <- qs::qread('data/visium/seurat_objects/pan-skin_merged.visium_data.harmony_integrated.seurat_object.QS')

meta <- obj@meta.data
visium_meta <- unique(obj@meta.data[, c("sample_id", "donor_id", "study_id", "tissue_type", "anatomic_site", "donor_age", "donor_sex")])


write.csv(obj@meta.data, file = "data/visium/metadata/pan-skin_merged.visium_data.harmony_integrated.spot_metadata.csv")
saveRDS(obj@meta.data, file = "data/visium/metadata/pan-skin_merged.visium_data.harmony_integrated.spot_metadata.RDS")

library(SingleCellExperiment)
sce <- SingleCellExperiment(
  assays = list(
    counts = obj[['Spatial']]$counts,
    logcounts = obj[['Spatial']]$data
  ),
  colData = obj@meta.data,
  reducedDims = SimpleList(
    PCA = Embeddings(obj, reduction = 'pca'),
    HARMONY = Embeddings(obj, reduction = 'harmony'),
    UMAP = Embeddings(obj, reduction = 'harmony_umap'),
    SPATIAL = Embeddings(obj, reduction = 'coords')
  )
)

saveRDS(sce, file = "data/visium/seurat_objects/pan-skin_merged.visium_data.harmony_integrated.sce_object.RDS")
qs::qsave(sce, file = "data/visium/seurat_objects/pan-skin_merged.visium_data.harmony_integrated.sce_object.QS")
rm(sce); gc()

adata <- AnnData(
  X = t(obj[['Spatial']]$data),
  obs = obj@meta.data,
  layers = list(
    counts = t(obj[['Spatial']]$counts),
    log1p = t(obj[['Spatial']]$data)
    ),
  obsm = list(
    pca = Embeddings(obj, reduction = 'pca'),
    harmony = Embeddings(obj, reduction = 'harmony'),
    umap = Embeddings(obj, reduction = 'harmony_umap'),
    spatial = Embeddings(obj, reduction = 'coords')
    )
  )
write_h5ad(adata, filename = "data/visium/seurat_objects/pan-skin_merged.visium_data.harmony_integrated.anndata.h5ad")

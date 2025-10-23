library(tidyverse)
library(Seurat)
library(Matrix)
library(rhdf5)
library(reticulate)
library(anndata)
setwd("D:/Dropbox/_projects/NormalSkinAtlas")


obj <- qs::qread("data/scrna/seurat_objects/normal_skin.scRNA.harmony.integrated.seurat_object.reclustered.annotated.filtered.QS")

adata <- AnnData(
  X = t(obj[['RNA']]$data),
  obs = obj@meta.data,
  layers = list(
    counts = t(obj[['RNA']]$counts),
    log1p = t(obj[['RNA']]$data)
  ),
  obsm = list(
    pca = Embeddings(obj, reduction = 'pca'),
    harmony = Embeddings(obj, reduction = 'harmony'),
    umap = Embeddings(obj, reduction = 'harmony_umap')
  )
)
write_h5ad(adata, filename = "data/scrna/seurat_objects/normal_skin.scRNA.harmony.integrated.reclustered.annotated.filtered.anndata.h5ad")
rm(adata, obj); gc()

## Stroma
obj <- qs::qread("data/scrna/reclustering/Stroma.reclustered.seurat_object.annotated.QS")
obj <- subset(obj, cell_type.detailed != "Doublet")
obj <- JoinLayers(obj)
obj[['RNA']] <- JoinLayers(obj[["RNA"]])
adata <- AnnData(
  X = t(obj[['RNA']]$data),
  obs = obj@meta.data,
  layers = list(
    counts = t(obj[['RNA']]$counts),
    log1p = t(obj[['RNA']]$data)
  ),
  obsm = list(
    pca = Embeddings(obj, reduction = 'pca'),
    harmony = Embeddings(obj, reduction = 'harmony'),
    umap = Embeddings(obj, reduction = 'harmony_umap')
  )
)
write_h5ad(adata, filename = "data/scrna/reclustering/Stroma.reclustered.seurat_object.annotated.anndata.h5ad")

## Immune
obj <- qs::qread("data/scrna/reclustering/Immune.reclustered.seurat_object.annotated.QS")
obj <- subset(obj, cell_type.detailed != "Doublet")
obj <- JoinLayers(obj)
obj[['RNA']] <- JoinLayers(obj[["RNA"]])
adata <- AnnData(
  X = t(obj[['RNA']]$data),
  obs = obj@meta.data,
  layers = list(
    counts = t(obj[['RNA']]$counts),
    log1p = t(obj[['RNA']]$data)
  ),
  obsm = list(
    pca = Embeddings(obj, reduction = 'pca'),
    harmony = Embeddings(obj, reduction = 'harmony'),
    umap = Embeddings(obj, reduction = 'harmony_umap')
  )
)
write_h5ad(adata, filename = "data/scrna/reclustering/Immune.reclustered.seurat_object.annotated.anndata.h5ad")

## Epithelia
obj <- qs::qread("data/scrna/reclustering/Epithelia.reclustered.seurat_object.annotated.QS")
obj <- subset(obj, cell_type.detailed != "Doublet")
obj <- JoinLayers(obj)
obj[['RNA']] <- JoinLayers(obj[["RNA"]])
adata <- AnnData(
  X = t(obj[['RNA']]$data),
  obs = obj@meta.data,
  layers = list(
    counts = t(obj[['RNA']]$counts),
    log1p = t(obj[['RNA']]$data)
  ),
  obsm = list(
    pca = Embeddings(obj, reduction = 'pca'),
    harmony = Embeddings(obj, reduction = 'harmony'),
    umap = Embeddings(obj, reduction = 'harmony_umap')
  )
)
write_h5ad(adata, filename = "data/scrna/reclustering/Epithelia.reclustered.seurat_object.annotated.anndata.h5ad")


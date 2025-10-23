root_path <- "D:/Dropbox/_projects/NormalSkinAtlas"
setwd(root_path)

library(tidyverse)
library(Matrix)
library(Seurat)
library(SeuratObject)
library(reticulate)
library(rhdf5)
library(anndata)

h5ad <- "data/merfish/BAYSOR/seurat_objects/ns-atlas.merfish_baysor.scanvi_integrated.cellcharter.anndata.annotated.h5ad"
adata <- anndata::read_h5ad(h5ad)
counts <- adata$layers['counts']
counts <- Matrix(t(counts), sparse = TRUE)
gc()
log1p <- adata$layers['log1p']
log1p <- Matrix(t(log1p), sparse=TRUE)
gc()

# Add observation metadata
metadata <- h5read(h5ad, name = 'obs')
metadata <- lapply(metadata, function(x){
  # print(x)

  if (class(x) == 'list') {
    res <- factor(x[['codes']], labels = setNames(x[['categories']], 0:(length(x[['categories']]) - 1)))
    return(res)
  }
  else {
    return(x)
  }
}) %>%
  bind_cols() %>%
  data.frame()
rownames(metadata) <- metadata$cell_barcode

barcodes <- intersect(rownames(metadata), colnames(counts))
metadata <- metadata[barcodes,]
counts <- counts[, barcodes]
log1p <- log1p[, barcodes]
gc()

# Create the Seurat object
obj <- CreateSeuratObject(counts = counts, meta.data = metadata)
DefaultAssay(obj) <- 'RNA'
obj[['RNA']]$data <- log1p

for (rep in names(adata$obsm)){
  mat <- as.matrix(adata$obsm[[rep]])
  rownames(mat) <- adata$obs_names
  colnames(mat) <- paste0(rep, "_", 1:ncol(mat))
  mat <- mat[barcodes, ]
  obj[[rep]] <- CreateDimReducObject(mat, key = paste0(gsub("_", "", rep), "_"), assay = 'RNA')
  rm(mat); gc()
}
obj <- UpdateSeuratObject(obj)
rm(counts, metadata, log1p, adata, barcodes, h5ad); gc()

qs::qsave(obj, file = "data/merfish/BAYSOR/seurat_objects/ns-atlas.merfish_baysor.scanvi_integrated.cellcharter.annotated.seurat_object.QS")

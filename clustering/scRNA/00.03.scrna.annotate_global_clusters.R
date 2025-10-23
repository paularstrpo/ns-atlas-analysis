library(tidyverse)
library(Matrix)
library(Seurat)
library(SeuratObject)
obj <- qs::qread("data/scrna/seurat_objects/normal_skin.scRNA.harmony.integrated.seurat_object.QS")
obj@meta.data$cell_barcode <- rownames(obj@meta.data)

chosen_res = 'harmony.snn_res.0.2'
celltype_annots <- readxl::read_xlsx("data/scrna/scrna.cell_type.annotations.xlsx", sheet = chosen_res) %>% data.frame()
celltype_annots[, chosen_res] <- as.character(celltype_annots[, chosen_res])
meta <- merge(obj@meta.data[, c("cell_barcode", chosen_res)], celltype_annots, by = chosen_res, all.x = TRUE, sort=FALSE)
rownames(meta) <- meta$cell_barcode
obj <- AddMetaData(obj, metadata = meta, col.name = 'cell_type.broad')
obj <- AddMetaData(obj, metadata = meta, col.name = 'cell_category')
obj <- UpdateSeuratObject(obj)
obj@meta.data$cell_type <- factor(obj@meta.data$cell_type.broad)
obj@meta.data$cell_type.broad.res_0.2 <- obj@meta.data$cell_type
obj@meta.data$cell_category <- factor(obj@meta.data$cell_category)
qs::qsave(obj, file = "data/scrna/seurat_objects/normal_skin.scRNA.harmony.integrated.seurat_object.annotated.QS")
saveRDS(obj, file = "data/scrna/seurat_objects/normal_skin.scRNA.harmony.integrated.seurat_object.annotated.RDS")

spatial_theme <- theme(
  rect = element_blank(),
  line = element_blank(),
  text = element_text(color = 'black', size = 10),
  title = element_text(face = 'bold', size = 12, color = 'black', hjust = 0.5),
  strip.background = element_rect(fill = 'black', color = 'black'),
  strip.text = element_text(color = 'white', face = 'bold', size = 12),
  plot.title = element_text(size = 14, color = 'black', face = 'bold', hjust = 0.5),
  axis.text = element_blank()
)

plt <- DimPlot(obj, reduction = 'harmony_umap', group.by = "cell_type.broad", alpha = 0.25, label = TRUE) +
  spatial_theme + NoLegend() + NoAxes() + coord_fixed()



png("figures/scrna/annotated/scrna.global_integration.cell_type.broad.umap.with_labels.png", width=14, height=8, units = 'in', res=300)
print(plt)
dev.off()

svglite::svglite("figures/scrna/annotated/scrna.global_integration.cell_type.broad.umap.with_labels.svg", width=14, height=8)
print(plt)
dev.off()


plt <- DimPlot(obj, reduction = 'harmony_umap', group.by = "cell_type.broad", alpha = 0.25, label = FALSE, raster = FALSE) +
  spatial_theme + NoLegend() + NoAxes() + coord_fixed()

png("figures/scrna/annotated/scrna.global_integration.cell_type.broad.umap.no_labels.png", width=14, height=8, units = 'in', res=300)
print(plt)
dev.off()


obj.filt <- subset(obj, subset = cell_type.broad != "Doublet")
qs::qsave(obj.filt, file = "data/scrna/seurat_objects/normal_skin.scRNA.harmony.integrated.seurat_object.annotated.filtered.QS")

plt <- DimPlot(obj.filt, reduction = 'harmony_umap', group.by = "cell_type.broad", alpha = 0.25, label = TRUE) +
  spatial_theme + NoLegend() + NoAxes() + coord_fixed()


png("figures/scrna/annotated/scrna.global_integration.cell_type.broad.umap.with_labels.no_doublets.png", width=14, height=8, units = 'in', res=300)
print(plt)
dev.off()

svglite::svglite("figures/scrna/annotated/scrna.global_integration.cell_type.broad.umap.with_labels.no_doublets.svg", width=14, height=8)
print(plt)
dev.off()


plt <- DimPlot(obj.filt, reduction = 'harmony_umap', group.by = "cell_type.broad", alpha = 0.25, label = FALSE, raster = FALSE) +
  spatial_theme + NoLegend() + NoAxes() + coord_fixed()

png("figures/scrna/annotated/scrna.global_integration.cell_type.broad.umap.no_labels.no_doublets.png", width=14, height=8, units = 'in', res=300)
print(plt)
dev.off()
obj <- qs::qread("data/scrna/seurat_objects/normal_skin.scRNA.harmony.integrated.seurat_object.reclustered.annotated.filtered.QS")
chosen_res = 'harmony.snn_res.1'
celltype_annots <- readxl::read_xlsx("data/scrna/scrna.cell_type.annotations.xlsx", sheet = chosen_res) %>% data.frame()
celltype_annots[, chosen_res] <- as.character(celltype_annots[, chosen_res])
meta <- merge(obj@meta.data[, c("cell_barcode", chosen_res)], celltype_annots, by = chosen_res, all.x = TRUE, sort=FALSE)
rownames(meta) <- meta$cell_barcode
obj <- AddMetaData(obj, metadata = meta, col.name = 'cell_type.broad')
obj <- UpdateSeuratObject(obj)
obj@meta.data$cell_type.broad.res_1 <- obj@meta.data$cell_type.broad

DimPlot(obj, reduction = 'harmony_umap', group.by = "cell_type.broad", alpha = 0.25, label = TRUE)

qs::qsave(obj, file = "data/scrna/seurat_objects/normal_skin.scRNA.harmony.integrated.seurat_object.reclustered.annotated.filtered.QS")

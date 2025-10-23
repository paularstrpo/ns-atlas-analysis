library(Seurat)
library(SeuratObject)
library(Matrix)
library(tidyverse)
library(harmony)
library(presto)
library(tidyverse)
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


## 1. Epithelia
scrna <- qs::qread("data/scrna/seurat_objects/normal_skin.scRNA.harmony.integrated.seurat_object.annotated.filtered.QS")
scrna <- subset(scrna, cell_category == 'Epithelia')
gc()

scrna[['RNA']] <- JoinLayers(scrna[['RNA']])
samples <- table(scrna$donor_id) > 30
table(samples)
names(samples[!samples])
samples_to_keep <- names(samples[samples])
Idents(scrna) <- "donor_id"
scrna <- subset(scrna, idents = samples_to_keep)
gc()
scrna[['RNA']] <- split(scrna[['RNA']], f = scrna$donor_id)
scrna <- process_reclustering(scrna)
qs::qsave(scrna, file = "data/scrna/reclustering/Epithelia.reclustered.seurat_object.QS")
gc()


merfish <- qs::qread("data/merfish/reclustering/Epithelia.reclustered.seurat_object.annotated.QS")
merfish <- subset(merfish, cell_type.detailed != "Plasma" & cell_type.detailed != "LC" & cell_type.detailed != "Doublet")
merfish.meta <- merfish@meta.data

merfish.umap <- DimPlot(merfish, reduction  = "reclust.umap", group.by = "cell_type.detailed", label = TRUE) + NoLegend() + NoAxes() + labs(title = "MERFISH")
scrna.umap <- DimPlot(scrna, reduction  = "reclust.umap", group.by = "cell_type.broad", label = TRUE) + NoLegend() + NoAxes() + labs(title = "scRNA")

merfish.umap | scrna.umap

overlap.genes <- intersect(rownames(merfish), rownames(scrna))
scrna[['RNA']] <- JoinLayers(scrna[['RNA']])

transfer.anchors <- FindTransferAnchors(reference = merfish, query = scrna, dims = 1:20,  features = overlap.genes)
saveRDS(transfer.anchors, file = "data/scrna/reclustering/Epithelia.merfish_cca_predictions.transfer_anchors.RDS")

predictions <- TransferData(anchorset = transfer.anchors, refdata = merfish.meta[colnames(merfish), "cell_type.detailed"], dims = 1:20)
saveRDS(predictions, file = "data/scrna/reclustering/Epithelia.merfish_cca_predictions.scrna.cell_type.detailed.RDS")

scrna <- AddMetaData(scrna, metadata = predictions$predicted.id, col.name = "cell_type.detailed.reclustered.cca.predicted_id")
scrna <- AddMetaData(scrna, metadata = predictions$prediction.score.max, col.name = "cell_type.detailed.reclustered.cca.prediction_score")
scrna <- AddMetaData(scrna, metadata = predictions)
qs::qsave(scrna, file = "data/scrna/reclustering/Epithelia.reclustered.seurat_object.QS")

DimPlot(scrna, reduction = 'reclust.umap',
        group.by = c("cell_type.detailed.reclustered.cca.predicted_id", "reclust.louvain_0.8"),
        label=TRUE, alpha = 0.25) & NoLegend() & NoAxes() & coord_fixed()


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
  write.csv(markerdf, file = paste0('data/scrna/cluster_markers/Epithelia/Epithelia.cluster_markers.', res, '.csv'))
}

saveRDS(cluster.marker.list, file =  'data/scrna/cluster_markers/Epithelia/Epithelia.reclustering.cluster_marker_list.RDS')


scrna <- qs::qread("data/scrna/reclustering/Epithelia.reclustered.seurat_object.QS")

## now add the harmonized annotations
chosen_res = 'reclust.louvain_0.6'
celltype_annots <- readxl::read_xlsx("data/scrna/scrna.cell_type.annotations.xlsx", sheet = "Epithelia")%>% data.frame()
celltype_annots[, chosen_res] <- as.character(celltype_annots[, chosen_res])
meta <- merge(scrna@meta.data[, c("cell_barcode", chosen_res)], celltype_annots, by = chosen_res, all.x = TRUE, sort=FALSE)
rownames(meta) <- meta$cell_barcode
scrna <- AddMetaData(scrna, metadata = meta, col.name = 'cell_type.reclustered')
scrna <- UpdateSeuratObject(scrna)
scrna@meta.data$cell_type.reclustered <- factor(scrna@meta.data$cell_type.reclustered)
DimPlot(scrna, reduction = 'reclust.umap',
        group.by = "cell_type.reclustered",
        label=TRUE, alpha = 0.25) + NoAxes()
qs::qsave(scrna, file = "data/scrna/reclustering/Epithelia.reclustered.seurat_object.annotated.QS")


#
#
# ## recluster the Spn and Grn KCs
#
# scrna <- qs::qread("data/scrna/reclustering/Epithelia.reclustered.seurat_object.annotated.QS")
# Idents(scrna) <- "cell_type.reclustered"
# scrna[['RNA']] <- JoinLayers(scrna[['RNA']])
#
# diff.ife <- subset(scrna, idents = "Diff KC")
# samples <- table(diff.ife$donor_id)
# samples <- names(samples[samples>30])
# Idents(diff.ife) <- "donor_id"
# diff.ife <- subset(diff.ife, idents = samples)
# diff.ife[['RNA']] <- split(diff.ife[['RNA']], f = diff.ife$donor_id)
# diff.ife <- process_reclustering(diff.ife)
# qs::qsave(diff.ife, file = "data/scrna/reclustering/Epithelia.DiffIFE.reclustered.seurat_object.QS")
# gc()
#
# diff.ife <- qs::qread("data/scrna/reclustering/Epithelia.DiffIFE.reclustered.seurat_object.QS")
#
# DimPlot(diff.ife, reduction = 'reclust.umap',label = TRUE,
#         group.by = c("cell_type.detailed.reclustered.cca.predicted_id","reclust.louvain_0.4",
#                      "cell_type.reclustered", "anatomic_site", "sample_barcode", "study_id")) &
#   NoLegend() & NoAxes()
#
# ump <- DimPlot(diff.ife, reduction = 'reclust.umap', group.by = "reclust.louvain_0.4", label = TRUE, alpha = 0.25) & NoLegend() & NoAxes()
#
# bar <- ggplot(diff.ife@meta.data) +
#   aes(x = reclust.louvain_0.4, fill = study_id) +
#   geom_bar(position='fill', color='black') +
#   scale_y_continuous(labels = scales::percent) +
#   theme_bw()
# bar2 <- ggplot(diff.ife@meta.data) +
#   aes(x = reclust.louvain_0.4, fill = anatomic_site) +
#   geom_bar(position='fill', color='black') +
#   scale_y_continuous(labels = scales::percent) +
#   theme_bw()
#
# ump | bar | bar2
#
#
# cells <- diff.ife$anatomic_site == 'sole'
# cells <- names(cells[cells])
# sole <- DimPlot(diff.ife, reduction = 'reclust.umap', cells.highlight = cells) + NoAxes() + NoLegend() + labs(title = "Sole")
# # sole
#
# cells <- diff.ife$anatomic_site == 'palm'
# cells <- names(cells[cells])
# palm <- DimPlot(diff.ife, reduction = 'reclust.umap', cells.highlight = cells) + NoAxes() + NoLegend() + labs(title = "Palm")
# # palm
#
# palm | sole
#
# diff.markers <- c("KRT1", "KRT10", "SOX9", "GJB2", "GJB6", "S100A8", "S100A9", "FOSL1", "GRHL3")
# diff.featurepanel <- FeaturePlot(diff.ife, reduction = 'reclust.umap', features = diff.markers, order = TRUE) & NoAxes() & NoLegend()
# diff.featurepanel
#
# diff.celltypes <- c("Cheng 2018", "Spn KC I", "Spn KC II", "Grn KC", "Spn KC I", "Spn KC I", "Ji 2020")
# names(diff.celltypes) <- c(0:6)
# diff.celltypes
# diff.ife$cell_type.reclustered <- factor(diff.ife$reclust.louvain_0.4, levels = names(diff.celltypes), labels = diff.celltypes)
#
# diff.anno <- DimPlot(diff.ife, reduction = "reclust.umap", group.by = "cell_type.reclustered", label=TRUE, alpha = 0.25) + NoAxes() + NoLegend() + labs(title = "Diff IFE")
#
# qs::qsave(diff.ife, file = "data/scrna/reclustering/Epithelia.DiffIFE.reclustered.seurat_object.annotated.QS")
#
#
# Idents(scrna) <- "cell_type.reclustered"
# bas.ife <- subset(scrna, idents = c("Bas KC"))
# samples <- table(bas.ife$donor_id)
# samples <- names(samples[samples>30])
# Idents(bas.ife) <- "donor_id"
# bas.ife <- subset(bas.ife, idents = samples)
# bas.ife[['RNA']] <- split(bas.ife[['RNA']], f = bas.ife$donor_id)
# bas.ife <- process_reclustering(bas.ife)
# qs::qsave(bas.ife, file = "data/scrna/reclustering/Epithelia.BasIFE.reclustered.seurat_object.QS")
#
# gc()
# # rm(scrna); gc()
# bas.ife <- qs::qread("data/scrna/reclustering/Epithelia.BasIFE.reclustered.seurat_object.QS")
#
# pan <- DimPlot(bas.ife, reduction = 'reclust.umap', label = TRUE, alpha=0.25,
#                group.by =  c("cell_type.detailed.reclustered.cca.predicted_id","reclust.louvain_0.4",
#                              "cell_type.reclustered", "anatomic_site", "sample_barcode", "study_id")) &
#   NoLegend() & NoAxes() & coord_fixed()
#
# pan
#
#
# cells <- bas.ife$anatomic_site == 'sole'
# cells <- names(cells[cells])
# sole <- DimPlot(bas.ife, reduction = 'reclust.umap', cells.highlight = cells) + NoAxes() + NoLegend() + labs(title = "Sole")
# # sole
#
# cells <- bas.ife$anatomic_site == 'palm'
# cells <- names(cells[cells])
# palm <- DimPlot(bas.ife, reduction = 'reclust.umap', cells.highlight = cells) + NoAxes() + NoLegend() + labs(title = "Palm")
# # palm
#
# palm | sole
#
# library(patchwork)
#
# ump <- DimPlot(bas.ife, reduction = 'reclust.umap', group.by = "reclust.louvain_1", label = TRUE, alpha=0.25) & NoLegend() & NoAxes()
#
# bar <- ggplot(bas.ife@meta.data) +
#   aes(x = reclust.louvain_0.8, fill = study_id) +
#   geom_bar(position='fill', color='black') +
#   scale_y_continuous(labels = scales::percent) +
#   theme_bw()
# bar2 <- ggplot(bas.ife@meta.data) +
#   aes(x = reclust.louvain_0.8, fill = anatomic_site) +
#   geom_bar(position='fill', color='black') +
#   scale_y_continuous(labels = scales::percent) +
#   theme_bw()
#
# (ump | bar | bar2) + plot_annotation(title = "Bas IFE")
# bas.markers <- c("KRT5", "KRT14",  "COL17A1", "ITGA6", "ITGB4", "SOX9", "S100A8", "S100A9", "GJB2", "GJB6")
# bas.featurepanel <- FeaturePlot(bas.ife, reduction = 'reclust.umap', features = bas.markers, order=TRUE) & NoAxes() & NoLegend()
# bas.featurepanel
#
#
# bas.celltypes <- c("Cheng 2018", "Bas KC I", "Bas KC", "Bas KC", "Cheng 2018", "Bas KC", "Bas KC", "Cheng 2018","Cheng 2018")
# names(bas.celltypes) <- c(0:8)
# bas.celltypes
# bas.ife$cell_type.reclustered <- factor(bas.ife$reclust.louvain_0.6, levels = names(bas.celltypes), labels = bas.celltypes)
#
# bas.anno <- DimPlot(bas.ife, reduction = "reclust.umap", group.by = "cell_type.reclustered", label=TRUE, alpha=0.25)+ NoAxes() + NoLegend() + labs(title = "Bas IFE")
# bas.anno
# bas.anno | diff.anno
#
# qs::qsave(bas.ife, file = "data/scrna/reclustering/Epithelia.BasIFE.reclustered.seurat_object.annotated.QS")
#
# meta <- rbind(bas.ife@meta.data, diff.ife@meta.data)
# scrna@meta.data$cell_type.reclustered.old <- as.character(scrna@meta.data$cell_type.reclustered)
# scrna@meta.data$cell_type.reclustered.round2 <- as.character(scrna@meta.data$cell_type.reclustered)
# scrna@meta.data[rownames(meta),]$cell_type.reclustered.round2 <- as.character(meta$cell_type.reclustered)
# scrna@meta.data$cell_type.reclustered.round2[scrna@meta.data$cell_type.reclustered.round2 == "Bas KC"] <- "Bas KC I"
# scrna@meta.data$cell_type.reclustered.round2[scrna@meta.data$cell_type.reclustered.round2 == "Diff KC"] <- "Spn KC I"
# scrna@meta.data$cell_type.reclustered <- as.character(scrna@meta.data$cell_type.reclustered.round2)
# DimPlot(scrna, reduction = "reclust.umap", group.by = 'cell_type.reclustered.round2', label=TRUE)

scrna@meta.data$cell_type.detailed <- scrna@meta.data$cell_type.reclustered
scrna@meta.data$cell_type.detailed[scrna@meta.data$cell_type.detailed %in% c("Cheng 2018", "Ji 2020", "Doublet")] <- "Doublet"
DimPlot(subset(scrna, cell_type.detailed != "Doublet"), reduction = "reclust.umap", group.by = 'cell_type.detailed', label=TRUE)

qs::qsave(scrna, file = "data/scrna/reclustering/Epithelia.reclustered.seurat_object.annotated.QS")
meta <- scrna@meta.data[, c("cell_barcode", "cell_type.reclustered", "cell_type.detailed", "cell_type.detailed.reclustered.cca.predicted_id")]
saveRDS(meta, file =  "data/scrna/reclustering/Epithelia.reclustered.metadata.RDS")


scrna <- qs::qread("data/scrna/reclustering/Epithelia.reclustered.seurat_object.annotated.QS")
scrna <- subset(scrna, cell_type.detailed != "Doublet")

# libraries
library(tidyverse)
library(Matrix)
library(Seurat)
library(SeuratObject)
library(circlize)
library(RColorBrewer)
library(ggrepel)
library(ggpubr)
library(patchwork)
library(ragg)
library(ggrastr)
library(paletteer)

## ggplot themes
theme <- theme(
  legend.position = 'right',
  text = element_text(color = 'black', size = 10),
  title = element_text(face = 'bold', size = 12, color = 'black', hjust = 0.5),
  strip.background = element_rect(fill = 'black', color = 'black'),
  strip.text = element_text(color = 'white', face = 'bold', size = 12),
  plot.title = element_text(size = 14, color = 'black', face = 'bold', hjust = 0.5),
  panel.background = element_rect(color = 'black', fill = NA),
  plot.background = element_blank(),
  legend.background = element_blank(),
  legend.key = element_blank()
)


global.colors <- c("lightcyan", "royalblue2")
global.cell_type.colors <- readRDS("data/reference/cell_type.lvl1.color_palette.final.RDS")
detailed.cell_type.colors <- readRDS("data/reference/cell_type.detailed.color_palette.final.RDS")

detailed.cell_types <- names(detailed.cell_type.colors)
broad.cell_types <- names(global.cell_type.colors)

global.dotplot.marker_genes <- list(
  "Cyc" = c("MKI67", "TOP2A"),
  "Basal KC" = c("COL17A1", "COL7A1"),
  "Spn KC" = c("RHOV", "KLF5"),
  "Grn KC" = c("DSC1", "FLG"),
  "Seb" = c("KRT79", "TKT"),
  "Ecc Gland" = c("SFRP1", "DNER"),
  "Bas HFE" = c("GJB2", "GJB6"),
  "Diff HFE" = c("S100A8", "SOX9"),
  "DP/DS" = c("COCH", "POSTN", "COL12A1"),
  "Fib" = c("COL1A1","PDGFRA", "PDGFRB"),
  "Peri" = c("ACTA2", "NOTCH3","RGS5"),
  "EC" = c("VWF", "ACKR1"),
  "LEC" = c("CCL21", "KDR", "LYVE1"),
  "SM" = c("DES", "MCAM", "MYH11"),
  "Adipo" = c("ADIPOQ", "G0S2",  "LPL"),
  "Schwann" = c("SOX10", "MPZ", "DKK3"),
  "T Cell"=c("IL32", "CD3E", "CD4", "CD8A"),
  "DC" = c("CD83", "ITGAX", "CLEC10A"),
  "Mac" = c("MRC1", "C1QA", "C1QB"),
  "Mast" = c("CPA3", "IL1RL1", "KIT")
)

epi.markers.final <-  c("TOP2A", "CDC20",
                        "COL17A1", "DST",
                        "IGFBP3","ITGA3", "PTN",
                        "GATA3", "KLF5",
                        "NOTCH3", "S100A8", "DSC1",
                        "KLK7","IL37","PRDM1", "FLG",
                        "SEMA3C", "INHBB", "SFRP1",  "DNER",
                        "THBS1", "KRT79", "TKT",
                        "MLANA", "MITF", "SOX9","IGFL2", "COL18A1",
                        "FGFR1", "KRT28", "KRT74")

epi.filt <- subset(scrna, cell_type.detailed != "Doublet")

cell_types <- c("Cyc KC", "Bas KC I", "Bas KC II", "Spn KC I", "Spn KC II", "Grn KC I", "Grn KC II", "Ecc Duct", "Ecc Gland", "Seb", "Melano", "Bas Inf", "Diff Inf", "Bulge/ORS Basal", "ORS Suprabasal", "IRS/HS")
epi.filt$cell_type.detailed <- factor(epi.filt$cell_type.detailed, levels = cell_types)

epi.dotplot <- DotPlot(epi.filt, group.by = "cell_type.detailed", features = epi.markers.final,
                       cols = global.colors, col.min = 0, dot.min = 0.05) +
  labs(title = NULL, x = NULL, y = NULL) + theme + coord_fixed()+
  scale_size_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100)) +
  guides(color = guide_colorbar(title = "Average\nExpression", frame.colour = 'black'),
         size = guide_legend(title = "Percent\nExpressed")) +
  theme(panel.grid.minor = element_line(color = "black", size = 0.5),
        panel.margin.x=unit(0.25, "lines"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12, color = 'black'),
        axis.text.y = element_text(hjust = 1, vjust = 0.6, size = 12, color = 'black')
  )
epi.dotplot

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

## Stroma

scrna <- qs::qread("data/scrna/seurat_objects/normal_skin.scRNA.harmony.integrated.seurat_object.annotated.filtered.QS")
scrna <- subset(scrna, cell_category %in% c("Stroma", "Fibroblast"))
gc()
scrna[['RNA']] <- JoinLayers(scrna[['RNA']])

samples <- table(scrna$donor_id)
samples <- names(samples[samples>30])
Idents(scrna) <- "donor_id"
scrna <- subset(scrna, idents=samples)

scrna[['RNA']] <- split(scrna[['RNA']], f = scrna$donor_id)
scrna <- process_reclustering(scrna)
qs::qsave(scrna, file = "data/scrna/reclustering/Stroma.reclustered.seurat_object.QS")

merfish <- qs::qread("data/merfish/reclustering/Stroma.reclustered.seurat_object.annotated.QS")
merfish <- subset(merfish, cell_type.detailed != "Doublet" & cell_type.detailed != "Adipo")
merfish.meta <- merfish@meta.data

overlap.genes <- intersect(rownames(merfish), rownames(scrna))
scrna[['RNA']] <- JoinLayers(scrna[['RNA']])

transfer.anchors <- FindTransferAnchors(reference = merfish, query = scrna, dims = 1:20,  features = overlap.genes)
saveRDS(transfer.anchors, file = "data/scrna/reclustering/Stroma.merfish_cca_predictions.transfer_anchors.RDS")

predictions <- TransferData(anchorset = transfer.anchors, refdata = merfish.meta[colnames(merfish), "cell_type.detailed"], dims = 1:20)
saveRDS(predictions, file = "data/scrna/reclustering/Stroma.merfish_cca_predictions.scrna.cell_type.detailed.RDS")

scrna <- AddMetaData(scrna, metadata = predictions$predicted.id, col.name = "cell_type.detailed.reclustered.cca.predicted_id")
scrna <- AddMetaData(scrna, metadata = predictions$prediction.score.max, col.name = "cell_type.detailed.reclustered.cca.prediction_score")
scrna <- AddMetaData(scrna, metadata = predictions)
qs::qsave(scrna, file = "data/scrna/reclustering/Stroma.reclustered.seurat_object.QS")

DimPlot(scrna, reduction = 'reclust.umap',
        group.by = c("cell_type.detailed.reclustered.cca.predicted_id","cell_type.broad", "reclust.louvain_0.8", "reclust.louvain_1.5", "reclust.louvain_2.5"),
        label=TRUE, alpha = 0.25, ncol=5) & NoLegend() & NoAxes() & coord_fixed()

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
  write.csv(markerdf, file = paste0('data/scrna/cluster_markers/Stroma/Stroma.cluster_markers.', res, '.csv'))
}

saveRDS(cluster.marker.list, file =  'data/scrna/cluster_markers/Stroma/Stroma.reclustering.cluster_marker_list.RDS')



scrna <- qs::qread("data/scrna/reclustering/Stroma.reclustered.seurat_object.QS")


## now add the harmonized annotations
chosen_res = 'reclust.louvain_0.8'
celltype_annots <- readxl::read_xlsx("data/scrna/scrna.cell_type.annotations.xlsx", sheet = "Stroma")%>% data.frame()
celltype_annots[, chosen_res] <- as.character(celltype_annots[, chosen_res])

meta <- merge(scrna@meta.data[, c("cell_barcode", chosen_res)], celltype_annots, by = chosen_res, all.x = TRUE, sort=FALSE)
rownames(meta) <- meta$cell_barcode
scrna$cell_type.reclustered <- NULL
scrna <- AddMetaData(scrna, metadata = meta, col.name = 'cell_type.reclustered')
scrna <- UpdateSeuratObject(scrna)

scrna@meta.data$cell_type.reclustered <- factor(scrna@meta.data$cell_type.reclustered)
qs::qsave(scrna, file = "data/scrna/reclustering/Stroma.reclustered.seurat_object.annotated.QS")


scrna <- qs::qread("data/scrna/reclustering/Stroma.reclustered.seurat_object.annotated.QS")

# get the fibs out and recluster
# include the schwann cells cause the signature for them is messy
Idents(scrna) <- "cell_type.reclustered"
fib <- subset(scrna, idents = c("DP/DS", "Papil Fib", "Perivasc Fib", "Retic Fib", "Schwann"))
fib[['RNA']] <- split(fib[['RNA']], f = fib$donor_id)
fib <- process_reclustering(fib)
qs::qsave(fib, file = "data/scrna/reclustering/Stroma.Fib.reclustered.seurat_object.QS")


pan <- DimPlot(fib, reduction = 'reclust.umap',
        group.by = c("cell_type.detailed.reclustered.cca.predicted_id", "cell_type.reclustered", "reclust.louvain_0.8"), label = TRUE) & NoAxes() & NoLegend()

pan

cells <- fib$anatomic_site == 'sole'
cells <- names(cells[cells])
sole <- DimPlot(fib, reduction = 'reclust.umap', cells.highlight = cells) + NoAxes() + NoLegend() + labs(title = "Sole")
sole

cells <- fib$anatomic_site == 'palm'
cells <- names(cells[cells])
palm <- DimPlot(fib, reduction = 'reclust.umap', cells.highlight = cells) + NoAxes() + NoLegend() + labs(title = "Palm")
palm

palm | sole

cells <- fib$anatomic_site == 'hip'
cells <- names(cells[cells])
hip <- DimPlot(fib, reduction = 'reclust.umap', cells.highlight = cells) + NoAxes() + NoLegend() + labs(title = "Hip")
hip

library(patchwork)
# (pan | sol) + plot_layout(widths = c(2, 1))

fib.markers <- c("MT1X", "APOE","SFRP1","SFRP2","WIF1","APCDD1","CXCL12","CCL19","CCDC80","PI16","PLA2G2A","MFAP5","MMP2", "COMP", "COL18A1", "ANGPTL7", "PRG4")

retic.fib.iii.markers <- c("THBS4", "ANGPTL7", "PRG4", "COMP", "CLU", "LUM", "ASPN", "BGN", "CRTAC1")
retic.fib.iv.markers <- c("CLDN1", "ITGA6", "TM4SF1", "NR2F2", "C2orf40", "FOXD1", "TNNC1", "TAGLN", "PTGDS", "APOD")

FeaturePlot(fib, reduction = 'reclust.umap', features = retic.fib.iii.markers, order=TRUE, raster=TRUE) & NoAxes()
FeaturePlot(fib, reduction = 'reclust.umap', features = retic.fib.iv.markers, order=TRUE, raster=TRUE) & NoAxes()

FeaturePlot(fib, reduction = 'reclust.umap', features = fib.markers, order=TRUE, raster=TRUE) & NoAxes()


library(ggpubr)
DotPlot(fib, features = unique(c(fib.markers, retic.fib.iii.markers, retic.fib.iV.markers)), group.by = "reclust.louvain_0.8") + coord_flip()

# VlnPlot(fib, features =  unique(c(fib.markers, retic.fib.iii.markers, retic.fib.iV.markers)), group.by = "reclust.louvain_0.8")


fib <- JoinLayers(fib)
cluster.marker.list <- list()
for (res in res.cols){
  print(res)
  Idents(fib) <- res
  markerdf <- presto::wilcoxauc(fib, res, assay = 'data', seurat_assay = 'RNA')
  markerdf$resolution <- res
  markerdf$diff.pct <- markerdf$pct_in - markerdf$pct_out
  markerdf$is.signif <- (markerdf$padj < 0.05)
  markerdf$is.positive <- markerdf$logFC > 0
  print(table(markerdf$is.signif & markerdf$is.positive, markerdf$group))

  cluster.marker.list[[res]] <- markerdf

  markerdf <- markerdf %>% dplyr::filter(is.positive & is.signif)
  write.csv(markerdf, file = paste0('data/scrna/cluster_markers/Stroma/Stroma.Fib.cluster_markers.', res, '.csv'))
}

saveRDS(cluster.marker.list, file =  'data/scrna/cluster_markers/Stroma/Stroma.Fib.reclustering.cluster_marker_list.RDS')



fib <- qs::qread("data/scrna/reclustering/Stroma.Fib.reclustered.seurat_object.QS")

## now add the harmonized annotations
chosen_res = 'reclust.louvain_0.8'
celltype_annots <- readxl::read_xlsx("data/scrna/scrna.cell_type.annotations.xlsx", sheet = "Stroma.Fib")%>% data.frame()
celltype_annots[, chosen_res] <- as.character(celltype_annots[, chosen_res])

meta <- merge(fib@meta.data[, c("cell_barcode", chosen_res)], celltype_annots, by = chosen_res, all.x = TRUE, sort=FALSE)
rownames(meta) <- meta$cell_barcode
fib$cell_type.reclustered <- NULL
fib <- AddMetaData(fib, metadata = meta, col.name = 'cell_type.reclustered')
fib <- UpdateSeuratObject(fib)

fib@meta.data$cell_type.reclustered <- factor(fib@meta.data$cell_type.reclustered)
qs::qsave(fib, file = "data/scrna/reclustering/Stroma.Fib.reclustered.seurat_object.annotated.QS")



scrna <- qs::qread("data/scrna/reclustering/Stroma.reclustered.seurat_object.annotated.QS")

meta <- fib@meta.data
scrna@meta.data$cell_type.reclustered.old <- as.character(scrna@meta.data$cell_type.reclustered)
scrna@meta.data$cell_type.reclustered.round2 <- as.character(scrna@meta.data$cell_type.reclustered)
scrna@meta.data[rownames(meta),]$cell_type.reclustered.round2 <- as.character(meta$cell_type.reclustered)

scrna@meta.data$cell_type.reclustered <- as.character(scrna@meta.data$cell_type.reclustered.round2)
# DimPlot(scrna, reduction = "reclust.umap", group.by = 'cell_type.reclustered.round2', label=TRUE)


scrna@meta.data$cell_type.detailed <- scrna@meta.data$cell_type.reclustered
DimPlot(scrna, reduction = "reclust.umap", group.by = 'cell_type.detailed', label=TRUE)

## save final objects
qs::qsave(scrna, file = "data/scrna/reclustering/Stroma.reclustered.seurat_object.annotated.QS")
meta <- scrna@meta.data[, c("cell_barcode", "cell_type.detailed")]
saveRDS(meta, file =  "data/scrna/reclustering/Stroma.reclustered.metadata.RDS")


scrna <- qs::qread("data/scrna/reclustering/Stroma.reclustered.seurat_object.annotated.QS")

fib.markers <- c("APOE","SFRP1","SFRP2","WIF1","APCDD1","CXCL12","CCL19","CCDC80","PI16","PLA2G2A","MFAP5","MMP2", "COMP", "COL18A1", "ANGPTL7")

str.markers.final <- c("COL1A1", "LRRC15", "POSTN",
                       "COMP", "WIF1", "COL18A1","THBS2", "COL6A1", "PDGFRA", "APOE",
                       "CXCL12", "CCL19", "C3", "C1R", "MMP2","PRG4", "PDGFRB",
                       "CCDC80", "PLA2G2A","SFRP1", "ANGPTL7", "PTGDS",  "KLF5", "ITGA6",
                       "DES", "ITGA8", "MYH11", "ACTA2", "JAG1", "NOTCH3",
                       "RGS5", "VWF", "ACKR1","KDR", "CCL21", "TFPI",
                       "MPZ", "SOX2", "SOX10")
str.filt <- subset(scrna, cell_type.detailed != "Doublet")
str.filt@meta.data$cell_type.detailed <- factor(str.filt@meta.data$cell_type.detailed,
                                                levels = c("DP/DS", "Papil Fib", "Perivasc Fib I", "Perivasc Fib II",
                                                           "Retic Fib I", "Retic Fib II", "Retic Fib III", "Retic Fib IV",
                                                           "SM", "Peri", "HEC", "VEC", "LEC", "Schwann"))
global.colors <- c("lightcyan", "royalblue2")

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

str.dotplot <- DotPlot(str.filt, group.by = "cell_type.detailed", features = str.markers.final,
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
svglite::svglite("figures/scrna/annotated/stroma.reclustered.dot_plot.annotated.svg", height = 6, width = 14)
str.dotplot & coord_fixed()
dev.off()

svglite::svglite("figures/scrna/annotated/stroma.reclustered.umap.annotated.svg", height = 8, width = 10)
DimPlot(str.filt, group.by = "cell_type.detailed", reduction = "reclust.umap", label=TRUE, alpha = 0.25) & NoAxes() & coord_fixed()
dev.off()

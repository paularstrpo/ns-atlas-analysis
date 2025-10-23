library(tidyverse)
library(Matrix)
library(Seurat)
library(SeuratObject)
library(circlize)
library(RColorBrewer)
library(ggrepel)
library(ggpubr)
library(patchwork)
library(ggrastr)
library(svglite)

global_colors <- c("lightcyan", "royalblue2")
cluster_colors <- unique(c(pals::polychrome(n=36), pals::glasbey(n=32)))
names(cluster_colors) = as.character(c(1:length(cluster_colors))-1)
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


cluster_prefix <- "scvi.leiden_"
run_prefix <- "ns-atlas.merfish"
data_path <- "data/merfish/BAYSOR/"
resolutions <-  c(0.2, 0.4, 0.5, 0.6, 0.8, 1, 1.5, 2, 2.5)
setwd("D:/Dropbox/_projects/NormalSkinAtlas")

obj <- qs::qread("data/merfish/BAYSOR/seurat_objects/ns-atlas.merfish_baysor.scvi_intergated.seurat_object.QS")

cluster.marker.list <- readRDS("data/merfish/BAYSOR/cluster_markers/ns-atlas.merfish_baysor.scvi.cluster_marker_list.seurat.RDS")

CellTypeMarkerList <- list(
  Adipo = c("ADIPOQ", "LRP1"),
  APM = c("DES", "MCAM", "MYH11"),
  Schwann = c("NRXN1", "MPZ"),

  LEC = c("LYVE1", "ANGPT2", "CCL21"),
  EC = c("VWF", "ACKR1", "FLT1", "KDR"),
  Peri = c("ACTA2","RGS5", "TCIM",  "RGS16"),

  Fib = c("PDGFRA","PDGFRB","APOE", "COL1A1", "DPP4", "SFRP2", "FMO1", "MMP2", "CCL19"),

  DP = c("COCH", "CRABP1", "FIBIN", "PLXDC1"),
  DS = c("ACAN", "COL11A1", "DPEP1", "LRRC15", "MYL4", "TNN", "ITGA8"),

  HF_HS = c("KRT83", "LGR5", "EPCAM", "LMCD1"),
  HS = c("KRT32", "KRT35"),

  IRF_IRS_IST = c("PTN", "KRT28", "KRT74", "FST", "PTHLH"),
  Bulge = c("CD34", "DIO2", "RCAN1", "WIF1"),
  ORS.B = c("ANGPTL7", "THBS1"),

  Seb = c("KRT79", "CTSB", "APOC1", "ASS1", "SCARB1", "TKT"),
  Ecc.Duct = c("KRT77", "SEMA3C", "GJB6", "TBX3","FGF7"),
  Ecc.Gland = c("SFRP1", "DNER", "SOX9", "TSPAN8", "PLA2G2A","CLDN3", "CLDN10", "ABI3BP"),

  Cyc = c("MKI67", "PCNA", "POLA1", "TOP2A"),

  Bas.KC = c("COL17A1", "COL7A1", "DST", "ITGA6", "ITGB1", "ITGB4"),
  Spn.KC = c("MT1X", "MT1E", "S100A8", "DSC1", "GSTA3", "HES1", "NOTCH3", "RHOV", "SERPINB4"),
  Grn.KC = c("FLG","KLK7", "IL1RN"),
  Melano = c("FRZB", "MITF", "MLANA", "SOX10", "TYR"),

  T.Cell=c("CD3E","CD3D","CD3G"), # control
  NK.T=c("CD8A", 'NKG7', "NCR1", "XCL1", "KLRK1", "GZMB","GZMA","GZMK"),
  APC=c("CD4", "HLA-DRB5"),
  T.Reg=c("FOXP3", "CTLA4"),
  T.Help = c("IL9", "RORC", "STAT2", "STAT6", "TBX21"),
  B.Cell=c("MS4A1", "CD19", "JCHAIN", "CD22", "CD79A", "CD27", "CD40"),

  DC = c("CCR7","CD80", "BATF3", "FLT3", "IRF4", "IL3RA", "TCL1A", "IRF7",
         "CCL17", "CCL4", "CD200", "IDO1", "CD1C", "CLEC10A", "CLEC9A", "IRF8"),
  LC = c("CD207", "CD1A", "S100B"),

  Mono = c("CD14", "C5AR1"),
  Mac = c("C1QA", "C1QB", "CD68", "CSF1R", "MRC1")
)


### clustering plots
for (r in resolutions){
  res <- paste0("scvi.leiden_", r)

  ## first make the plot itself
  plt <- DimPlot(obj, reduction = 'umap', group.by = res, alpha=0.25, label = TRUE, raster.dpi=c(1200,1200)) + NoAxes() + NoLegend()

  svglite(paste0("figures/merfish/BAYSOR/cluster_umap/", run_prefix, '.cluster_umap.', res, ".svg"), width=7, height=7)
  print(plt)
  dev.off()

  png(paste0("figures/merfish/BAYSOR/cluster_umap/", run_prefix, '.cluster_umap.', res, ".png"), width=7, height=7, units = 'in', res=300)
  print(plt)
  dev.off()

 #  rm(plt)

  ### Proportion barplots
  obj@meta.data$seurat_clusters <- obj@meta.data[, res]

  # donor id
  plt <- ggplot(obj@meta.data) +
    aes(x = seurat_clusters, fill = donor_id) +
    geom_bar(position='fill', color='black') +
    theme_classic() + theme + labs(title = res, x = "cluster") +
    scale_y_continuous(labels = scales::percent)

  png(paste0("figures/merfish/BAYSOR/proportion_barplots/", run_prefix, '.proportion_barplot.donor_id.', res, ".png"), width=14, height=8, units = 'in', res=300)
  print(plt)
  dev.off()

  svglite::svglite(paste0("figures/merfish/BAYSOR/proportion_barplots/", run_prefix, '.proportion_barplot.donor_id.', res, ".svg"), width=14, height=8)
  print(plt)
  dev.off()


  # sample barcode
  plt <- ggplot(obj@meta.data) +
    aes(x = seurat_clusters, fill = sample_barcode) +
    geom_bar(position='fill', color='black') +
    theme_classic() + theme +
    labs(title = res, x = "cluster") +
    scale_y_continuous(labels = scales::percent)

  png(paste0("figures/merfish/BAYSOR/proportion_barplots/", run_prefix, '.proportion_barplot.sample_barcode.', res, ".png"), width=14, height=8, units = 'in', res=300)
  print(plt)
  dev.off()

  svglite::svglite(paste0("figures/merfish/BAYSOR/proportion_barplots/", run_prefix, '.proportion_barplot.sample_barcode.', res, ".svg"), width=14, height=8)
  print(plt)
  dev.off()


  # imaging run batch
  plt <- ggplot(obj@meta.data) +
    aes(x = seurat_clusters, fill = batch) +
    geom_bar(position='fill', color='black') +
    theme_classic() + theme +
    labs(title = res, x = "cluster") +
    scale_y_continuous(labels = scales::percent)


  png(paste0("figures/merfish/BAYSOR/proportion_barplots/", run_prefix, '.proportion_barplot.imaging_batch.', res, ".png"), width=14, height=8, units = 'in', res=300)
  print(plt)
  dev.off()

  svglite::svglite(paste0("figures/merfish/BAYSOR/proportion_barplots/", run_prefix, '.proportion_barplot.imaging_batch.', res, ".png"), width=14, height=8)
  print(plt)
  dev.off()

  # gene panel
  plt <- ggplot(obj@meta.data) +
    aes(x = seurat_clusters, fill = gene_panel) +
    geom_bar(position='fill', color='black') +
    theme_classic() + theme + labs(title = res, x = "cluster", fill = "Gene\nPanel") +
    scale_y_continuous(labels = scales::percent)


  png(paste0("figures/merfish/BAYSOR/proportion_barplots/", run_prefix, '.proportion_barplot.gene_panel.', res, ".png"), width=14, height=8, units = 'in', res=300)
  print(plt)
  dev.off()

  svglite(paste0("figures/merfish/BAYSOR/proportion_barplots/", run_prefix, '.proportion_barplot.gene_panel.', res, ".svg"), width=14, height=8)
  print(plt)
  dev.off()

  # sex
  plt <- ggplot(obj@meta.data) +
    aes(x = seurat_clusters, fill = donor_sex) +
    geom_bar(position='fill', color='black') +
    theme_classic() + theme +
    labs(title = res, x = "cluster") +
    scale_fill_manual(values = c("Male"="dodgerblue", "Female"="maroon2")) +
    scale_y_continuous(labels = scales::percent)


  png(paste0("figures/merfish/BAYSOR/proportion_barplots/", run_prefix, '.proportion_barplot.donor_sex.', res, ".png"), width=14, height=8, units = 'in', res=300)
  print(plt)
  dev.off()


  svglite::svglite(paste0("figures/merfish/BAYSOR/proportion_barplots/", run_prefix, '.proportion_barplot.donor_sex.', res, ".svg"), width=14, height=8)
  print(plt)
  dev.off()

  # anatomic site
  plt <- ggplot(obj@meta.data) +
    aes(x = seurat_clusters, fill = anatomic_site) +
    geom_bar(position='fill', color='black') +
    theme_classic() + theme +
    labs(title = res, x = "cluster") +
    scale_y_continuous(labels = scales::percent)


  png(paste0("figures/merfish/BAYSOR/proportion_barplots/", run_prefix, '.proportion_barplot.anatomic_site.', res, ".png"), width=14, height=8, units = 'in', res=300)
  print(plt)
  dev.off()


  svglite::svglite(paste0("figures/merfish/BAYSOR/proportion_barplots/", run_prefix, '.proportion_barplot.anatomic_site.', res, ".svg"), width=14, height=8)
  print(plt)
  dev.off()


  # collection strategy
  plt <- ggplot(obj@meta.data) +
    aes(x = seurat_clusters, fill = collection_source) +
    geom_bar(position='fill', color='black') +
    theme_classic() + theme +
    labs(title = res, x = "cluster") +
    scale_y_continuous(labels = scales::percent)


  png(paste0("figures/merfish/BAYSOR/proportion_barplots/", run_prefix, '.proportion_barplot.collection_type.', res, ".png"), width=14, height=8, units = 'in', res=300)
  print(plt)
  dev.off()


  svglite::svglite(paste0("figures/merfish/BAYSOR/proportion_barplots/", run_prefix, '.proportion_barplot.collection_type.', res, ".svg"), width=14, height=8)
  print(plt)
  dev.off()


  plt <- DotPlot(obj, group.by = res, features = CellTypeMarkerList,
                 cols = global_colors, col.min = 0) +
    labs(title = "Cell Type Markers", x = NULL, y = res) +  theme +
    rotate_x_text() +
    geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
    scale_size_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100)) +
    guides(color = guide_colorbar(title = "Average\nExpression", frame.colour = 'black'),
           size = guide_legend(title = "Percent\nExpressed")) +
    theme(legend.position = 'bottom',
          panel.grid.minor = element_line(color = "black", size = 0.5))


  png(paste0("figures/merfish/BAYSOR/dot_plots/", run_prefix, '.cell_type_markers.', res, ".png"), width=35, height=10, units = 'in', res=300)
  print(plt)
  dev.off()

  svglite::svglite(paste0("figures/merfish/BAYSOR/dot_plots/", run_prefix, '.cell_type_markers.', res, ".svg"), width=14, height=8)
  print(plt)
  dev.off()

  markerdf <- cluster.marker.list[[res]]


  top.markers <- markerdf %>%
    filter(is.signif) %>%
    group_by(cluster) %>%
    top_n(5, wt=avg_log2FC) %>%
    ungroup() %>%
    arrange(as.numeric(as.character(cluster))) %>%
    dplyr::select(gene) %>%
    deframe() %>%
    unique()

  plt <- DotPlot(obj, group.by = res, features = top.markers, cols = global_colors, col.min = 0) +
    labs(title = res, x = NULL, y = res) + theme + rotate_x_text(angle=45) +
    geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
    scale_size_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100)) +
    guides(color = guide_colorbar(title = "Average\nExpression", frame.colour = 'black'),
           size = guide_legend(title = "Percent\nExpressed")) +
    theme(legend.position = 'bottom', panel.grid.minor = element_line(color = "black", size = 0.5)) +
    coord_equal()

  png(paste0("figures/merfish/BAYSOR/dot_plots/", run_prefix, '.cluster_markers.top_logfc.', res, ".png"), width=36, height=10, units = 'in', res=300)
  print(plt)
  dev.off()


  svglite::svglite(paste0("figures/merfish/BAYSOR/dot_plots/", run_prefix, '.cluster_markers.top_logfc.', res, ".svg"), width=36, height=10)
  print(plt)
  dev.off()

  rm(plt, top.markers, markerdf)
}


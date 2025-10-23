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
setwd("D:/Dropbox/_projects/NormalSkinAtlas")

cluster_prefix <- "reclust.louvain_"
data_path <- "data/scrna/"

resolutions <- c(0.2, 0.4, 0.5, 0.6, 0.8, 1, 1.5, 2, 2.5)
res.cols <- paste0(cluster_prefix, resolutions)
detailed.cell_type.colors <- readRDS("data/reference/cell_type.detailed.color_palette.final.RDS")
global_colors <- c("mistyrose", "maroon2")
cluster_colors <- unique(c(pals::polychrome(n=36), pals::glasbey(n=32)))

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

CellTypeMarkerList <- list(
  Adipo = c("ADIPOQ", "LRP1"),
  Schwann = c("NRXN1", "MPZ"),
  Neuron = c("CACNA1B", "NEFH", "SNAP25", "TAC1", "SLITRK6",'RBFOX1','RBFOX2','RBFOX3','CAMK2A'),
  LEC = c("LYVE1", "ANGPT2", "CCL21"),
  EC = c("VWF", "ACKR1", "FLT1", "KDR"),
  Peri = c("ACTA2","RGS5", "TCIM",  "RGS16","PDGFRB"),
  SM = c("DES", "MCAM", "MYH11"),

  Pan.Fib = c("PDGFRA", "COL1A1"),
  Papil.Fib = c("COMP", "COL18A1", "WIF1"),
  Perivasc.Fib = c("APOE", "CCL19", "PI16", "FN1", "CCDC80", "CXCL12"),
  Retic.Fib.I.II = c("MMP2","SFRP2", "PRG4"),
  Retic.Fib.III.IV = c("THBS4", "CLU", "LUM","ANGPTL7","SFRP1", "APOD"),

  DP = c("COCH", "CRABP1", "FIBIN", "PLXDC1"),
  DS = c("ACAN", "COL11A1", "DPEP1", "LRRC15", "MYL4", "TNN", "ITGA8"),

  HF_HS = c("KRT83", "LGR5", "EPCAM", "LMCD1","KRT32", "KRT35"),
  IRF_IRS_IST = c("PTN", "KRT28", "KRT74", "FST", "PTHLH"),
  ORS.Bulge = c("SOX9", "CD34", "DIO2", "RCAN1", "THBS1"),

  Seb = c("KRT79", "CTSB", "APOC1", "ASS1", "SCARB1", "TKT"),
  Ecc.Duct = c("KRT77", "SEMA3C", "GJB6", "TBX3","FGF7"),
  Ecc.Gland = c("DNER", "TSPAN8", "PLA2G2A","CLDN3", "CLDN10", "ABI3BP"),

  Cyc = c("MKI67", "PCNA", "POLA1", "TOP2A"),

  Bas.KC = c("COL17A1", "COL7A1", "DST", "ITGA6", "ITGB1", "ITGB4"),
  Spn.KC = c("MT1X", "MT1E", "S100A8", "DSC1", "GSTA3", "HES1", "NOTCH3", "RHOV", "SERPINB4"),
  Grn.KC = c("FLG","KLK7", "IL1RN"),
  Melano = c("FRZB", "MITF", "MLANA", "SOX10", "TYR"),

  T.Cell=c("CD3E","CD3D","CD3G","HLA-DRB5", "CD4", "CD8B","CD8A"), # control
  NK=c('NKG7', "NCR1", "GNLY"),
  Cytotox=c("GZMK","GZMA","GZMH","GZMB","GZMM"), #activation
  Exh = c('PDCD1','CTLA4'),
  Naive =c('IL2RG','SELL'), # IL2RG=CD132; CD62L=SELL
  Tolerance=c("NR4A2","NR4A3","NR4A1"),
  T.Naive=c("IL7R","CCR7","S100A4"), # Naive CD4+ T # Memory CD4+
  T.Reg = c('FOXP3','IL10'),
  T.Help = c("IL9", "RORC", "STAT2", "STAT6", "TBX21"),

  B.Cell=c("CD19","MS4A1","CD27","CD38"), # B cells
  B.Naive=c("IL4R","IGHD","CD69","IGHM"), # Naive
  Plasma=c("MZB1","IGHG1","IGHG2","IGHG3","IGHG4","IGHA1","IGHA2","JCHAIN",'XBP1'), #plasma

  ILC1=c('IL6R','CXCR3','BCL11B'),
  ILC2=c('KLRG1','GATA3','PTGDR2','SLAMF1'),
  ILC3=c('IL23R','XCL1','XCL2'), #IL3 when CD69+ and IL7R+

  DC = c("CD80", "BATF3", "FLT3", "IRF4", "IL3RA", "TCL1A", "IRF7",
         "CCL17", "CCL4", "CD200", "IDO1", "CD1C", "CLEC10A", "CLEC9A", "IRF8"),
  LC = c("CD207", "CD1A", "S100B"),

  CD14.Mono = c("CD14", "C5AR1", "LYZ"),
  FCGR3A.Mono = c("FCGR3A", "MS4A7"),
  Mac = c("C1QA", "C1QB", "CD68", "CSF1R", "MRC1")
)
# obj <- qs::qread("data/scrna/seurat_objects/normal_skin.scRNA.harmony.integrated.seurat_object.merfish.annotated.QS")
# obj[['RNA']] <- JoinLayers(obj[['RNA']])
# meta <- obj@meta.data
gc()

categories <- c("Epithelia", "Stroma", "Immune")
print("---- STARTING PLOTS ----")
for (category in categories){


  print(paste("PLOTS FOR:", category))

  reclust.obj <- readRDS(paste0("data/scrna/reclustering/", category, ".reclustered.seurat_object.RDS"))
  run_prefix <- make.names(category)

  # reclust.obj <- AddMetaData(reclust.obj, obj@meta.data, col.name = "sample_barcode")


  DefaultAssay(reclust.obj) <- "RNA"
  reclust.obj[['RNA']] <- JoinLayers(reclust.obj[['RNA']])
  reclust.obj[['RNA']]$data <- as(reclust.obj[['RNA']]$data, Class = "dgCMatrix")

  data_path <- paste0("data/scrna/cluster_markers/", category, "/")
  cluster.marker.list <- readRDS(paste0(data_path, category, '.reclustering.cluster_marker_list.RDS'))

  figure_path <- paste0("figures/scrna/reclustering/", make.names(category), "/")
  dir.create(figure_path, showWarnings = FALSE, recursive = TRUE)
  dir.create(paste0(figure_path, "cluster_umap/"), showWarnings = FALSE, recursive = TRUE)
  dir.create(paste0(figure_path, "proportion_barplots/"), showWarnings = FALSE, recursive = TRUE)
  dir.create(paste0(figure_path, "dot_plots/"), showWarnings = FALSE, recursive = TRUE)
  dir.create(paste0(figure_path, "summary_panel/"), showWarnings = FALSE, recursive = TRUE)

  umap.panel <- DimPlot(reclust.obj, group.by = res.cols, ncol = ceiling(length(res.cols)/2), alpha = 0.25, reduction = 'reclust.umap', raster = TRUE, label = TRUE) & NoAxes() & theme & spatial_theme & coord_fixed()
  png(paste0(figure_path, "summary_panel/", run_prefix, ".cluster_umap.all_resolutions.summary_panel.png"), width=24, height=10, units = 'in', res=600)
  print(umap.panel)
  dev.off()

  prev.annot.panel <- DimPlot(reclust.obj, group.by = c("cell_type.broad", "cell_type.reclustered"), reduction = 'reclust.umap', label = TRUE, ncol=2) & NoAxes() & theme & spatial_theme & coord_fixed()
  png(paste0(figure_path, "summary_panel/", run_prefix, ".cluster_umap.previous_annotations.summary_panel.png"), width=14, height=6, units = 'in', res=600)
  print(prev.annot.panel)
  dev.off()

  ### clustering plots
  for (res in res.cols){

    ## first make the plot itself
    umap.plt <- DimPlot(reclust.obj, reduction = 'reclust.umap', group.by = res, label = FALSE, raster = FALSE, alpha = 0.25)

    ## then add the cluster labels, theme and color scheme
    umap.plt <- LabelClusters(umap.plt, id = res, size = 4, fontface = 'bold', repel = FALSE) +
      theme + spatial_theme


    png(paste0(figure_path, "cluster_umap/", run_prefix, '.cluster_umap.', res, ".png"), width=10, height=8, units = 'in', res=600)
    print(umap.plt)
    dev.off()

    svglite::svglite(paste0(figure_path, "cluster_umap/", run_prefix, '.cluster_umap.', res, ".svg"), width=10, height=8)
    print(umap.plt)
    dev.off()

    ct.plt <- DotPlot(reclust.obj, group.by = res, features = CellTypeMarkerList,
                      cols = global_colors, col.min = 0) +
      labs(title = "Cell Type Markers", x = NULL, y = res) +  theme +
      rotate_x_text() +
      geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
      scale_size_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100)) +
      guides(color = guide_colorbar(title = "Average\nExpression", frame.colour = 'black'),
             size = guide_legend(title = "Percent\nExpressed")) +
      theme(legend.position = 'bottom',
            panel.grid.minor = element_line(color = "black", size = 0.5))

    svglite::svglite(paste0(figure_path, "dot_plots/", run_prefix, '.cell_type_markers.', res, ".svg"), width=35, height=10)
    print(ct.plt)
    dev.off()

    png(paste0(figure_path, "dot_plots/",run_prefix, '.cell_type_markers.', res, ".png"), width=35, height=10, units = 'in', res=600)
    print(ct.plt)
    dev.off()

    markerdf <- cluster.marker.list[[res]] %>%  filter(is.signif & is.positive)
    top.markers <- markerdf %>%
      filter(is.signif & is.positive) %>%
      group_by(group) %>%
      top_n(5, wt=logFC) %>%
      ungroup() %>%
      arrange(as.numeric(as.character(group))) %>%
      dplyr::select(feature) %>%
      deframe() %>%
      unique()

    tm.plt <- DotPlot(reclust.obj, group.by = res, features = top.markers, cols = global_colors, col.min = 0) +
      labs(title = res, x = NULL, y = res) + theme + rotate_x_text() +
      geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
      scale_size_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100)) +
      guides(color = guide_colorbar(title = "Average\nExpression", frame.colour = 'black'),
             size = guide_legend(title = "Percent\nExpressed")) +
      theme(legend.position = 'bottom',
            panel.grid.minor = element_line(color = "black", size = 0.5))

    svglite::svglite(paste0(figure_path, "dot_plots/", run_prefix, '.cluster_markers.top_auc.', res, ".svg"), width=28, height=10)
    print(tm.plt)
    dev.off()

    png(paste0(figure_path, "dot_plots/", run_prefix, '.cluster_markers.top_auc.', res, ".png"), width=28, height=10, units = 'in', res=600)
    print(tm.plt)
    dev.off()

    pn <- DimPlot(reclust.obj, group.by = c("cell_type",'cell_type.reclustered',  res), reduction = 'reclust.umap', label = TRUE, ncol = 3) & NoAxes() & theme & spatial_theme

    recl.panel <- (pn | (tm.plt + NoLegend()))
    recl.panel <- (recl.panel / ct.plt)


    png(paste0(figure_path, "summary_panel/", run_prefix, '.', res, ".summary_panel.png"), width=36, height=12, units = 'in', res=600)
    print(recl.panel)
    dev.off()


    svglite::svglite(paste0(figure_path, "summary_panel/", run_prefix, '.', res, ".summary_panel.svg"), width=36, height=12)
    print(recl.panel)
    dev.off()


    ### Proportion barplots
    reclust.obj@meta.data$seurat_clusters <- reclust.obj@meta.data[, res]

    # donor id
    plt <- ggplot(reclust.obj@meta.data) +
      aes(x = seurat_clusters, fill = donor_id) +
      geom_bar(position='fill', color='black') +
      theme_classic() + theme + labs(title = res, x = "cluster") +
      scale_y_continuous(labels = scales::percent)

    svglite::svglite(paste0(figure_path, "proportion_barplots/", run_prefix, '.proportion_barplot.donor_id.', res, ".svg"), width=14, height=8)
    print(plt)
    dev.off()

    png(paste0(figure_path, "proportion_barplots/", run_prefix, '.proportion_barplot.donor_id.', res, ".png"), width=14, height=8, units = 'in', res=600)
    print(plt)
    dev.off()



    # sample barcode
    plt <- ggplot(reclust.obj@meta.data) +
      aes(x = seurat_clusters, fill = sample_barcode) +
      geom_bar(position='fill', color='black') +
      theme_classic() + theme +
      labs(title = res, x = "cluster") +
      scale_y_continuous(labels = scales::percent)


    svglite::svglite(paste0(figure_path, "proportion_barplots/", run_prefix, '.proportion_barplot.sample_barcode.', res, ".svg"), width=14, height=8)
    print(plt)
    dev.off()

    png(paste0(figure_path, "proportion_barplots/", run_prefix, '.proportion_barplot.sample_barcode.', res, ".png"), width=14, height=8, units = 'in', res=600)
    print(plt)
    dev.off()

#
#
#     # semisupervised cell type predicted ID
#     plt <- ggplot(reclust.obj@meta.data) +
#       aes(x = seurat_clusters, fill = cell_type.detailed.reclustered.cca.predicted_id) +
#       geom_bar(position='fill', color='black') +
#       theme_classic() + theme +
#       labs(title = res, x = "cluster") +
#       scale_fill_manual(values = detailed.cell_type.colors) +
#       scale_y_continuous(labels = scales::percent)
#
#     svglite::svglite(paste0(figure_path, "proportion_barplots/", run_prefix, '.proportion_barplot.predicted_id.semisupervised.', res, ".svg"), width=14, height=8)
#     print(plt)
#     dev.off()
#
#
#     png(paste0(figure_path, "proportion_barplots/", run_prefix, '.proportion_barplot.predicted_id.semisupervised.', res, ".png"), width=14, height=8, units = 'in', res=600)
#     print(plt)
#     dev.off()
#


    # imaging run batch
    plt <- ggplot(reclust.obj@meta.data) +
      aes(x = seurat_clusters, fill = study_id) +
      geom_bar(position='fill', color='black') +
      theme_classic() + theme +
      labs(title = res, x = "cluster") +
      scale_y_continuous(labels = scales::percent)

    svglite::svglite(paste0(figure_path, "proportion_barplots/", run_prefix, '.proportion_barplot.batch.', res, ".svg"), width=14, height=8)
    print(plt)
    dev.off()


    png(paste0(figure_path, "proportion_barplots/", run_prefix, '.proportion_barplot.batch.', res, ".png"), width=14, height=8, units = 'in', res=600)
    print(plt)
    dev.off()

    # sex
    plt <- ggplot(reclust.obj@meta.data) +
      aes(x = seurat_clusters, fill = donor_sex) +
      geom_bar(position='fill', color='black') +
      theme_classic() + theme +
      labs(title = res, x = "cluster") +
      scale_fill_manual(values = c("Male"="dodgerblue", "Female"="maroon2")) +
      scale_y_continuous(labels = scales::percent)

    svglite::svglite(paste0(figure_path, "proportion_barplots/", run_prefix, '.proportion_barplot.donor_sex.', res, ".svg"), width=14, height=8)
    print(plt)
    dev.off()

    png(paste0(figure_path, "proportion_barplots/", run_prefix, '.proportion_barplot.donor_sex.', res, ".png"), width=14, height=8, units = 'in', res=600)
    print(plt)
    dev.off()

    # anatomic site
    plt <- ggplot(reclust.obj@meta.data) +
      aes(x = seurat_clusters, fill = anatomic_site) +
      geom_bar(position='fill', color='black') +
      theme_classic() + theme +
      labs(title = res, x = "cluster") +
      scale_y_continuous(labels = scales::percent)

    svglite::svglite(paste0(figure_path, "proportion_barplots/", run_prefix, '.proportion_barplot.anatomic_site.', res, ".svg"), width=14, height=8)
    print(plt)
    dev.off()

    png(paste0(figure_path, "proportion_barplots/", run_prefix, '.proportion_barplot.anatomic_site.', res, ".png"), width=14, height=8, units = 'in', res=600)
    print(plt)
    dev.off()

    rm(plt, top.markers, markerdf)

   }
  rm(reclust.obj); gc()
}


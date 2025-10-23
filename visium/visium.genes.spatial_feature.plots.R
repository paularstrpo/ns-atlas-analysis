library(tidyverse)
library(Matrix)
library(Seurat)
library(SeuratObject)
library(BPCells)
library(circlize)
library(RColorBrewer)
library(ggrepel)
library(ggpubr)
library(patchwork)
library(ggrastr)
library(svglite)

setwd("D:/Dropbox/_projects/NormalSkinAtlas")

# dotplot colors
global_colors <- c("lightcyan", "royalblue2")

# cluster colors
cluster_colors <- unique(c(pals::polychrome(n=36), pals::glasbey(n=32)))
names(cluster_colors) <- as.character(c(1:length(cluster_colors))-1)

# i/o settings
data_path <- "data/visium/"
run_prefix <- "pan-skin_merged.visium_data"

# cluster settings
resolutions <- c(0.2, 0.4, 0.5, 0.6, 0.8, 1, 1.2, 1.5, 2, 2.5)
cluster_prefix <- "harmony.snn_res."
res.cols <- paste0(cluster_prefix, resolutions)

# load object
obj <- qs::qread('data/visium/seurat_objects/pan-skin_merged.visium_data.harmony_integrated.seurat_object.QS')

# get sample level info
sample_barcodes <- unique(obj@meta.data$sample_id)
sample_table <- unique(obj@meta.data[, c("sample_id","image_id", "donor_id", "study_id", "tissue_type", "condition", "disease", "anatomic_site", "donor_age", "donor_sex")])
sample_table$image_id[is.na(sample_table$image_id)] <- sample_table$sample_id[is.na(sample_table$image_id)]
# select sample subset types to plot
ns.samples <- sample_table$sample_id[sample_table$tissue_type == 'NS']
ad.samples <- sample_table$sample_id[sample_table$tissue_type == 'AD']
pp.samples <- sample_table$sample_id[sample_table$tissue_type == 'PP']
hs.samples <- sample_table$sample_id[sample_table$tissue_type == 'HS']

inf.samples <- c(pp.samples, ad.samples, hs.samples)

scc.samples <- sample_table$sample_id[sample_table$tissue_type == 'SCC']
bcc.samples <- sample_table$sample_id[sample_table$tissue_type == 'BCC']
tumor.samples <- c(scc.samples, bcc.samples)

# spatial_genes <- c("KRT6A","KRT6B","KRT6C","S100A8","S100A9","IGLC1","IGHM","JCHAIN","CD79A","CXCL12","CXCL13","CCL19","CCL5","CCL13","CCL18","CXCL9",
#                   "CXCL10","PTCH1","PTCH2","SMO","EPCAM","TRAC","TRBC1","TRBC2","CD2","CCR7","CCL22","CCL17","CD3E","TNC")
# spatial_genes <- c("KRT13", "MYC", "EHF", "KLK8", "ST6GALNAC2", "HIST1H4C", "IGHG2", "IGLC1")
spatial_genes <- c("MZB1", "IGLV3-1", "IGLV6-57", "IGHG2")
### do all neighborhoods
for (gene in spatial_genes){
  print(gene)


  plot_list <- lapply(sample_barcodes, function(sample){
    image_id = make.names(as.character(sample_table[sample_table$sample_id == sample, "image_id"]))
    SpatialFeaturePlot(obj, images = image_id, features = gene, crop = TRUE) +
      theme(title = element_text(color = 'black', hjust = 0.5), legend.position='right') +
      labs(title = sample)
  })
  names(plot_list) <- sample_barcodes

  spatial_cluster_panel <- cowplot::plot_grid(plotlist = plot_list, ncol=10, align = 'hv')

  svglite(paste0("figures/visium/spatial_feature/", run_prefix, '.all_samples.', gene, ".spatial_featureplot.svg"), width=48, height=32)
  print(spatial_cluster_panel)
  dev.off()

  png(paste0("figures/visium/spatial_feature/", run_prefix, '.all_samples.', gene, ".spatial_featureplot.png"), width=48, height=32, units = 'in', res=300)
  print(spatial_cluster_panel)
  dev.off()


  ## normal skin only
  spatial_cluster_panel <- cowplot::plot_grid(plotlist = plot_list[ns.samples], ncol=6, align = 'hv')

  svglite(paste0("figures/visium/spatial_feature/", run_prefix, '.normal_skin.', gene, ".spatial_featureplot.svg"), width=48, height=32)
  print(spatial_cluster_panel)
  dev.off()

  png(paste0("figures/visium/spatial_feature/", run_prefix, '.normal_skin.', gene, ".spatial_featureplot.png"), width=48, height=32, units = 'in', res=300)
  print(spatial_cluster_panel)
  dev.off()


  ## AD Lesional Samples
  spatial_cluster_panel <- cowplot::plot_grid(plotlist = plot_list[ad.samples], ncol=3, align = 'hv')

  svglite(paste0("figures/visium/spatial_feature/", run_prefix, '.atopic-ls.', gene, ".spatial_featureplot.svg"), width=32, height=16)
  print(spatial_cluster_panel)
  dev.off()

  png(paste0("figures/visium/spatial_feature/", run_prefix, '.atopic-ls.', gene, ".spatial_featureplot.png"), width=32, height=16, units = 'in', res=300)
  print(spatial_cluster_panel)
  dev.off()


  ## PS Lesional Samples
  spatial_cluster_panel <- cowplot::plot_grid(plotlist = plot_list[pp.samples], ncol=2, align = 'hv')

  svglite(paste0("figures/visium/spatial_feature/", run_prefix, '.psoriasis-ls.', gene, ".spatial_featureplot.svg"), width=32, height=16)
  print(spatial_cluster_panel)
  dev.off()

  png(paste0("figures/visium/spatial_feature/", run_prefix, '.psoriasis-ls.', gene, ".spatial_featureplot.png"), width=32, height=16, units = 'in', res=300)
  print(spatial_cluster_panel)
  dev.off()

  ## HS Samples
  spatial_cluster_panel <- cowplot::plot_grid(plotlist = plot_list[hs.samples], ncol=3, align = 'hv')

  svglite(paste0("figures/visium/spatial_feature/", run_prefix, '.cluster_spatial.hs.', gene, ".spatial_featureplot.svg"), width=32, height=16)
  print(spatial_cluster_panel)
  dev.off()

  png(paste0("figures/visium/spatial_feature/", run_prefix, '.cluster_spatial.hs.', gene, ".spatial_featureplot.png"), width=32, height=16, units = 'in', res=300)
  print(spatial_cluster_panel)
  dev.off()

  ## all inf disease samples
  spatial_cluster_panel <- cowplot::plot_grid(plotlist = plot_list[inf.samples], ncol=4, align = 'hv')

  svglite(paste0("figures/visium/spatial_feature/", run_prefix, '.all_inflm.', gene, ".spatial_featureplot.svg"), width=32, height=16)
  print(spatial_cluster_panel)
  dev.off()

  png(paste0("figures/visium/spatial_feature/", run_prefix, '.all_inflm.', gene, ".spatial_featureplot.png"), width=32, height=16, units = 'in', res=300)
  print(spatial_cluster_panel)
  dev.off()

  ## BCC Lesional Samples
  spatial_cluster_panel <- cowplot::plot_grid(plotlist = plot_list[bcc.samples], ncol=4, align = 'hv')


  svglite(paste0("figures/visium/spatial_feature/", run_prefix, '.bcc-tumor.', gene, ".spatial_featureplot.svg"), width=32, height=16)
  print(spatial_cluster_panel)
  dev.off()

  png(paste0("figures/visium/spatial_feature/", run_prefix, '.bcc-tumor.', gene, ".spatial_featureplot.png"), width=32, height=16, units = 'in', res=300)
  print(spatial_cluster_panel)
  dev.off()


  ## SCC Lesional Samples
  spatial_cluster_panel <- cowplot::plot_grid(plotlist =  plot_list[scc.samples], ncol=4, align = 'hv')


  svglite(paste0("figures/visium/spatial_feature/", run_prefix, '.scc-tumor.', gene, ".spatial_featureplot.svg"), width=32, height=16)
  print(spatial_cluster_panel)
  dev.off()

  png(paste0("figures/visium/spatial_feature/", run_prefix, '.scc-tumor.', gene, ".spatial_featureplot.png"), width=32, height=16, units = 'in', res=300)
  print(spatial_cluster_panel)
  dev.off()


  ## All Tumor Samples
  spatial_cluster_panel <- cowplot::plot_grid(plotlist = plot_list[tumor.samples], ncol=4, align = 'hv')

  svglite(paste0("figures/visium/spatial_feature/", run_prefix, '.all_tumor.', gene, ".spatial_featureplot.svg"), width=32, height=16)
  print(spatial_cluster_panel)
  dev.off()

  png(paste0("figures/visium/spatial_feature/", run_prefix, '.all_tumor.', gene, ".spatial_featureplot.png"), width=32, height=16, units = 'in', res=300)
  print(spatial_cluster_panel)
  dev.off()

  rm(plot_list, spatial_cluster_panel); gc()
}

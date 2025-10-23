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

# neighborhood colors
neighborhood_colors <- readRDS("data/reference/neighborhood.color_palette.final.RDS")
neighborhoods <- paste0("N", (1:length(neighborhood_colors)-1))
names(neighborhood_colors) <- neighborhoods

# i/o settings
data_path <- "data/visium/"
run_prefix <- "pan-skin_merged.visium_data"
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
### do all spatial cluster plots
for (r in resolutions){
  res <- paste0(cluster_prefix, r)
  print(res)

  clusters <- levels(obj@meta.data[, res])
  obj@meta.data$seurat_clusters <- obj@meta.data[, res]
  Idents(obj) <- "seurat_clusters"


  for (cluster in clusters){
    print(cluster)

    dir.create(path = paste0("figures/visium/cluster_spatial/spatial_highlight/", res), showWarnings = FALSE, recursive = TRUE)

    cluster_cells <- rownames(obj@meta.data[obj@meta.data$seurat_clusters == cluster,])

    plot_list <- lapply(sample_barcodes, function(sample){
      image_id = make.names(as.character(sample_table[sample_table$sample_id == sample, "image_id"]))
      SpatialDimPlot(obj, images = image_id, cells.highlight = cluster_cells, crop = TRUE, cols.highlight = c(deframe(cluster_colors[cluster]), "gray90")) +
        labs(title = sample) +
        theme(title = element_text(face = 'bold', color = 'black', hjust = 0.5)) +
        NoLegend()
    })
    names(plot_list) <- sample_barcodes

    spatial_cluster_panel <- cowplot::plot_grid(plotlist = plot_list, ncol=10, align = 'hv')

    file_prefix <- paste0("figures/visium/cluster_spatial/spatial_highlight/", res, "/", run_prefix, ".", res, "_C", cluster, ".spatial_highlight.")

    svglite(paste0(file_prefix, 'all_samples', ".svg"), width=48, height=32)
    print(spatial_cluster_panel)
    dev.off()

    png(paste0(file_prefix, 'all_samples', ".png"), width=48, height=32, units = 'in', res=300)
    print(spatial_cluster_panel)
    dev.off()


    ## normal skin only
    spatial_cluster_panel <- cowplot::plot_grid(plotlist = plot_list[ns.samples], ncol=6, align = 'hv')

    svglite(paste0(file_prefix, 'normal_skin', ".svg"), width=48, height=32)
    print(spatial_cluster_panel)
    dev.off()

    png(paste0(file_prefix, 'normal_skin', ".png"), width=48, height=32, units = 'in', res=300)
    print(spatial_cluster_panel)
    dev.off()


    ## HS Samples
    spatial_cluster_panel <- cowplot::plot_grid(plotlist = plot_list[hs.samples], ncol=3, align = 'hv')

    svglite(paste0(file_prefix, 'hs', ".svg"), width=32, height=16)
    print(spatial_cluster_panel)
    dev.off()

    png(paste0(file_prefix, 'hs', ".png"), width=32, height=16, units = 'in', res=300)
    print(spatial_cluster_panel)
    dev.off()

    ## AD Lesional Samples
    spatial_cluster_panel <- cowplot::plot_grid(plotlist = plot_list[ad.samples], ncol=3, align = 'hv')

    svglite(paste0(file_prefix, 'atopic-ls', ".svg"), width=32, height=16)
    print(spatial_cluster_panel)
    dev.off()

    png(paste0(file_prefix, 'atopic-ls', ".png"), width=32, height=16, units = 'in', res=300)
    print(spatial_cluster_panel)
    dev.off()


    ## PS Lesional Samples
    spatial_cluster_panel <- cowplot::plot_grid(plotlist = plot_list[pp.samples], ncol=2, align = 'hv')

    svglite(paste0(file_prefix, 'psoriasis-ls', ".svg"), width=18, height=18)
    print(spatial_cluster_panel)
    dev.off()

    png(paste0(file_prefix, 'psoriasis-ls', ".png"), width=18, height=18, units = 'in', res=300)
    print(spatial_cluster_panel)
    dev.off()


    ## all inf disease samples
    spatial_cluster_panel <- cowplot::plot_grid(plotlist = plot_list[inf.samples], ncol=4, align = 'hv')

    svglite(paste0(file_prefix, 'all_inflm', ".svg"), width=36, height=18)
    print(spatial_cluster_panel)
    dev.off()

    png(paste0(file_prefix, 'all_inflm', ".png"), width=36, height=18, units = 'in', res=300)
    print(spatial_cluster_panel)
    dev.off()

    ## BCC Lesional Samples
    spatial_cluster_panel <- cowplot::plot_grid(plotlist = plot_list[bcc.samples], ncol=4, align = 'hv')

    svglite(paste0(file_prefix, 'bcc-tumor', ".svg"), width=36, height=16)
    print(spatial_cluster_panel)
    dev.off()

    png(paste0(file_prefix, 'bcc-tumor', ".png"), width=36, height=16, units = 'in', res=300)
    print(spatial_cluster_panel)
    dev.off()

    ## SCC Lesional Samples
    spatial_cluster_panel <- cowplot::plot_grid(plotlist =  plot_list[scc.samples], ncol=4, align = 'hv')

    svglite(paste0(file_prefix, 'scc-tumor', ".svg"), width=36, height=16)
    print(spatial_cluster_panel)
    dev.off()

    png(paste0(file_prefix, 'scc-tumor', ".png"), width=36, height=16, units = 'in', res=300)
    print(spatial_cluster_panel)
    dev.off()


    ## All Tumor Samples
    spatial_cluster_panel <- cowplot::plot_grid(plotlist = plot_list[tumor.samples], ncol=4, align = 'hv')

    svglite(paste0(file_prefix, 'all_tumor', ".svg"), width=36, height=36)
    print(spatial_cluster_panel)
    dev.off()

    png(paste0(file_prefix, 'all_tumor', ".png"), width=36, height=36, units = 'in', res=300)
    print(spatial_cluster_panel)
    dev.off()

    rm(plot_list, spatial_cluster_panel); gc()
  }

}

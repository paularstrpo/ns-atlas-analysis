root_path <- "D:/Dropbox/_projects/NormalSkinAtlas"
setwd(root_path)

library(tidyverse)
library(Matrix)
library(Seurat)
library(SeuratObject)

library(ggpubr)
library(ggrepel)
library(ggbeeswarm)

library(viridis)
library(circlize)
library(RColorBrewer)
library(paletteer)
library(patchwork)

obj <- qs::qread("data/merfish/BAYSOR/seurat_objects/ns-atlas.merfish_baysor.scanvi_intergated.cellcharter.seurat_object.QS", nthreads=8)
compute_angle <- function(perc) {
  angle = -1
  if(perc < 0.25) # 1st q [90,0]
    angle = 90 - (perc/0.25) * 90
  else if(perc < 0.5) # 2nd q [0, -90]
    angle = (perc-0.25) / 0.25 * -90
  else if(perc < 0.75) # 3rd q [90, 0]
    angle = 90 - ((perc-0.5) / 0.25 * 90)
  else if(perc < 1.00) # last q [0, -90]
    angle = ((perc -0.75)/0.25) * -90
  # Or even more compact, but less readable
  if(perc < 0.5) # 1st half [90, -90]
    angle = (180 - (perc/0.5) * 180) - 90
  else # 2nd half [90, -90]
    angle = (90 - ((perc - 0.5)/0.5) * 180)
  return(angle)
}

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

broad_celltype.colors <- Polychrome::light.colors(n = length(levels(obj@meta.data$cell_type.broad)))
names(broad_celltype.colors) <- levels(obj@meta.data$cell_type.broad)
detailed_celltype.colors <- readRDS("data/reference/cell_type.detailed.color_palette.varibow.non_shuffled.RDS")

metadf <- obj@meta.data

root_path <- "figures/merfish/BAYSOR/cellcharter/"
dir.create(root_path, showWarnings = FALSE, recursive = TRUE)

maxima <- c(2, 3, 10)
res.cols <- paste0("cellcharter_cluster.k_", maxima)

for (res in res.cols){
  print(res)
  figure_path <- paste0(root_path, res, '/')
  dir.create(figure_path, showWarnings = FALSE, recursive = TRUE)
  dir.create(paste0(figure_path, "proportion_barplots/"), showWarnings = FALSE, recursive = TRUE)
  dir.create(paste0(figure_path, "spatial_highlight/"), showWarnings = FALSE, recursive = TRUE)

  obj$cellcharter_domain <- obj@meta.data[, res]
  metadf$cellcharter_domain <- metadf[, res]
  cell_counts <- table(metadf$cellcharter_domain)

  nhood.cols <- colorway::varibow(length(unique(obj$cellcharter_domain)))
  names(nhood.cols) <- 0:(length(unique(obj$cellcharter_domain))-1)


  domain_spatial <- ggplot(metadf) +
    aes(x = spatial_1, y = spatial_2, color = cellcharter_domain) +
    geom_point(size = 0.25) +
    facet_wrap(~ sample_barcode, scales = 'free') +
    theme_classic() + spatial_theme +
    scale_color_manual(values = nhood.cols) +
    guides(color = guide_legend(override.aes = list(size = 5), title = 'CellCharter\nDomain'))

  file_prefix <- paste0(figure_path, "ns-atlas.merfish.cellcharter.", res, ".neighborhood_spatial.sample_barcode")

  png(filename = paste0(file_prefix, ".png"), height = 35, width=30, units = 'in', res = 150)
  print(domain_spatial)
  dev.off()

  ## proportion barplots

  ## domain by sample
  plt <- ggplot(metadf) +
    aes(x = sample_barcode, fill = cellcharter_domain) +
    geom_bar(position='fill', color='black') +
    theme_classic() +
    scale_y_continuous(labels = scales::percent) +
    scale_fill_manual(values = nhood.cols) +
    guides(fill = guide_legend(title = 'CellCharter\nDomain')) +
    spatial_theme + theme(axis.text=element_text(size=10, color='black')) +
    rotate_x_text()

  file_prefix <- paste0(figure_path, "proportion_barplots/ns-atlas.merfish.cellcharter.proportion_barplot.sample_barcode.", res)
  svglite::svglite(paste0(file_prefix, ".svg"), width=36, height=8)
  print(plt)
  dev.off()

  png(paste0(file_prefix, ".png"), width=36, height=8, units = 'in', res = 150)
  print(plt)
  dev.off()


  ## domain by sample
  plt <- ggplot(metadf) +
    aes(x = donor_id, fill = cellcharter_domain) +
    geom_bar(position='fill', color='black') +
    theme_classic() +
    scale_y_continuous(labels = scales::percent) +
    scale_fill_manual(values = nhood.cols) +
    guides(fill = guide_legend(title = 'CellCharter\nDomain')) +
    spatial_theme + theme(axis.text=element_text(size=10, color='black')) +
    rotate_x_text()

  file_prefix <- paste0(figure_path, "proportion_barplots/ns-atlas.merfish.cellcharter.proportion_barplot.donor_id.", res)
  svglite::svglite(paste0(file_prefix, ".svg"), width=36, height=8)
  print(plt)
  dev.off()

  png(paste0(file_prefix, ".png"), width=36, height=8, units = 'in', res = 150)
  print(plt)
  dev.off()

  ## gene panel
  plt <- ggplot(metadf) +
    aes(x = cellcharter_domain, fill = gene_panel) +
    geom_bar(position='fill', color='black') +
    theme_classic() + labs(x = NULL) +
    scale_y_continuous(labels = scales::percent) +
    spatial_theme + theme(axis.text=element_text(size=10, color='black'))

  file_prefix <- paste0(figure_path, "proportion_barplots/ns-atlas.merfish.cellcharter.proportion_barplot.gene_panel.", res)
  svglite::svglite(paste0(file_prefix, ".svg"), width=14, height=8)
  print(plt)
  dev.off()

  png(paste0(file_prefix, ".png"), width=18, height=8, units = 'in', res = 150)
  print(plt)
  dev.off()


  ## batch
  plt <- ggplot(metadf) +
    aes(x = cellcharter_domain, fill = batch) +
    geom_bar(position='fill', color='black') +
    theme_classic() + labs(x = NULL) +
    scale_y_continuous(labels = scales::percent) +
    spatial_theme + theme(axis.text=element_text(size=10, color='black'))

  file_prefix <- paste0(figure_path, "proportion_barplots/ns-atlas.merfish.cellcharter.proportion_barplot.imaging_batch.", res)
  svglite::svglite(paste0(file_prefix, ".svg"), width=14, height=8)
  print(plt)
  dev.off()

  png(paste0(file_prefix, ".png"), width=18, height=8, units = 'in', res = 150)
  print(plt)
  dev.off()


  ## collection source
  plt <- ggplot(metadf) +
    aes(x = cellcharter_domain, fill = collection_source) +
    geom_bar(position='fill', color='black') +
    theme_classic() + labs(x = NULL) +
    scale_y_continuous(labels = scales::percent) +
    spatial_theme + theme(axis.text=element_text(size=10, color='black'))

  file_prefix <- paste0(figure_path, "proportion_barplots/ns-atlas.merfish.cellcharter.proportion_barplot.collection_source.", res)
  svglite::svglite(paste0(file_prefix, ".svg"), width=14, height=8)
  print(plt)
  dev.off()

  png(paste0(file_prefix, ".png"), width=18, height=8, units = 'in', res = 150)
  print(plt)
  dev.off()

  ## cell type
  plt <- ggplot(metadf) +
    aes(x = cellcharter_domain, fill = cell_type.broad) +
    geom_bar(position='fill', color='black') +
    theme_minimal() +  NoGrid() +
    scale_y_continuous(labels = scales::percent) +
    coord_flip() + labs(x = NULL, y = NULL) +
    theme(axis.text.x = element_text(color = 'black', size=12),
          axis.text.y = element_text(color = 'black', size=12, hjust = 1, vjust = 0.5, face='bold'),
          legend.title = element_text(color = 'black', size=12, face='bold', vjust = 0.5, hjust = 1),
          legend.text = element_text(color = 'black', size=10))

  file_prefix <- paste0(figure_path, "proportion_barplots/ns-atlas.merfish.cellcharter.proportion_barplot.cell_type.broad.", res)
  svglite::svglite(paste0(file_prefix, ".svg"), width=18, height=8)
  print(plt)
  dev.off()

  png(paste0(file_prefix, ".png"), width=18, height=8, units = 'in', res = 150)
  print(plt)
  dev.off()

  ## cell category
  plt <- ggplot(metadf) +
    aes(x = cellcharter_domain, fill =  cell_category) +
    geom_bar(position='fill', color='black') +
    theme_minimal() + NoGrid() +
    scale_y_continuous(labels = scales::percent) +
    coord_flip() + labs(x = NULL, y = NULL) +
    theme(axis.text.x = element_text(color = 'black', size=12),
          axis.text.y = element_text(color = 'black', size=12, hjust = 1, vjust = 0.5, face='bold'),
          legend.title = element_text(color = 'black', size=12, face='bold', hjust = 0.5),
          legend.text = element_text(color = 'black', size=10))


  file_prefix <- paste0(figure_path, "proportion_barplots/ns-atlas.merfish.cellcharter.proportion_barplot.cell_category.", res)
  svglite::svglite(paste0(file_prefix, ".svg"), width=18, height=8)
  print(plt)
  dev.off()

  png(paste0(file_prefix, ".png"), width=18, height=8, units = 'in', res = 150)
  print(plt)
  dev.off()



  ## predicted cell type
  plt <- ggplot(metadf) +
    aes(x = cellcharter_domain, fill= scrna_predicted_id) +
    geom_bar(position='fill', color='black') +
    theme_minimal() + NoGrid() +
    scale_y_continuous(labels = scales::percent) +
    coord_flip() + labs(x = NULL, y = NULL) +
    theme(axis.text.x = element_text(color = 'black', size=12),
          axis.text.y = element_text(color = 'black', size=12, hjust = 1, vjust = 0.5, face='bold'),
          legend.title = element_text(color = 'black', size=12, face='bold', hjust = 0.5),
          legend.text = element_text(color = 'black', size=10))


  file_prefix <- paste0(figure_path, "proportion_barplots/ns-atlas.merfish.cellcharter.proportion_barplot.scrna_predicted_id.", res)
  svglite::svglite(paste0(file_prefix, ".svg"), width=18, height=8)
  print(plt)
  dev.off()

  png(paste0(file_prefix, ".png"), width=18, height=8, units = 'in', res = 150)
  print(plt)
  dev.off()

  ## predicted cell type
  plt <- ggplot(metadf) +
    aes(x = cellcharter_domain, fill= cell_type.detailed) +
    geom_bar(position='fill', color='black') +
    theme_minimal() + NoGrid() +
    scale_y_continuous(labels = scales::percent) +
    scale_fill_manual(values=detailed_celltype.colors)+
    coord_flip() + labs(x = NULL, y = NULL) +
    theme(axis.text.x = element_text(color = 'black', size=12),
          axis.text.y = element_text(color = 'black', size=12, hjust = 1, vjust = 0.5, face='bold'),
          legend.title = element_text(color = 'black', size=12, face='bold', hjust = 0.5),
          legend.text = element_text(color = 'black', size=10))


  file_prefix <- paste0(figure_path, "proportion_barplots/ns-atlas.merfish.cellcharter.proportion_barplot.cell_type.detailed.", res)
  svglite::svglite(paste0(file_prefix, ".svg"), width=18, height=8)
  print(plt)
  dev.off()

  png(paste0(file_prefix, ".png"), width=18, height=8, units = 'in', res = 150)
  print(plt)
  dev.off()


  ## predicted cell type
  plt <- ggplot(metadf) +
    aes(x = cell_type.detailed, fill= cellcharter_domain) +
    geom_bar(position='fill', color='black') +
    theme_minimal() + NoGrid() +
    scale_y_continuous(labels = scales::percent) +
    scale_fill_manual(values=nhood.cols)+
    coord_flip() + labs(x = NULL, y = NULL) +
    theme(axis.text.x = element_text(color = 'black', size=12),
          axis.text.y = element_text(color = 'black', size=12, hjust = 1, vjust = 0.5, face='bold'),
          legend.title = element_text(color = 'black', size=12, face='bold', hjust = 0.5),
          legend.text = element_text(color = 'black', size=10))


  file_prefix <- paste0(figure_path, "proportion_barplots/ns-atlas.merfish.cellcharter.proportion_barplot.cell_type.detailed.flipped.", res)
  svglite::svglite(paste0(file_prefix, ".svg"), width=18, height=8)
  print(plt)
  dev.off()

  png(paste0(file_prefix, ".png"), width=18, height=8, units = 'in', res = 150)
  print(plt)
  dev.off()


  ## anatomic sites
  plt <- ggplot(metadf) +
    aes(x = anatomic_site, fill = cellcharter_domain) +
    geom_bar(position='fill', color='black', alpha=0.7) +
    theme_minimal() + NoGrid() +
    scale_y_continuous(labels = scales::percent) +
    scale_fill_manual(values = nhood.cols) +
    coord_flip() + labs(x = NULL, y = NULL) +
    theme(axis.text.x = element_text(color = 'black', size=12),
          axis.text.y = element_text(color = 'black', size=12, hjust = 1, vjust = 0.5, face='bold'),
          legend.title = element_text(color = 'black', size=12, face='bold', hjust = 0.5),
          legend.text = element_text(color = 'black', size=10))

  file_prefix <- paste0(figure_path, "proportion_barplots/ns-atlas.merfish.cellcharter.proportion_barplot.anatomic_site.", res)

  svglite::svglite(paste0(file_prefix, ".svg"), width=18, height=8)
  print(plt)
  dev.off()

  png(paste0(file_prefix, ".png"), width=18, height=8, units = 'in', res = 150)
  print(plt)
  dev.off()

  reclust.meta <- metadf[, c("cell_barcode", 'sample_barcode', 'batch', 'spatial_1','spatial_2', 'cellcharter_domain')]
  neighborhoods <- levels(metadf$cellcharter_domain)

  ## Neighborhood Composition Donut Plots
  donut_list <- list()
  for (nhood in neighborhoods){
    df <- metadf %>%
      filter(cellcharter_domain == nhood) %>%
      group_by(cellcharter_domain,  cell_category,  cell_type.broad) %>%
      summarise(no.cells = n()) %>%
      ungroup()

    sum_total_cells <- sum(df$no.cells)
    firstLevel <- df %>% group_by(cellcharter_domain) %>% summarize(total.cells=sum(no.cells))

    secondLevel  <- df %>%
      group_by(cell_type.broad) %>%
      summarize(no.cells=sum(no.cells)) %>%
      mutate(cell_type.broad = reorder(cell_type.broad, no.cells),
             pct.cells = no.cells / sum_total_cells) %>%
      arrange(desc(pct.cells)) %>%
      ungroup() %>%
      mutate(running=cumsum(pct.cells), pos = running - pct.cells/2) %>%
      group_by(1:n()) %>% # to compute row by row
      mutate(angle=compute_angle((running - pct.cells/2) / 1)) %>%
      ungroup()


    donut_list[[nhood]] <- ggplot(firstLevel)  +
      geom_bar(data=firstLevel, aes(x=1, y=1), fill='white', stat='identity') +
      geom_bar(data = secondLevel, aes(x = 2, y = pct.cells, fill = cell_type.broad),
               position = 'stack', stat = 'identity', color = 'black') +
      geom_text(data = secondLevel[secondLevel$pct.cells>0.01,],
                aes(label = cell_type.broad, x = 2, y = pos, angle = angle),
                position = 'identity', color = 'black', fontface = 'bold') +
      coord_polar('y') + theme_minimal() + spatial_theme +
      scale_fill_manual(values = broad_celltype.colors)+
      NoAxes() + labs(title = nhood) + NoLegend()
  }

  neighborhood_composition.donut_panel <- patchwork::wrap_plots(donut_list, guides = 'collect', nrow=2)

  file_prefix <- paste0(figure_path, "proportion_barplots/ns-atlas.merfish.cellcharter.donut_panel.cell_type.broad.", res)

  svglite::svglite(paste0(file_prefix, ".svg"), height = 14, width = 36)
  print(neighborhood_composition.donut_panel)
  dev.off()

  png(paste0(file_prefix, ".png"), height = 14, width = 36, units = 'in', res = 150)
  print(neighborhood_composition.donut_panel)
  dev.off()


  ## Neighborhood Composition Donut Plots
  donut_list <- list()
  for (nhood in neighborhoods){
    df <- metadf %>%
      filter(cellcharter_domain == nhood) %>%
      group_by(cellcharter_domain,  cell_category,  cell_type.detailed) %>%
      summarise(no.cells = n()) %>%
      ungroup()

    sum_total_cells <- sum(df$no.cells)
    firstLevel <- df %>% group_by(cellcharter_domain) %>% summarize(total.cells=sum(no.cells))

    secondLevel  <- df %>%
      group_by(cell_type.detailed) %>%
      summarize(no.cells=sum(no.cells)) %>%
      mutate(cell_type.detailed = reorder(cell_type.detailed, no.cells),
             pct.cells = no.cells / sum_total_cells) %>%
      arrange(desc(pct.cells)) %>%
      ungroup() %>%
      mutate(running=cumsum(pct.cells), pos = running - pct.cells/2) %>%
      group_by(1:n()) %>% # to compute row by row
      mutate(angle=compute_angle((running - pct.cells/2) / 1)) %>%
      ungroup()


    donut_list[[nhood]] <- ggplot(firstLevel)  +
      geom_bar(data=firstLevel, aes(x=1, y=1), fill='white', stat='identity') +
      geom_bar(data = secondLevel, aes(x = 2, y = pct.cells, fill =cell_type.detailed),
               position = 'stack', stat = 'identity', color = 'black') +
      geom_text(data = secondLevel[secondLevel$pct.cells>0.01,],
                aes(label = cell_type.detailed, x = 2, y = pos, angle = angle),
                position = 'identity', color = 'black', fontface = 'bold') +
      coord_polar('y') + theme_minimal() + spatial_theme +
      scale_fill_manual(values = detailed_celltype.colors)+
      NoAxes() + labs(title = nhood) + NoLegend()
  }

  neighborhood_composition.donut_panel <- patchwork::wrap_plots(donut_list, guides = 'collect', nrow=2)

  file_prefix <- paste0(figure_path, "proportion_barplots/ns-atlas.merfish.cellcharter.donut_panel.cell_type.detailed.", res)

  svglite::svglite(paste0(file_prefix, ".svg"), height = 10, width = 36)
  print(neighborhood_composition.donut_panel)
  dev.off()

  png(paste0(file_prefix, ".png"), height = 10, width = 36, units = 'in', res = 150)
  print(neighborhood_composition.donut_panel)
  dev.off()

}

library(clustree)
meta$cellcharter_cluster.k_1 <- 0
metadf <- meta[, c("cell_barcode", "sample_barcode", "cell_type.broad","cellcharter_cluster.k_1", "cellcharter_cluster.k_2", "cellcharter_cluster.k_3", "cellcharter_cluster.k_10")]
treeplt <- clustree(metadf, prefix="cellcharter_cluster.k_", prop_filter=0.25)
treeplt

## spatial highlight plots
for (res in res.cols){
  print(res)
  figure_path <- paste0(root_path, res, '/')
  dir.create(figure_path, showWarnings = FALSE, recursive = TRUE)
  dir.create(paste0(figure_path, "proportion_barplots/"), showWarnings = FALSE, recursive = TRUE)
  dir.create(paste0(figure_path, "spatial_highlight/"), showWarnings = FALSE, recursive = TRUE)

  obj$cellcharter_domain <- obj@meta.data[, res]
  metadf$cellcharter_domain <- metadf[, res]
  cell_counts <- table(metadf$cellcharter_domain)

  nhood.cols <- colorway::varibow(length(unique(obj$cellcharter_domain)))
  names(nhood.cols) <- 0:(length(unique(obj$cellcharter_domain))-1)

  reclust.meta <- metadf[, c("cell_barcode", 'sample_barcode', 'batch','spatial_1','spatial_2', 'cellcharter_domain')]
  clusters <- levels(metadf$cellcharter_domain)

  for (cluster in clusters){

    print(cluster)
    nhood <- paste0("N", cluster)
    reclust.meta$to_highlight <- factor(ifelse(reclust.meta$cellcharter_domain == cluster, nhood, "Other"), levels = c("Other", nhood))
    reclust.meta <- reclust.meta %>% arrange(to_highlight)

    cols <- c("gray90", nhood.cols[names(nhood.cols) == cluster])
    names(cols) <- c("Other", nhood)

    plt <- ggplot(reclust.meta) +
      aes(x = spatial_1, y = spatial_2, color = to_highlight) +
      facet_wrap(~ sample_barcode, scales = 'free') +
      geom_point(size = 0.25, show.legend = FALSE) +
      scale_color_manual(values = cols) +
      labs(title = paste0("Neighborhood Spatial Highlight: ", nhood)) +
      theme_classic() + spatial_theme + NoLegend()


    file_prefix <- paste0(figure_path, "spatial_highlight/ns-atlas.merfish.cellcharter.", res, ".spatial_highlight.cluster_", cluster, ".sample_barcode")

    png(paste0(file_prefix, ".png"), height=35, width=30, units = 'in', res = 300)
    print(plt)
    dev.off()

  }
}

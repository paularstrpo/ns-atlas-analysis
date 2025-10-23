# working directory
setwd("D:/Dropbox/_projects/NormalSkinAtlas")

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

## functions

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



## plot colors
global.colors <- c("lightcyan", "royalblue2")

global.cell_type.colors <- c(
  "Bas KC" = "dodgerblue",
  "Spn KC" = "darkorchid1",
  "Grn KC" = "gold",
  "Seb" =  "salmon",
  "Ecc Duct" = "cyan1",
  "Ecc Gland" = "blue",
  "HFE" = "lawngreen",
  "Fibro" = "maroon2",
  "Peri" = "yellow",
  "EC" = "orange",
  "LEC" ="sienna",
  "Adipo" = "mistyrose1",
  "Schwann" = "purple2",
  "Lym" = "darkblue",
  "Plasma" = "darkseagreen1",
  "LC" = "red",
  "Mac/DC" = "darkgreen",
  "Mast" = "maroon4"
)
## load data
obj <- qs::qread("data/merfish/BAYSOR/seurat_objects/ns-atlas.merfish_baysor.scanvi_integrated.cellcharter.annotated.seurat_object.QS")
obj <- UpdateSeuratObject(obj)

obj@meta.data$cell_type.broad <- factor(obj@meta.data$cell_type.broad, levels = names(global.cell_type.colors))
obj@meta.data$neighborhood.detailed <- obj@meta.data$neighborhood
metadf <- obj@meta.data
neighborhood.colors <- pals::glasbey(n=length(levels(obj@meta.data$neighborhood.detailed)))
names(neighborhood.colors) <- levels(obj@meta.data$neighborhood.detailed)
neighborhoods <- levels(obj@meta.data$neighborhood.detailed)

meta <- obj@meta.data %>% filter(cell_type.detailed != "Doublet" & cell_type.detailed != 'Unknown')

d165.obj <- subset(obj, donor_id == "D165")
d165.metadf <- d165.obj@meta.data

anatomic_sites <- unique(d165.metadf$anatomic_site)

d165.postaur.obj <- subset(obj, subset = donor_id == "D165" & anatomic_site == "Postauricular")
d165.occip.obj <- subset(obj, subset = donor_id == "D165" & anatomic_site == "Occipital Scalp")
spatial.feature_plot.genes <-  c("COL17A1", "FLG", "KLF5", "SOX9","KRT79","COL1A1","VWF", "ADIPOQ", "MRC1")

global.dotplot.marker_genes <- list(
    "Basal KC" = c("COL17A1", "COL7A1"),
    "Spn KC" = c("RHOV", "KLF5"),
    "Grn KC" = c("DSC1", "FLG"),
    "Seb" = c("KRT79", "TKT"),
    "Ecc Duct" = c("KRT77", "SEMA3C"),
    "Ecc Gland" = c("SFRP1", "DNER"),
    "HFE" = c("GJB6", "SOX9"),
    "Fibro" = c("COL1A1","PDGFRA"),
    "Peri" = c("ACTA2", "RGS5"),
    "EC" = c("VWF", "ACKR1"),
    "LEC" = c("CCL21", "LYVE1"),
    "Adipo" = c("ADIPOQ", "G0S2"),
    "Schwann" = c("SOX10", "MPZ"),
    "Lym"=c("CD3E", "CD4", "CD8A"),
    "Plasma"=c("JCHAIN"),
    "LC" = c("CD207", "CD1A"),
    "Mac/DC" = c('CD83', "MRC1", "C1QA"),
    "Mast" = c("IL1RL1", "KIT")
  )

global.dotplot.marker_genes <- unique(unlist(global.dotplot.marker_genes))


## PLOTS
## 1. Global / Broad Cell Type UMAP Plot
global.umap <- DimPlot(obj, reduction = 'scvi_umap', alpha = 0.25, group.by = 'cell_type.broad', label = TRUE, raster.dpi = c(1200, 1200)) +
  spatial_theme + NoAxes() +
  scale_color_manual(values = global.cell_type.colors) +
  coord_fixed() + labs(title = NULL) +
    guides(color = guide_legend(override.aes = list(size = 5, alpha = 1))) +
    theme(legend.text = element_text(size = 12, color='black'))

global.umap

svglite::svglite("figures/merfish/BAYSOR/annotated/cell_type.broad.scvi_umap.svg", height = 6, width = 8)
print(global.umap)
dev.off()

global.umap <- DimPlot(obj, reduction = 'scvi_umap', alpha = 0.25, group.by = 'cell_type.broad', label = FALSE,  raster.dpi = c(1200, 1200)) +
  spatial_theme + NoAxes() +
  scale_color_manual(values = global.cell_type.colors) +
  coord_fixed() + labs(title = NULL) +
  guides(color = guide_legend(override.aes = list(size = 5, alpha = 1))) +
  theme(legend.text = element_text(size = 12, color='black')) + NoLegend()
png("figures/merfish/BAYSOR/annotated/cell_type.broad.scvi_umap.png", height = 6, width = 8, units = 'in', res = 300)
print(global.umap)
dev.off()

metadf$collection_type <- factor(ifelse(as.character(metadf$collection_source) =='autopsy', 'autopsy', 'surgical discards'))

collection_plot1 <- ggplot(metadf) +
  aes(x = cell_type.broad, fill = collection_type) +
  geom_bar(position='fill', color='black') +
  theme_classic() +
  labs(x = "Cell Type", y = "Proportion (%)") +
  scale_fill_manual(values = c("dodgerblue", 'goldenrod1'), name = "Collection Strategy") +
  scale_y_continuous(labels = scales::percent) +
  rotate_x_text(angle=90) +
  theme(axis.text.x=element_text(color='black', size=12),
        axis.text.y=element_text(color='black', size=12),
        axis.title.x=element_text(color='black', size=12, face='bold'),
        axis.title.y = element_text(color='black', size=12, face='bold'),
        legend.text = element_text(color='black', size=12),
        legend.title=element_text(color='black', size=12, face='bold'))

svglite::svglite("figures/merfish/BAYSOR/annotated/broad.cell_type.proportions_by_collection_strategy.svg", height=5,width=7)
collection_plot1
dev.off()

collection_plot2 <- ggplot(metadf) +
  aes(x = collection_source, fill = cell_type.broad) +
  geom_bar(position='fill', color='black') +
  theme_classic() + theme + labs(x = "Collection Strategy") +
  scale_fill_manual(values = global.cell_type.colors) +
  scale_y_continuous(labels = scales::percent)

collection_plot2

compartment_plot <- ggplot(metadf) +
  aes(x = spatial_1, y = spatial_2, color = tissue_compartment) +
  geom_point(size = 0.25, alpha=0.5) +
  facet_wrap(~ sample_barcode, scales = 'free', nrow=9) +
  scale_color_manual(values = c(epidermis="orchid1", dermis='goldenrod1', subcutis='dodgerblue')) +
  spatial_theme + NoLegend() + NoAxes()

png("figures/merfish/BAYSOR/annotated/compartment_plot.png", width=48, height=36, res=150, units='in')
print(compartment_plot)
dev.off()



## 2. Cell Type Spatial - D165
anatomic_site.spatial.global_celltype.plot_list <- lapply(anatomic_sites, function(x){

  df <- d165.metadf[d165.metadf$anatomic_site == x,]

  min1 <- min(df$spatial_1)
  min2 <- min(df$spatial_2)
  diff <- max((max(df$spatial_1)-min(df$spatial_1)), max(df$spatial_1)-min(df$spatial_1))

  plt <- ggplot(df) +
    aes(x = spatial_1, y = spatial_2, color = cell_type.broad) +
    geom_point(size = 0.25) +
    theme_classic() + spatial_theme + coord_fixed() +
    scale_color_manual(values = global.cell_type.colors) +
    annotate('segment', x = min1-100, y = min2-100, xend = min1 + 1000 - 100, yend = min2-100,
             color = 'black', linewidth = 1.5 , lineend='square') +
    theme(axis.text = element_blank(), axis.ticks = element_blank()) +
    NoLegend() + NoAxes() + labs(title = x)

  return(plt)
})

d165.cell_type.spatial <- cowplot::plot_grid(plotlist = anatomic_site.spatial.global_celltype.plot_list, nrow = 2)

png(paste0("figures/merfish/BAYSOR/annotated/D165.cell_type.spatial.fixed.png"), width = 36, height = 18, units = "in", res = 600)
print(d165.cell_type.spatial)
dev.off()

svglite::svglite("figures/merfish/BAYSOR/annotated/D165.cell_type.spatial.fixed.svg", height = 16, width = 36, bg = 'transparent')
print(d165.cell_type.spatial)
dev.off()

## global spatial feature plots -- use occipital scalp
df <- d165.occip.obj@meta.data
min1 <- min(df$spatial_1)
min2 <- min(df$spatial_2)
diff <- max((max(df$spatial_1)-min(df$spatial_1)), max(df$spatial_1)-min(df$spatial_1))

global.spatial.feature_plots <- FeaturePlot(d165.occip.obj, features = spatial.feature_plot.genes, reduction = 'spatial',
                   ncol = length(spatial.feature_plot.genes), cols = c("gray0", "#3f0fff"), raster = FALSE, order = TRUE) &
  spatial_theme & theme & NoAxes() & coord_fixed() & NoLegend()

global.spatial.feature_plots <- global.spatial.feature_plots +
  annotate('segment', x = min1, y = min2, xend = min1 + 1000, yend = min2, color = 'black', linewidth = 1.5 , lineend='square')

svglite::svglite("figures/merfish/BAYSOR/annotated/D165.feature.spatial.fixed.svg", height = 6, width = 36)
print(global.spatial.feature_plots)
dev.off()


png(paste0("figures/merfish/BAYSOR/annotated/D165.feature.spatial.png"), height = 6, width = 36, units = "in", res = 600)
print(global.spatial.feature_plots)
dev.off()


global.marker.dotplot <- DotPlot(obj, group.by = "cell_type.broad", features = global.dotplot.marker_genes,
                   cols = global.colors, col.min = 0, dot.min = 0.05) +
  labs(title = NULL, x = NULL, y = NULL) + theme +
  scale_size_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100)) +
  guides(color = guide_colorbar(title = "Average\nExpression", frame.colour = 'black'),
         size = guide_legend(title = "Percent\nExpressed")) +
  theme(panel.grid.minor = element_line(color = "black", size = 0.5),
        panel.margin.x=unit(0.25, "lines"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12, color = 'black'),
        axis.text.y = element_text(hjust = 1, vjust = 0.6, size = 12, color = 'black')
        )
# global.marker.dotplot
svglite::svglite("figures/merfish/BAYSOR/annotated/cell_type.marker_dotplot.svg", width = 18, height=6)
global.marker.dotplot
dev.off()


### neighborhoods -- use postauricular
df <- d165.postaur.obj@meta.data
min1 <- min(df$spatial_1)
min2 <- min(df$spatial_2)
diff <- max((max(df$spatial_1)-min(df$spatial_1)), max(df$spatial_1)-min(df$spatial_1))

## full neighborhood spatial plot
neighborhood.spatial_plot <- ggplot(df) +
  aes(x = spatial_1, y = spatial_2, color = neighborhood.detailed) +
  geom_point() +
  theme_classic() +  spatial_theme + coord_fixed() +
  scale_color_manual(values = neighborhood.colors) +
  annotate('segment', x = min1-100, y = min2-100, xend = min1 + 1000 - 100, yend = min2-100,
           color = 'black', linewidth = 1.5 , lineend='square') +
  theme(axis.text = element_blank(), axis.ticks = element_blank()) +
  NoLegend() + NoAxes()

svg("figures/merfish/BAYSOR/annotated/neighborhood.detailed.D165_postauricular.spatial.svg", height=8, width=8)
neighborhood.spatial_plot
dev.off()


## panel with spatial highlights for each neighborhood
neighborhood.highlight.pltlist <- lapply(neighborhoods, function(cluster){
  cells <- colnames(d165.postaur.obj)[d165.postaur.obj$neighborhood.detailed == cluster]
  DimPlot(d165.postaur.obj, group.by = "neighborhood.detailed", cols.highlight = c(neighborhood.colors[cluster]), cells.highlight = cells,
          reduction = 'spatial', pt.size = 0.5, sizes.highlight = 0.5) +
    spatial_theme + theme + NoAxes() + coord_fixed() + NoLegend() + labs(title = cluster) +
    annotate('segment', x = min1-100, y = min2-100, xend = min1 + 1000 - 100, yend = min2-100,
             color = 'black', linewidth = 1.5 , lineend='square')
})

neighborhood.highlight.spatial_panel <- cowplot::plot_grid(plotlist = neighborhood.highlight.pltlist, nrow=2)
file_prefix <- paste0("figures/merfish/BAYSOR/annotated/D165.neighborhood_highlight.postauricular_final")

svglite::svglite(paste0(file_prefix, ".svg"), height=12, width=36)
print(neighborhood.highlight.spatial_panel)
dev.off()

donors <- c("D107", "D151", "D145", "D149", "D165", "D077", "D082")
for (donor in donors){
  donor.obj <- subset(obj, donor_id == donor)
  sites <- unique(as.character(donor.obj$anatomic_site))

  site.list <- lapply(sites, function(x){
    sample.obj <- subset(donor.obj, anatomic_site == x)
    df <- sample.obj@meta.data
    min1 <- min(df$spatial_1)
    min2 <- min(df$spatial_2)
    diff <- max((max(df$spatial_1)-min(df$spatial_1)), max(df$spatial_1)-min(df$spatial_1))
    ## full neighborhood spatial plot
    res <- ggplot(df) +
      aes(x = spatial_1, y = spatial_2, color = neighborhood.detailed) +
      geom_point(size=0.25) +
      theme_classic() +  spatial_theme + coord_fixed() +
      scale_color_manual(values = neighborhood.colors) +
      annotate('segment', x = min1-100, y = min2-100, xend = min1 + 1000 - 100, yend = min2-100,
               color = 'black', linewidth = 1.5 , lineend='square') +
      theme(axis.text = element_blank(), axis.ticks = element_blank()) +
      NoLegend() + NoAxes() +labs(title=x)
    return(res)
  })

  donor.panel <- cowplot::plot_grid(plotlist = site.list, nrow=2)

  png(paste0("figures/merfish/BAYSOR/annotated/neighborhood_spatial.", donor, ".png"), height=12, width=24, units = 'in', res = 300)
  print(donor.panel)
  dev.off()

  svglite::svglite(paste0("figures/merfish/BAYSOR/annotated/neighborhood_spatial.", donor, ".svg"), height=12, width=24)
  print(donor.panel)
  dev.off()
}



donors <- c("D107", "D151", "D145", "D149", "D165", "D077", "D082")
for (donor in donors){
  donor.obj <- subset(obj, donor_id == donor)
  sites <- unique(as.character(donor.obj$anatomic_site))

  site.list <- lapply(sites, function(x){
    sample.obj <- subset(donor.obj, anatomic_site == x)
    df <- sample.obj@meta.data
    ## full neighborhood spatial plot
    res <- ggplot(df) +
      aes(x = spatial_1, y = spatial_2, color = neighborhood.detailed) +
      geom_point(size=0.25) +
      theme_classic() +  spatial_theme + coord_fixed() +
      scale_color_manual(values = neighborhood.colors) +
      theme(axis.text = element_blank(), axis.ticks = element_blank()) +
      NoLegend() + NoAxes() +labs(title=x)
    return(res)
  })

  donor.panel <- cowplot::plot_grid(plotlist = site.list, nrow=2)

  png(paste0("figures/merfish/BAYSOR/annotated/neighborhood_spatial.", donor, ".noscalebar.png"), height=12, width=24, units = 'in', res = 300)
  print(donor.panel)
  dev.off()

  svglite::svglite(paste0("figures/merfish/BAYSOR/annotated/neighborhood_spatial.", donor, ".noscalebar.svg"), height=12, width=24)
  print(donor.panel)
  dev.off()
}



# spatial plots for each donor
donors <- c("D107", "D151", "D145", "D149", "D165", "D077", "D082")
for (donor in donors){
  donor.obj <- subset(obj, donor_id == donor)
  sites <- unique(as.character(donor.obj$anatomic_site))

  site.list <- lapply(sites, function(x){
    sample.obj <- subset(donor.obj, anatomic_site == x)
    df <- sample.obj@meta.data
    min1 <- min(df$spatial_1)
    min2 <- min(df$spatial_2)
    diff <- max((max(df$spatial_1)-min(df$spatial_1)), max(df$spatial_1)-min(df$spatial_1))


    ## panel with spatial highlights for each neighborhood
    res <- lapply(neighborhoods, function(cluster){
      cells <- colnames(sample.obj)[sample.obj$neighborhood.detailed == cluster]
      DimPlot(sample.obj, group.by = "neighborhood.detailed",
              cols = "gray90",
              cols.highlight = c(neighborhood.colors[cluster]),
              cells.highlight = cells, order = cluster,
              reduction = 'spatial', pt.size = 0.5, sizes.highlight = 0.5) +
        spatial_theme + NoAxes() + coord_fixed() + NoLegend() + labs(title = cluster) +
        annotate('segment', x = min1-100, y = min2-100, xend = min1 + 1000 - 100, yend = min2-100,
                 color = 'black', linewidth = 1.5 , lineend='square')
    })

    res <- cowplot::plot_grid(plotlist = res, nrow=1)
    return(res)
  })

  donor.panel <- cowplot::plot_grid(plotlist = site.list, ncol=1)

  png(paste0("figures/merfish/BAYSOR/annotated/neighborhood.spatial_highlight.", donor, ".png"), height=36, width=24, units = 'in', res = 300)
  print(donor.panel)
  dev.off()

  donor.panel <- cowplot::plot_grid(plotlist = site.list, ncol=1, labels = sites, vjust = 0.5, hjust = 0, label_x = 0, label_y = 0.5)

  svglite::svglite(paste0("figures/merfish/BAYSOR/annotated/neighborhood.spatial_highlight.", donor, ".svg"), height=36, width=24)
  print(donor.panel)
  dev.off()
}

library(tidyverse)
library(DirichletReg)

# make a table of samples by clusters/cell-types -> number of cells per cluster per sample
cell_counts <- data.frame(table(meta$cell_type.detailed, meta$neighborhood))
colnames(cell_counts) <- c("cell_type.detailed", "neighborhood", "cell_count")

x <- table(meta$cell_type.detailed, meta$neighborhood)
class(x) <- 'matrix'

# makes cell type proportions
drmatrix <- DR_data(x)

detailed_celltype.colors <- readRDS("data/reference/cell_type.detailed.color_palette.varibow.non_shuffled.RDS")
detailed.cell_type.colors <- detailed_celltype.colors
###
cell_prop <- as.matrix(drmatrix)
class(cell_prop) <- 'matrix'
###
library(ComplexHeatmap)
row_ha = rowAnnotation(" " = names(neighborhood.colors),
                          col = list(" "=neighborhood.colors))
## make the top annotation
colors <- detailed_celltype.colors

column_ha <- HeatmapAnnotation(
  " " = names(detailed_celltype.colors[rownames(cell_prop)]),
  col = list(" " = colors[rownames(cell_prop)]),
  na_col = "grey"
)

svglite::svglite("figures/merfish/BAYSOR/annotated/neighborhood_composition_proportion_heatmap.svg", height = 5, width = 18)
Heatmap(t(cell_prop), name='prop',
        col = circlize::colorRamp2(c(0,1), c('white', 'blue')),
        row_names_side = 'left', row_dend_side = 'right',
        bottom_annotation = column_ha,
        left_annotation = row_ha,
        cluster_rows = TRUE, cluster_columns = TRUE)
dev.off()

## Neighborhood Composition Donut Plots
donut_list <- list()
for (nhood in neighborhoods){
  df <- meta %>%
    filter(neighborhood.detailed == nhood) %>%
    group_by(neighborhood.detailed, cell_type.broad, cell_type.detailed) %>%
    summarise(no.cells = n()) %>%
    ungroup()

  sum_total_cells <- sum(df$no.cells)
  firstLevel <- df %>% group_by(neighborhood.detailed) %>% summarize(total.cells=sum(no.cells))

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
    geom_bar(data = secondLevel, aes(x = 2, y = pct.cells, fill = cell_type.detailed),
             position = 'stack', stat = 'identity', color = 'black') +
    geom_text(data = secondLevel[secondLevel$pct.cells>0.01,],
              aes(label = cell_type.detailed, x = 2, y = pos, angle = angle),
              position = 'identity', color = 'black', fontface = 'bold') +
    scale_fill_manual(values = detailed.cell_type.colors) +
    coord_polar('y') + theme_minimal() + spatial_theme + NoAxes() + labs(title = nhood) + NoLegend()
}

neighborhood_composition.donut_panel <- patchwork::wrap_plots(donut_list, nrow=2, guides = 'collect')
neighborhood_composition.donut_panel

svglite::svglite("figures/merfish/BAYSOR/annotated/donut_panel.cell_type.detailed.svg", height = 14, width = 36)
neighborhood_composition.donut_panel
dev.off()

## Neighborhood Spatial - ALL SAMPLES
sample_barcodes <- sort(unique(meta$sample_barcode))
all.neighborhood.spatial.list <- lapply(sample_barcodes, function(x){

  df <- meta[meta$sample_barcode == x,]

  min1 <- min(df$spatial_1)
  min2 <- min(df$spatial_2)
  diff <- max((max(df$spatial_1)-min(df$spatial_1)), max(df$spatial_1)-min(df$spatial_1))

  plt <- ggplot(df) +
    aes(x = spatial_1, y = spatial_2, color = neighborhood.detailed) +
    geom_point(size = 0.1) +
    theme_classic() + spatial_theme + coord_fixed() +
    scale_color_manual(values = neighborhood.colors) +
    annotate('segment', x = min1-100, y = min2-100, xend = min1 + 1000 - 100, yend = min2-100,
             color = 'black', linewidth = 1.5 , lineend='square') +
    theme(axis.text = element_blank(), axis.ticks = element_blank()) +
    NoLegend() + NoAxes(keep.ticks = FALSE, keep.text = FALSE) + labs(title = NULL, x=NULL, y=NULL)

  return(plt)
})
names(all.neighborhood.spatial.list) <- sample_barcodes


all.neighborhood.spatial_panel <- cowplot::plot_grid(plotlist = all.neighborhood.spatial.list, nrow = 4)

png(paste0("figures/merfish/BAYSOR/annotated/neighborhood_spatial.all_samples.png"), width = 46, height = 10, units = "in", res = 600, bg = 'transparent')
print(all.neighborhood.spatial_panel)
dev.off()

png(paste0("figures/merfish/BAYSOR/annotated/neighborhood_spatial.all_samples.skinny.png"), width = 46, height = 10, units = "in", res = 600, bg = 'transparent')
print(all.neighborhood.spatial_panel)
dev.off()

all.neighborhood.spatial_panel <- cowplot::plot_grid(plotlist = all.neighborhood.spatial.list, nrow = 8)
png(paste0("figures/merfish/BAYSOR/annotated/neighborhood_spatial.all_samples.tall.png"), width = 46, height = 14, units = "in", res = 600, bg = 'transparent')
print(all.neighborhood.spatial_panel)
dev.off()


## make fibroblast spatial plots
fibroblasts <- detailed.cell_types[grepl("Fib", detailed.cell_types)]
donors <- c("D107", "D151", "D145", "D149", "D165", "D077", "D082")
for (donor in donors){
  donor.obj <- subset(obj, donor_id == donor)
  sites <- unique(as.character(donor.obj$anatomic_site))

  site.list <- lapply(sites, function(x){
    sample.obj <- subset(donor.obj, anatomic_site == x)
    df <- sample.obj@meta.data
    min1 <- min(df$spatial_1)
    min2 <- min(df$spatial_2)
    diff <- max((max(df$spatial_1)-min(df$spatial_1)), max(df$spatial_1)-min(df$spatial_1))


    ## panel with spatial highlights for each neighborhood
    res <- lapply(fibroblasts, function(cluster){
      cells <- colnames(sample.obj)[sample.obj$cell_type.detailed == cluster]
      DimPlot(sample.obj, group.by = "cell_type.detailed",
              cols = "gray90",
              cols.highlight = c(detailed.cell_type.colors[cluster]),
              cells.highlight = cells, order = cluster,
              reduction = 'spatial', pt.size = 0.5, sizes.highlight = 0.5) +
        spatial_theme + NoAxes() + coord_fixed() + NoLegend() + labs(title = cluster) +
        annotate('segment', x = min1-100, y = min2-100, xend = min1 + 1000 - 100, yend = min2-100,
                 color = 'black', linewidth = 1.5 , lineend='square')
    })

    res <- cowplot::plot_grid(plotlist = res, nrow=1)
    return(res)
  })

  donor.fib.panel <- cowplot::plot_grid(plotlist = site.list, ncol=1, labels = sites, vjust = 0.5, hjust = 0, label_x = 0, label_y = 0.5)

  png(paste0("figures/merfish/BAYSOR/annotated/fibroblast_spatial_highlight.", donor, ".png"), height=36, width=24, units = 'in', res = 300)
  print(donor.fib.panel)
  dev.off()

  svglite::svglite(paste0("figures/merfish/BAYSOR/annotated/fibroblast_spatial_highlight.", donor, ".svg"), height=36, width=24)
  print(donor.fib.panel)
  dev.off()
}

## plot k stability
kstabdf <- read.csv("data/merfish/BAYSOR/cellcharter/ns-atlas.merfish_baysor.scanvi_integrated.cellcharter.k_stability_matrix.csv")
iters <- colnames(kstabdf)[2:ncol(kstabdf)]
kstabdf <- kstabdf %>% pivot_longer(cols = all_of(iters), names_to = 'iter', values_to = 'stability')
colnames(kstabdf) <- c("k.clusters", "iter", "stability")
stabmeans <- kstabdf %>% group_by(k.clusters) %>%
  summarise(standard_error = sd(stability) / sqrt(length(stability)),
    stability = mean(stability))

kvalues <- seq(min(stabmeans$k.clusters), max(stabmeans$k.clusters), by = 1)

kstab.plot <- ggplot(stabmeans) + aes(x = k.clusters, y = stability) +
  geom_ribbon(aes(ymin=stability - standard_error, ymax=stability + standard_error), alpha=0.25, fill = "gray50") +
  geom_line() + geom_point() + theme_bw() + theme + labs(x = "K Clusters", y = "Cluster Stability", title = "CellCharter K Stability") +
  scale_x_continuous(labels = kvalues, breaks = kvalues,
                     limits = c(min(stabmeans$k.clusters), max(stabmeans$k.clusters))) +
  theme(axis.text = element_text(color='black', size=12), axis.title=element_text(color='black', size=14))

svglite::svglite("figures/merfish/BAYSOR/annotated/merfish.cellcharter.kstability_plot.svg", height=3, width=7)
kstab.plot
dev.off()

kmeans.elbow <- read.csv("data/merfish/BAYSOR/neighborhood_analysis/cellcharter_mbkmeans_qc_metrics.csv")

ggplot(kmeans.elbow) + aes(x=K, y = silhouette) + geom_point() + geom_line()

wcss.elbow <- ggplot(kmeans.elbow) + aes(x=K, y = ssd) + geom_point() + geom_line() +
  theme_bw() + theme + labs(x = "K Clusters", y = "Within-Cluster Sum of Squared Distances", title = "CellCharter MiniBatchKMeans Elbow")



gmm.elbow <- read.csv("data/merfish/BAYSOR/neighborhood_analysis/cellcharter_gmm_qc_metrics.csv")

aic.elbow <- ggplot(gmm.elbow) + aes(x=K, y = aic) + geom_point() + geom_line() +
  theme_bw() + theme + labs(x = "K Clusters", y = "AIC", title = "CellCharter GMM Elbow")

bic.elbow <- ggplot(gmm.elbow) + aes(x=K, y = bic) + geom_point() + geom_line()+
  theme_bw() + theme + labs(x = "K Clusters", y = "BIC", title = "CellCharter GMM Elbow")

wcss.elbow | aic.elbow | bic.elbow

ggplot(kmeans.elbow) + aes(x=K, y = ssd) + geom_point() + geom_line() +
  theme_bw() + theme + labs(x = "K Clusters", y = "SSD", title = "CellCharter MiniBatchKMeans Elbow")


nhood.compartment.plot <- ggplot(meta) +
  aes(x = neighborhood.detailed, fill = tissue_compartment) +
  geom_bar(position='fill', color='black') +
  theme_classic() +
  labs(x = NULL, y = "Proportion (%)") +
  scale_fill_manual(values = c(epidermis="orchid1", dermis='goldenrod1', subcutis='dodgerblue'), name = "Tissue Compartment") +
  scale_y_continuous(labels = scales::percent) + coord_flip() +
  theme(axis.text.x=element_text(color='black', size=12),
        axis.text.y=element_text(color='black', size=12),
        axis.title.x=element_text(color='black', size=12, face='bold'),
        axis.title.y = element_text(color='black', size=12, face='bold'),
        legend.text = element_text(color='black', size=12),
        legend.title=element_text(color='black', size=12, face='bold'))

svglite::svglite("figures/annotated/neighborhood.donor_id.tissue_compartment.proportion_plot.svg", height=4, width=12)
nhood.compartment.plot
dev.off()

nhood.collection.plot <- ggplot(metadf) +
  aes(x = neighborhood.detailed, fill = collection_type) +
  geom_bar(position='fill', color='black') +
  theme_classic() +
  labs(x = NULL, y = "Proportion (%)") +
  scale_fill_manual(values = c("dodgerblue", 'goldenrod1'), name = "Collection Strategy") +
  scale_y_continuous(labels = scales::percent) + coord_flip() +
  theme(axis.text.x=element_text(color='black', size=12),
        axis.text.y=element_text(color='black', size=12),
        axis.title.x=element_text(color='black', size=12, face='bold'),
        axis.title.y = element_text(color='black', size=12, face='bold'),
        legend.text = element_text(color='black', size=12),
        legend.title=element_text(color='black', size=12, face='bold'))
nhood.collection.plot

nhood.donor.plot <- ggplot(metadf) +
  aes(x = neighborhood.detailed, fill = donor_id) +
  geom_bar(position='fill', color='black') +
  theme_classic() +
  labs(x = NULL, y = "Proportion (%)") +
  scale_y_continuous(labels = scales::percent) + coord_flip() +
  theme(axis.text.x=element_text(color='black', size=12),
        axis.text.y=element_text(color='black', size=12),
        axis.title.x=element_text(color='black', size=12, face='bold'),
        axis.title.y = element_text(color='black', size=12, face='bold'),
        legend.text = element_text(color='black', size=12),
        legend.title=element_text(color='black', size=12, face='bold'))

svglite::svglite("figures/merfish/BAYSOR/annotated/neighborhood.donor_id.collection_strategy.proportion_plot.svg", height=6, width=18)
nhood.collection.plot | nhood.donor.plot
dev.off()


cell_type.propdf.dermis <- readRDS("data/merfish/differential_abundance/cell_type.detailed.anatomic_site.cell_proportions.dermis.RDS")
cell_type.propdf <- cell_type.propdf.dermis %>% filter(anatomic_site %in% sites & cell_type.detailed %in% fibroblasts)
cell_type.wilcoxdf <- readRDS("data/merfish/differential_abundance/cell_type.detailed.no_doublets.proportions.wilcox_test.one_vs_all.dermis.RDS")
pvals <- cell_type.wilcoxdf %>% filter(p.adj < 0.05 & cell_type.detailed %in% fibroblasts)

fib.abundance <- ggplot(data = cell_type.propdf)  +
  aes(x = anatomic_site, y = cell_proportion, fill = anatomic_site) +
  geom_point(size = -1, aes(fill = anatomic_site)) +
  facet_wrap(~ cell_type.detailed, scales = 'free', nrow=1) +
  geom_boxplot(alpha = 0.7, show.legend=FALSE) +
  geom_text(data = pvals,
            aes(x = anatomic_site, y = Inf,
                label = padj.signif,
                color = direction),
            hjust = 0.5, vjust = 1, size = 8, fontface = 'bold',
            show.legend = FALSE) +
  scale_color_manual(values = c("+" = "tomato", "-" = "dodgerblue") ) +
  scale_fill_paletteer_d("Polychrome::light") +
  theme_bw() + theme + labs(x = NULL, y = "Cell Proportion per Sample") +
  theme(axis.text.x = element_text(size = 14, color='black'),
        legend.text = element_text(color='black', size = 14, hjust = 0)) +
  guides(fill = guide_legend(override.aes = list(shape = 22, size=10), title = NULL)) +
  NoLegend() + rotate_x_text(angle=90)

svglite::svglite("figures/merfish/differential_abundance/dermis_fib_abundance.boxplot.svg", height=5, width=18)
print(fib.abundance)
dev.off()



perivasc.fibs <- subset(obj, cell_type.detailed == "Perivasc Fib")


VlnPlot(perivasc.fibs, features = "CCL19", group.by = 'anatomic_site', log = TRUE)



## retic fib III and IV + schwann spatial examples


donors <- c("D077", "D145", "D149", "D165")
sites <- c("sole", "postauricular")

highlight.colors <- c('Retic Fib II'="seagreen1", "Retic Fib III"="dodgerblue", "Schwann"="purple", "Other"='gray90')
obj@meta.data$highlight <- factor(ifelse(as.character(obj@meta.data$cell_type.detailed) %in% names(highlight.colors),as.character(obj@meta.data$cell_type.detailed), "Other"), levels = names(highlight.colors))
table(obj@meta.data$highlight)




plot.list <- lapply(donors, function(donor){
  plt <- lapply(sites, function(site){
    sample.obj <- subset(obj, donor_id == donor & anatomic_site == site)
    df <- sample.obj@meta.data
    ggplot(df) +
      aes(x = center_x, y = center_y, color = highlight) +
      geom_point(size=1, alpha=0.5) + theme_minimal() + spatial_theme +
      NoAxes() + scale_color_manual(values = highlight.colors) +
      coord_fixed() +  NoLegend() # + labs(title = paste0(donor, site))
  })
  cowplot::plot_grid(plotlist = plt, ncol=1)
})

png("figures/merfish/BAYSOR/annotated/fib34panel.nolab.png",width = 18,height = 8, units='in', res=300)
cowplot::plot_grid(plotlist = plot.list, ncol=4)
dev.off()
## scalp vs sole bas and spn kc vs inf spatial



## Cell Type Detailed - ALL SAMPLES
sample_barcodes <- sort(unique(metadf$sample_barcode))
all.celltypes.spatial.list <- lapply(sample_barcodes, function(x){

  df <- metadf[metadf$sample_barcode == x,]

  min1 <- min(df$spatial_1)
  min2 <- min(df$spatial_2)
  diff <- max((max(df$spatial_1)-min(df$spatial_1)), max(df$spatial_1)-min(df$spatial_1))

  plt <- ggplot(df) +
    aes(x = spatial_1, y = spatial_2, color = cell_type.detailed) +
    geom_jitter(size = 0.25) +
    theme_classic() + spatial_theme + coord_fixed() +
    scale_color_manual(values = detailed.cell_type.colors) +
    annotate('segment', x = min1-100, y = min2-100, xend = min1 + 1000 - 100, yend = min2-100,
             color = 'black', linewidth = 1.5 , lineend='square') +
    theme(axis.text = element_blank(), axis.ticks = element_blank()) +
    NoLegend() + NoAxes(keep.ticks = FALSE, keep.text = FALSE) + labs(title = NULL, x=NULL, y=NULL)

  return(plt)
})
names(all.celltypes.spatial.list) <- sample_barcodes


all.ctd.spatial_panel <- cowplot::plot_grid(plotlist = all.celltypes.spatial.list, ncol = 11, labels = sample_barcodes)
png(paste0("figures/merfish/BAYSOR/annotated/celltypedetailed_spatial.all_samples.all.png"), width = 28, height = 36, units = "in", res = 600, bg = 'transparent')
print(all.ctd.spatial_panel)
dev.off()

sampledf <- unique(metadf[,c("sample_barcode", "donor_id", "collection_source", "anatomic_site")])
for (sample in sample_barcodes){
  site = sampledf[sampledf$sample_barcode == sample, "anatomic_site"]
  collection = sampledf[sampledf$sample_barcode == sample, "collection_source"]
  fname = paste0("figures/merfish/BAYSOR/annotated/sample_spatial/", make.names(collection), ".", sample, ".", site, ".detailed_celltype.spatial_plot.png")

  png(fname, width = 4, height = 4, units = "in", res = 300, bg = 'transparent')
  print(all.celltypes.spatial.list[[sample]])
  dev.off()
}

detailed.cell_types <- names(detailed.cell_type.colors)
sample_barcodes <- sort(as.character(sample_barcodes))
folder <- 'figures/merfish/BAYSOR/annotated/detailed_celltype.spatial_highlight/'

### spatial cluster highlight
for (cell_type in detailed.cell_types){
  print(paste0("Plotting Cell Type: ", cell_type, "!"))

  metadf$highlight <- ifelse(as.character(metadf$cell_type.detailed) == cell_type, cell_type, "Other")
  metadf$highlight <- factor(metadf$highlight, levels = c( "Other", cell_type), ordered = TRUE)
  metadf <- metadf %>% arrange(sample_barcode, highlight, center_x, center_y)
  cols <- c(detailed.cell_type.colors[cell_type], 'gray90')
  names(cols) <- c(cell_type, "Other")

  spatial.list <- lapply(sample_barcodes, function(x){

    df <- metadf[metadf$sample_barcode == x,]

    min1 <- min(df$spatial_1)
    min2 <- min(df$spatial_2)
    diff <- max((max(df$spatial_1)-min(df$spatial_1)), max(df$spatial_1)-min(df$spatial_1))

    plt <- ggplot(df) +
      aes(x = spatial_1, y = spatial_2, color = highlight) +
      geom_jitter(size = 0.25, alpha=0.7) +
      spatial_theme + coord_fixed() +
      scale_color_manual(values = cols) +
      annotate('segment',
               x = min1-100, xend = min1 + 1000 - 100,
               y = min2-100,  yend = min2-100,
               color = 'black', linewidth = 1.5 ,
               lineend='square') +
      theme(axis.text = element_blank(),
            axis.ticks = element_blank()) +
      NoLegend() +
      NoAxes(keep.ticks = FALSE, keep.text = FALSE) +
      labs(title = x, x = NULL, y = NULL)

    return(plt)
  })
  names(spatial.list) <- sample_barcodes


  spatial_panel <- patchwork::wrap_plots(plotlist = spatial.list, ncol = 9) +
    plot_annotation(title = cell_type)



  file_prefix <- paste0(folder,"/detailed.cell_type.", make.names(cell_type), ".spatial_highlight.all_samples.")

  png(filename = paste0(file_prefix, "png"), height = 40, width=30, units = 'in', res=300)
  print(spatial_panel)
  dev.off()

  file_prefix <- paste0(folder,"/detailed.cell_type.", make.names(cell_type), ".spatial_highlight.all_samples.nolabels.")

  png(filename = paste0(file_prefix, "png"), height = 40, width=30, units = 'in', res=300)
  print(spatial_panel & labs(title=NULL))
  dev.off()

  print(paste0(file_prefix, "png DONE!"))

}

x = "D042_SKIN_NS_S01_R01"
samples2plot <- c("D042_SKIN_NS_S01_R01", "D165_SKIN_NS_S01_R01", "D165_SKIN_NS_S03_R01", "D151_SKIN_NS_S02_R01", "D151_SKIN_NS_S12_R01", "D165_SKIN_NS_S12_R01")
metadf$spatial_1 <- metadf$center_x
metadf$spatial_2 <- metadf$center_y
for (x in samples2plot){
  df <- metadf[metadf$sample_barcode == x,]

  min1 <- min(df$spatial_1)
  min2 <- min(df$spatial_2)
  diff <- max((max(df$spatial_1)-min(df$spatial_1)), max(df$spatial_1)-min(df$spatial_1))

  detailed.cell_types <- levels(metadf$cell_type.detailed)

  spatial.list <- lapply(detailed.cell_types, function(cell_type){


    df$highlight <- ifelse(as.character(df$cell_type.detailed) == cell_type, cell_type, "Other")
    df$highlight <- factor(df$highlight, levels = c( "Other", cell_type), ordered = TRUE)
    df <- df %>% arrange(highlight, center_x, center_y)
    cols <- c(detailed.cell_type.colors[cell_type], 'gray90')
    names(cols) <- c(cell_type, "Other")

    plt <- ggplot(df) +
      aes(x = spatial_1, y = spatial_2, color = highlight) +
      geom_jitter(size = 0.25, alpha=0.7) +
      spatial_theme + coord_fixed() +
      scale_color_manual(values = cols) +
      annotate('segment',
               x = min1-100, xend = min1 + 1000 - 100,
               y = min2-100,  yend = min2-100,
               color = 'black', linewidth = 1.5 ,
               lineend='square') +
      theme(axis.text = element_blank(),
            axis.ticks = element_blank()) +
      NoLegend() +
      NoAxes(keep.ticks = FALSE, keep.text = FALSE) +
      labs(title = cell_type, x = NULL, y = NULL)

    return(plt)
  })
  names(spatial.list) <- detailed.cell_types
  spatial_panel <- cowplot::plot_grid(plotlist=spatial.list, greedy = TRUE, nrow=4)
  file_prefix <- paste0(folder,"/detailed_celltype.sample_highlights.", make.names(x), ".")

  png(filename = paste0(file_prefix, "png"), height = 20, width=40, units = 'in', res=300)
  print(spatial_panel)
  dev.off()

  file_prefix <- paste0(folder,"/detailed_celltype.sample_highlights.", make.names(x), ".nolabels.")

  png(filename = paste0(file_prefix, "png"), height = 20, width=40, units = 'in', res=300)
  print(spatial_panel & labs(title=NULL))
  dev.off()

  print(paste0(file_prefix, "png DONE!"))


  # Epithelial cell types
  epithelia.cell_types <- intersect(detailed.cell_types, as.character(metadf[metadf$cell_category=='Epithelia', 'cell_type.detailed']))
  epithelia.panel <- cowplot::plot_grid(plotlist=spatial.list[epithelia.cell_types], greedy = TRUE, nrow=2)
  file_prefix <- paste0(folder,"/detailed_celltype.epithelia.sample_highlights.", make.names(x), ".")
  png(filename = paste0(file_prefix, "png"), height = 8, width=24, units = 'in', res=300)
  print(epithelia.panel)
  dev.off()

  svglite::svglite(paste0(file_prefix, "svg"), height = 8, width=24)
  print(epithelia.panel)
  dev.off()


  stroma.cell_types <- intersect(detailed.cell_types, as.character(metadf[metadf$cell_category=='Stroma', 'cell_type.detailed']))
  stroma.panel <- cowplot::plot_grid(plotlist=spatial.list[stroma.cell_types], greedy = TRUE, nrow=2)
  file_prefix <- paste0(folder,"/detailed_celltype.stroma.sample_highlights.", make.names(x), ".")
  png(filename = paste0(file_prefix, "png"), height = 8, width=24, units = 'in', res=300)
  print(stroma.panel)
  dev.off()


  svglite::svglite(paste0(file_prefix, "svg"), height = 8, width=24)
  print(stroma.panel)
  dev.off()


  immune.cell_types <- intersect(detailed.cell_types, as.character(metadf[metadf$cell_category=='Immune', 'cell_type.detailed']))
  immune.panel <- cowplot::plot_grid(plotlist=spatial.list[immune.cell_types], greedy = TRUE, nrow=2)
  file_prefix <- paste0(folder,"/detailed_celltype.immune.sample_highlights.", make.names(x), ".")
  png(filename = paste0(file_prefix, "png"), height = 8, width=24, units = 'in', res=300)
  print(immune.panel)
  dev.off()


  svglite::svglite(paste0(file_prefix, "svg"), height = 8, width=24)
  print(immune.panel)
  dev.off()

}


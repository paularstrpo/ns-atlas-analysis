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
library(cowplot)
library(grid)
library(gridExtra)
library(circlize)
library(ComplexHeatmap)

set.seed(123)

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


stroma.colors <- c("#1C86EE", "#E31A1C", "#008B00", "#6A3D9A", "#FF7F00", "#8B4500",
                   "#FFD700", "#7EC0EE", "#FB9A99", "#90EE90", "#CAB2D6", "#FDBF6F",
                   "#EEE685", "#B03060", "#FF83FA")

names(stroma.colors) <- c("HEC","Retic Fib I","Perivasc Fib II","VEC","Papil Fib","Schwann",
                          "Peri","Perivasc Fib I","LEC","Retic Fib III","SM","Adipo",
                          "DP","Retic Fib II","DS")

immune.colors <-  c("#1C86EE", "#E31A1C", "#008B00", "#6A3D9A", "#FF7F00", "#8B4500",
                    "#96268D", "#7EC0EE", "#FB9A99", "#90EE90", "#CAB2D6", "#FDBF6F", "#EEE685")
names(immune.colors) <- c("Naïve T","CD4+ Th","CD4+ Treg","CD8+ Tc","Plasma","Cyc Imm",
                          "LC","CD1C+ DC","CLEC9A+ DC","CCR7+ DC","Mac","Mono",
                          "Mast")

## plot colors
global.colors <- c("lightcyan", "royalblue2")
detailed.cell_type.colors <- readRDS("data/reference/cell_type.detailed.color_palette.varibow.non_shuffled.RDS")
detailed.cell_types <- names(detailed.cell_type.colors)

## use this tutorial for being able to cluster the dotplots
## https://divingintogeneticsandgenomics.com/post/clustered-dotplot-for-single-cell-rnaseq/
lgd_list <- list(
  Legend( labels = c(0,0.25,0.5,0.75,1), title = "pt",
          graphics = list(
            function(x, y, w, h) grid.circle(x = x, y = y, r = 0  * unit(2, "mm"),
                                             gp = gpar(fill = "black")),
            function(x, y, w, h) grid.circle(x = x, y = y, r = sqrt(0.25) * unit(2, "mm"),
                                             gp = gpar(fill = "black")),
            function(x, y, w, h) grid.circle(x = x, y = y, r = sqrt(0.5) * unit(2, "mm"),
                                             gp = gpar(fill = "black")),
            function(x, y, w, h) grid.circle(x = x, y = y, r = sqrt(0.75) * unit(2, "mm"),
                                             gp = gpar(fill = "black")),
            function(x, y, w, h) grid.circle(x = x, y = y, r = 1 * unit(2, "mm"),
                                             gp = gpar(fill = "black")))
  ))



make.heatmap <- function(obj, marker.genes, assay, celltype.colors){

  dotplot <- DotPlot(obj, group.by = "cell_type.detailed",
                     features = marker.genes, dot.min = 0.05)

  df <- dotplot$data

  ### the matrix for the scaled expression
  expr.matrix <- df %>%
    select(-pct.exp, -avg.exp) %>%
    pivot_wider(names_from = id, values_from = avg.exp.scaled) %>%
    as.data.frame()

  row.names(expr.matrix) <- expr.matrix$features.plot
  expr.matrix <- expr.matrix[,-1] %>% as.matrix()

  ## reorder columns
  celltypes <- intersect(names(celltype.colors), colnames(expr.matrix))
  expr.matrix <- expr.matrix[, celltypes]


  ### the percent matrix
  percent.matrix <- df %>%
    select(-avg.exp, -avg.exp.scaled) %>%
    pivot_wider(names_from = id, values_from = pct.exp) %>%
    as.data.frame()

  row.names(percent.matrix) <- percent.matrix$features.plot
  percent.matrix <- percent.matrix[,-1] %>% as.matrix()

  percent.matrix <- percent.matrix[, celltypes]


  ## function for mapping the circle size
  layer_fun = function(j, i, x, y, w, h, fill){
    grid.rect(x = x, y = y, width = w, height = h,
              gp = gpar(col = NA, fill = NA))
    grid.circle(x=x,y=y,r= sqrt(pindex(percent.matrix, i, j)/100)  * unit(2, "mm"),
                gp = gpar(fill = col_fun(pindex(expr.matrix, i, j)), col = NA))}

  # quantiles <- quantile(expr.matrix, c(0.05, 0.5, 0.9, 0.95))

  ## any value that is greater than 2 will be mapped to the highest
  col_fun <- circlize::colorRamp2(c(-1, 0, 1, 1.5, 2), c("lightcyan", "#94efff" ,'dodgerblue', "blue1", "blue4"))

  ## make the top annotation
  colors <- celltype.colors[colnames(expr.matrix)]

  column_ha <- HeatmapAnnotation(
    " " = celltypes,
    col = list(" " = colors),
    na_col = "grey"
  )

  hm <- Heatmap(expr.matrix,
                column_title = assay,
                column_title_side = 'top',
                heatmap_legend_param=list(title="Expr"),
                col = col_fun,
                rect_gp = gpar(type = "none"),
                row_names_gp = gpar(fontsize = 10),
                layer_fun = layer_fun,
                row_names_side = "left",
                row_names_centered = FALSE,
                cluster_rows = FALSE,
                cluster_columns = FALSE,
                border = "black",
                bottom_annotation = column_ha
  )
  return(hm)
}



str.merfish <- qs::qread("data/merfish/BAYSOR/seurat_objects/Stroma.reclustered.seurat_object.annotated.filtered.QS")
str.scrna <- qs::qread("data/scrna/reclustering/Stroma.reclustered.seurat_object.annotated.QS")
str.scrna <- subset(str.scrna, cell_type.detailed != "Doublet")
gc()
str.merfish.barplot <- ggplot(str.merfish@meta.data) +
  aes(x = donor_id, fill = cell_type.detailed) +
  geom_bar(position='fill', color='black', show.legend=TRUE) +
  theme_classic() + theme +
  scale_fill_manual(values = stroma.colors) +
  scale_y_continuous(labels = scales::percent) +
   coord_fixed() +coord_flip()+
  labs(x = "MERFISH Donor ID", y = NULL, title = "Stroma")

str.scrna.barplot <- ggplot(str.scrna@meta.data) +
  aes(x = reorder(paste0(study_id, ".", donor_id), study_id), fill = cell_type.detailed) +
  geom_bar(position='fill', color='black', show.legend=TRUE) +
  theme_classic() + theme +
  scale_fill_manual(values = stroma.colors) +
  scale_y_continuous(labels = scales::percent) +
  coord_fixed() +coord_flip()+
  labs(x = "scRNA Donor ID", y = NULL, title = "Stroma")


str.panel <- (str.merfish.barplot / str.scrna.barplot) & NoLegend()

str.markers <- c(
  "ADIPOQ", "LPL", "G0S2", # ADIPO
  "ACTA2", "JAG1", "NOTCH3", "RGS5", #PERI
  "DES", "ITGA8", "MYH11",  #SM
  "VWF", "ACKR1", 'CD34', #EC
  "KDR", "CCL21", "TFPI", #LEC
  "MPZ", "SOX2", "SOX10", #Schwann
  "COL1A1", "PDGFRA",
  'COL11A1',"LRRC15","SFRP1", "POSTN", #DS
  "COCH","CRABP1",# DP
  "COMP", "COL18A1", "APCDD1", #PAPIL
  "SFRP2", "MMP2", "PRG4", "ANGPTL7", # RETIC
  "PTGDS",  "KLF5", "ITGA6", #RETIC CONT
  "APOE","CCL19", "CXCL12", "C3" #PERIVASC
  )

merfish.stroma_dotplot <- make.heatmap(str.merfish, assay = "MERFISH", marker.genes = str.markers, celltype.colors = stroma.colors)
scrna.stroma_dotplot <- make.heatmap(str.scrna, assay = 'scRNA', marker.genes = str.markers,  celltype.colors = stroma.colors)


svglite::svglite("figures/annotated/Stroma.MERFISH+scRNA.cell_type.detailed.reclust_dotplot.with_clustering.svg", height = 10, width = 8)
draw(merfish.stroma_dotplot + scrna.stroma_dotplot, annotation_legend_list = lgd_list, merge_legends = TRUE)
dev.off()

str.scrna.umap <- DimPlot(str.scrna, group.by = "cell_type.detailed", reduction = "reclust.umap", label=TRUE, alpha = 0.25, raster = TRUE) +
  NoAxes() +  scale_color_manual(values = stroma.colors) + labs(title=NULL) + NoLegend()


svglite::svglite("figures/annotated/scrna.stroma.reclustered.umap.annotated.svg", height = 6, width = 8)
str.scrna.umap
dev.off()


png("figures/annotated/scrna.stroma.reclustered.umap.annotated.png", height = 6, width = 8, units = 'in', res = 600)
DimPlot(str.scrna, group.by = "cell_type.detailed", reduction = "reclust.umap", label=FALSE, alpha = 0.25, raster = FALSE) +
  NoAxes() +  scale_color_manual(values = stroma.colors) + labs(title=NULL) + NoLegend()
dev.off()

## scRNA qc umaps
png("figures/annotated/scrna.stroma.reclustered.umap.sample_barcode.png", height = 6, width = 8, units = 'in', res = 600)
DimPlot(str.scrna, group.by = "sample_barcode", reduction = "reclust.umap", label=FALSE, alpha = 0.25, raster = FALSE) +
  NoAxes() + labs(title=NULL) + NoLegend()
dev.off()

png("figures/annotated/scrna.stroma.reclustered.umap.donor_id.png", height = 6, width = 8, units = 'in', res = 600)
DimPlot(str.scrna, group.by = "donor_id", reduction = "reclust.umap", label=FALSE, alpha = 0.25, raster = FALSE) +
  NoAxes() + labs(title=NULL) + NoLegend()
dev.off()


png("figures/annotated/scrna.stroma.reclustered.umap.study_id.png", height = 6, width = 8, units = 'in', res = 600)
DimPlot(str.scrna, group.by = "study_id", reduction = "reclust.umap", label=FALSE, alpha = 0.25, raster = FALSE) +
  NoAxes() + labs(title=NULL) + NoLegend()
dev.off()



## MERFISH qc umaps

str.merfish.umap <- DimPlot(str.merfish, group.by = "cell_type.detailed", reduction = "reclust_umap", label=TRUE, alpha = 0.25, raster = TRUE) +
  NoAxes() +  scale_color_manual(values = stroma.colors) + labs(title=NULL) + NoLegend()


svglite::svglite("figures/annotated/merfish.stroma.reclustered.umap.annotated.svg", height = 6, width = 8)
str.merfish.umap
dev.off()


png("figures/annotated/merfish.stroma.reclustered.umap.annotated.png", height = 6, width = 8, units = 'in', res = 600)
DimPlot(str.merfish, group.by = "cell_type.detailed", reduction = "reclust_umap", label=FALSE, alpha = 0.25, raster = FALSE) +
  NoAxes() +  scale_color_manual(values = stroma.colors) + labs(title=NULL) + NoLegend()
dev.off()

png("figures/annotated/merfish.stroma.reclustered.umap.sample_barcode.png", height = 6, width = 8, units = 'in', res = 600)
DimPlot(str.merfish, group.by = "sample_barcode", reduction = "reclust_umap", label=FALSE, alpha = 0.25, raster = FALSE) +
  NoAxes() + labs(title=NULL) + NoLegend()
dev.off()

png("figures/annotated/merfish.stroma.reclustered.umap.donor_id.png", height = 6, width = 8, units = 'in', res = 600)
DimPlot(str.merfish, group.by = "donor_id", reduction = "reclust_umap", label=FALSE, alpha = 0.25, raster = FALSE) +
  NoAxes() + labs(title=NULL) + NoLegend()
dev.off()


png("figures/annotated/merfish.stroma.reclustered.umap.batch.png", height = 6, width = 8, units = 'in', res = 600)
DimPlot(str.merfish, group.by = "run_name.short", reduction = "reclust_umap", label=FALSE, alpha = 0.25, raster = FALSE) +
  NoAxes() + labs(title=NULL) + NoLegend()
dev.off()


png("figures/annotated/merfish.stroma.reclustered.umap.anatomic_site.png", height = 6, width = 8, units = 'in', res = 600)
DimPlot(str.merfish, group.by = "anatomic_site", reduction = "reclust_umap", label=FALSE, alpha = 0.25, raster = FALSE) +
  NoAxes() + labs(title=NULL) + NoLegend()
dev.off()


epithelia.markers <-  c("TOP2A", "MKI67",
                        "COL17A1", "DST",
                        "GATA3", "KLF5",
                        "NOTCH3", "DSC1", "S100A8",
                        "FLG","PRDM1",'KLK7',
                        "MLANA","MITF",
                        "PTN", "IGFL2",'DIO2',
                        "FGFR1","COL18A1",
                        "THBS1","SOX9",
                        "TKT", "KRT79",
                        "KRT28", "KRT74",
                        "SEMA3C", "INHBB",'KRT77',
                        "SFRP1", "DNER", 'CLDN10',
                        'ACTA2','MYH11')


epi.merfish <- qs::qread("data/merfish/BAYSOR/seurat_objects/Epithelia.reclustered.seurat_object.annotated.filtered.QS")
epi.cell_types <- c("Cyc KC", "Bas KC", "Spn KC I", "Spn KC II", "Grn KC", "Melano", "Bas Inf", "Diff Inf", "Bulge",
                    "ORS Basal", "ORS Suprabasal", "Bas Seb", "Diff Seb","IRS/HS", "Ecc Duct", "Ecc Gland", "Ecc Myoepi" )
epi.merfish$cell_type.detailed <- factor(epi.merfish$cell_type.detailed, levels = epi.cell_types)

epi.scrna <- qs::qread("data/scrna/reclustering/Epithelia.reclustered.seurat_object.annotated.QS")
epi.scrna <- subset(epi.scrna, cell_type.detailed != "Doublet")
epithelia.colors <- detailed.cell_type.colors[levels(epi.merfish$cell_type.detailed)]

epi.merfish.barplot <- ggplot(epi.merfish@meta.data) +
  aes(x = donor_id, fill = cell_type.detailed) +
  geom_bar(position='fill', color='black', show.legend=TRUE) +
  theme_classic() + theme +
  scale_fill_manual(values = epithelia.colors) +
  scale_y_continuous(labels = scales::percent) +
   coord_fixed() +coord_flip()+
  labs(x = "MERFISH Donor ID", y = NULL, title = "Epithelia")

epi.scrna.barplot <- ggplot(epi.scrna@meta.data) +
  aes(x = reorder(paste0(study_id, ".", donor_id), study_id), fill = cell_type.detailed) +
  geom_bar(position='fill', color='black', show.legend=TRUE) +
  theme_classic() + theme +
  scale_fill_manual(values = epithelia.colors) +
  scale_y_continuous(labels = scales::percent) +
  coord_fixed() +coord_flip()+
  labs(x = "scRNA Donor ID", y = NULL, title = "Epithelia")


epi.panel <- (epi.merfish.barplot / epi.scrna.barplot) & NoLegend()


merfish.epi_dotplot <- make.heatmap(epi.merfish, assay = "MERFISH", marker.genes = epithelia.markers, celltype.colors = epithelia.colors)
scrna.epi_dotplot <- make.heatmap(epi.scrna, assay = 'scRNA', marker.genes = epithelia.markers,  celltype.colors = epithelia.colors)


svglite::svglite("figures/annotated/Epithelia.MERFISH+scRNA.cell_type.detailed.reclust_dotplot.with_clustering.svg", height = 10, width = 8)
draw(merfish.epi_dotplot + scrna.epi_dotplot, annotation_legend_list = lgd_list, merge_legends = TRUE)
dev.off()


png("figures/annotated/Epithelia.MERFISH+scRNA.cell_type.detailed.reclust_dotplot.with_clustering.png", height = 10, width = 8, res=300, units='in')
draw(merfish.epi_dotplot + scrna.epi_dotplot, annotation_legend_list = lgd_list, merge_legends = TRUE)
dev.off()


epi.scrna.umap <- DimPlot(epi.scrna, group.by = "cell_type.detailed",  reduction='reclust.umap', label = TRUE, alpha = 0.25, raster=TRUE) +
  scale_color_manual(values = epithelia.colors) +
  NoAxes() + NoLegend() + labs(title=NULL)


svglite::svglite("figures/annotated/scrna.epithelia.reclustered.umap.annotated.svg", height = 6, width = 8)
epi.scrna.umap
dev.off()

epi.scrna.umap.nolab <- DimPlot(epi.scrna, group.by = "cell_type.detailed",  reduction='reclust.umap', label = FALSE, alpha = 0.25, raster=FALSE) +
  scale_color_manual(values = epithelia.colors) +
  NoAxes() + NoLegend() + labs(title=NULL)

png("figures/annotated/scrna.epithelia.reclustered.umap.annotated.png", height = 6, width = 8, units='in', res = 600)
epi.scrna.umap.nolab
dev.off()



## scRNA qc umaps
png("figures/annotated/scrna.epithelia.reclustered.umap.sample_barcode.png", height = 6, width = 8, units = 'in', res = 600)
DimPlot(epi.scrna, group.by = "sample_barcode", reduction = "reclust.umap", label=FALSE, alpha = 0.25, raster = FALSE) +
  NoAxes() + labs(title=NULL) + NoLegend()
dev.off()

png("figures/annotated/scrna.epithelia.reclustered.umap.donor_id.png", height = 6, width = 8, units = 'in', res = 600)
DimPlot(epi.scrna, group.by = "donor_id", reduction = "reclust.umap", label=FALSE, alpha = 0.25, raster = FALSE) +
  NoAxes() + labs(title=NULL) + NoLegend()
dev.off()


png("figures/annotated/scrna.epithelia.reclustered.umap.study_id.png", height = 6, width = 8, units = 'in', res = 600)
DimPlot(epi.scrna, group.by = "study_id", reduction = "reclust.umap", label=FALSE, alpha = 0.25, raster = FALSE) +
  NoAxes() + labs(title=NULL) + NoLegend()
dev.off()



## MERFISH qc umaps


epi.merfish.umap <- DimPlot(epi.merfish, group.by = "cell_type.detailed",  reduction='reclust_umap', label = TRUE, alpha = 0.25, raster=TRUE) +
  scale_color_manual(values = epithelia.colors) +
  NoAxes() + NoLegend() + labs(title=NULL)


svglite::svglite("figures/annotated/merfish.epithelia.reclustered.umap.annotated.svg", height = 6, width = 8)
epi.merfish.umap
dev.off()

epi.merfish.umap.nolab <- DimPlot(epi.merfish, group.by = "cell_type.detailed",  reduction='reclust_umap', label = FALSE, alpha = 0.25, raster=FALSE) +
  scale_color_manual(values = epithelia.colors) +
  NoAxes() + NoLegend() + labs(title=NULL)

png("figures/annotated/merfish.epithelia.reclustered.umap.annotated.png", height = 6, width = 8, units='in', res = 600)
epi.merfish.umap.nolab
dev.off()


png("figures/annotated/merfish.epithelia.reclustered.umap.sample_barcode.png", height = 6, width = 8, units = 'in', res = 600)
DimPlot(epi.merfish, group.by = "sample_barcode", reduction = "reclust_umap", label=FALSE, alpha = 0.25, raster = FALSE) +
  NoAxes() + labs(title=NULL) + NoLegend()
dev.off()

png("figures/annotated/merfish.epithelia.reclustered.umap.donor_id.png", height = 6, width = 8, units = 'in', res = 600)
DimPlot(epi.merfish, group.by = "donor_id", reduction = "reclust_umap", label=FALSE, alpha = 0.25, raster = FALSE) +
  NoAxes() + labs(title=NULL) + NoLegend()
dev.off()


png("figures/annotated/merfish.epithelia.reclustered.umap.batch.png", height = 6, width = 8, units = 'in', res = 600)
DimPlot(epi.merfish, group.by = "run_name.short", reduction = "reclust_umap", label=FALSE, alpha = 0.25, raster = FALSE) +
  NoAxes() + labs(title=NULL) + NoLegend()
dev.off()


png("figures/annotated/merfish.epithelia.reclustered.umap.anatomic_site.png", height = 6, width = 8, units = 'in', res = 600)
DimPlot(epi.merfish, group.by = "anatomic_site", reduction = "reclust_umap", label=FALSE, alpha = 0.25, raster = FALSE) +
  NoAxes() + labs(title=NULL) + NoLegend()
dev.off()

rm(epi.merfish, epi.scrna, epi.scrna.umap, epi.scrna.umap.nolab); gc()


imm.merfish <- qs::qread("data/merfish/BAYSOR/seurat_objects/Immune.reclustered.seurat_object.annotated.filtered.QS")
imm.scrna <- qs::qread("data/scrna/reclustering/Immune.reclustered.seurat_object.annotated.QS")
imm.scrna <- subset(imm.scrna, cell_type.detailed != "Doublet")
immune.cell_types <- c( "Naïve T", "CD4+ Th", "CD4+ Treg", "CD8+ Tc", 'Plasma',
                        'Cyc Imm', "LC", "CD1C+ DC","CLEC9A+ DC","CCR7+ DC",
                        "Mac", 'Mono',"Mast")

imm.merfish.barplot <- ggplot(imm.merfish@meta.data) +
  aes(x = donor_id, fill = cell_type.detailed) +
  geom_bar(position='fill', color='black', show.legend=TRUE) +
  theme_classic() + theme +
  scale_fill_manual(values = immune.colors) +
  scale_y_continuous(labels = scales::percent) +
  coord_fixed() +coord_flip()+
  labs(x = "MERFISH Donor ID", y = NULL, title = "Immune")

imm.scrna.barplot <- ggplot(imm.scrna@meta.data) +
  aes(x = reorder(paste0(study_id, ".", donor_id), study_id), fill = cell_type.detailed) +
  geom_bar(position='fill', color='black', show.legend=TRUE) +
  theme_classic() + theme +
  scale_fill_manual(values = immune.colors) +
  scale_y_continuous(labels = scales::percent) +
  coord_fixed() +coord_flip()+
  labs(x = "scRNA Donor ID", y = NULL, title = "Immune")


imm.panel <- (imm.merfish.barplot / imm.scrna.barplot) & NoLegend()


imm.panel

imm.markers <- c(
  "CD3E", 'IL2RG',
  "CCR4","IL7R", "CD69",
  "CD40LG", "CD4","RORC",
  "CTLA4", "FOXP3","SELL", 'CD27',
  "CD8A", "CCL5",
  "GZMA", "NKG7",
  "JCHAIN",
  "MKI67", "TOP2A",
  "CD1A", "CD207",
  "CD1C", "CLEC10A","CLEC9A",
  'IRF8', 'CD40', "CCR7","CD83",
  "C1QA", "MRC1", "CD163",
  "CD14",  'C5AR1', "ITGAX",
  "IL1RL1", "KIT", "HDC")


merfish.imm_dotplot <-  make.heatmap(imm.merfish, assay = "MERFISH", marker.genes = imm.markers, celltype.colors = immune.colors)
scrna.imm_dotplot <-  make.heatmap(imm.scrna, assay = "scRNA", marker.genes = imm.markers, celltype.colors =  immune.colors)


svglite::svglite("figures/annotated/Immune.MERFISH+scRNA.cell_type.detailed.reclust_dotplot.with_clustering.svg", height = 10, width = 8)
draw(merfish.imm_dotplot + scrna.imm_dotplot, annotation_legend_list = lgd_list, merge_legends = TRUE)
dev.off()


png("figures/annotated/Immune.MERFISH+scRNA.cell_type.detailed.reclust_dotplot.with_clustering.png", height = 10, width = 8, units='in', res=300)
draw(merfish.imm_dotplot + scrna.imm_dotplot, annotation_legend_list = lgd_list, merge_legends = TRUE)
dev.off()


imm.scrna.umap <- DimPlot(imm.scrna, group.by = "cell_type.detailed",  reduction='reclust.umap', label = TRUE, alpha = 0.25, raster=TRUE) +
  scale_color_manual(values = immune.colors) +
  NoAxes() + NoLegend() + labs(title=NULL)


svglite::svglite("figures/annotated/scrna.immune.reclustered.umap.annotated.svg", height = 6, width = 8)
imm.scrna.umap
dev.off()

imm.scrna.umap.nolab <- DimPlot(imm.scrna, group.by = "cell_type.detailed",  reduction='reclust.umap', label = FALSE, alpha = 0.25, raster=FALSE) +
  scale_color_manual(values = immune.colors) +
  NoAxes() + NoLegend() + labs(title=NULL)

png("figures/annotated/scrna.immune.reclustered.umap.annotated.png", height = 6, width = 8, units='in', res = 600)
imm.scrna.umap.nolab
dev.off()


## scRNA qc umaps
png("figures/annotated/scrna.immune.reclustered.umap.sample_barcode.png", height = 6, width = 8, units = 'in', res = 600)
DimPlot(imm.scrna, group.by = "sample_barcode", reduction = "reclust.umap", label=FALSE, alpha = 0.25, raster = FALSE) +
  NoAxes() + labs(title=NULL) + NoLegend()
dev.off()

png("figures/annotated/scrna.immune.reclustered.umap.donor_id.png", height = 6, width = 8, units = 'in', res = 600)
DimPlot(imm.scrna, group.by = "donor_id", reduction = "reclust.umap", label=FALSE, alpha = 0.25, raster = FALSE) +
  NoAxes() + labs(title=NULL) + NoLegend()
dev.off()


png("figures/annotated/scrna.immune.reclustered.umap.study_id.png", height = 6, width = 8, units = 'in', res = 600)
DimPlot(imm.scrna, group.by = "study_id", reduction = "reclust.umap", label=FALSE, alpha = 0.25, raster = FALSE) +
  NoAxes() + labs(title=NULL) + NoLegend()
dev.off()



## MERFISH qc umaps



imm.merfish.umap <- DimPlot(imm.merfish, group.by = "cell_type.detailed",  reduction='reclust_umap', label = TRUE, alpha = 0.25, raster=TRUE) +
  scale_color_manual(values = immune.colors) +
  NoAxes() + NoLegend() + labs(title=NULL)


svglite::svglite("figures/annotated/merfish.immune.reclustered.umap.annotated.svg", height = 6, width = 8)
imm.merfish.umap
dev.off()

imm.merfish.umap.nolab <- DimPlot(imm.merfish, group.by = "cell_type.detailed",  reduction='reclust_umap', label = FALSE, alpha = 0.25, raster=FALSE) +
  scale_color_manual(values = immune.colors) +
  NoAxes() + NoLegend() + labs(title=NULL)

png("figures/annotated/merfish.immune.reclustered.umap.annotated.png", height = 6, width = 8, units='in', res = 600)
imm.merfish.umap.nolab
dev.off()



png("figures/annotated/merfish.immune.reclustered.umap.sample_barcode.png", height = 6, width = 8, units = 'in', res = 600)
DimPlot(imm.merfish, group.by = "sample_barcode", reduction = "reclust_umap", label=FALSE, alpha = 0.25, raster = FALSE) +
  NoAxes() + labs(title=NULL) + NoLegend()
dev.off()

png("figures/annotated/merfish.immune.reclustered.umap.donor_id.png", height = 6, width = 8, units = 'in', res = 600)
DimPlot(imm.merfish, group.by = "donor_id", reduction = "reclust_umap", label=FALSE, alpha = 0.25, raster = FALSE) +
  NoAxes() + labs(title=NULL) + NoLegend()
dev.off()


png("figures/annotated/merfish.immune.reclustered.umap.batch.png", height = 6, width = 8, units = 'in', res = 600)
DimPlot(imm.merfish, group.by = "run_name.short", reduction = "reclust_umap", label=FALSE, alpha = 0.25, raster = FALSE) +
  NoAxes() + labs(title=NULL) + NoLegend()
dev.off()


png("figures/annotated/merfish.immune.reclustered.umap.anatomic_site.png", height = 6, width = 8, units = 'in', res = 600)
DimPlot(imm.merfish, group.by = "anatomic_site", reduction = "reclust_umap", label=FALSE, alpha = 0.25, raster = FALSE) +
  NoAxes() + labs(title=NULL) + NoLegend()
dev.off()




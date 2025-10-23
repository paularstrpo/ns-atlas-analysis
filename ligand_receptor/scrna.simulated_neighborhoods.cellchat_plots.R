## scRNA CellChat Plots

library(Seurat)
library(Matrix)
library(tidyverse)
library(CellChat)

setwd("D:/Dropbox/_projects/NormalSkinAtlas")

# interaction_list <- qs::qread("data/scrna/ligand_receptor/scrna.neighborhood.cellchat.result_tables.list.QS")
# pathway_list <- qs::qread("data/scrna/ligand_receptor/scrna.neighborhood.cellchat.result_tables.pathways.list.QS")
cellchat_list <- qs::qread(file = "data/scrna/ligand_receptor/scrna.neighborhood.cellchat_object.list.QS")
neighborhoods <- names(cellchat_list)
neighborhood.colors <- readRDS("data/reference/neighborhood.color_palette.final.RDS")
names(neighborhoods) <- names(cellchat_list)
cell_type.detailed.colors <- readRDS("data/reference/cell_type.detailed.color_palette.final.RDS")
gc()

for (nhood in neighborhoods){

  file_prefix <- paste0("figures/scrna/ligand_receptor/scrna.neighborhood_cellchat.", make.names(nhood), ".pathway_circos")
  svglite::svglite(paste0(file_prefix, ".svg"), width = 16, height=14)
  netVisual_chord_gene(cellchat_list[[nhood]], slot.name ="netP", small.gap = 0.2, lab.cex = 0.6,
                       title.name = nhood, color.use = cell_type.detailed.colors, show.legend = FALSE)
  dev.off()

  png(paste0(file_prefix, ".png"), width = 16, height=14, units='in', res=300)
  netVisual_chord_gene(cellchat_list[[nhood]], slot.name ="netP", small.gap = 0.2, lab.cex = 0.6,
                       title.name = nhood)
  dev.off()
}


for (nhood in neighborhoods){

  ## plot top signaling pathways per neighborhood
  file_prefix <- paste0("figures/scrna/ligand_receptor/scrna.neighborhood_cellchat.", make.names(nhood), ".pathway_circos")
  svglite::svglite(paste0(file_prefix, ".svg"), width = 16, height=14)
  netVisual_chord_gene(cellchat_list[[nhood]], slot.name ="netP", small.gap = 0.2, lab.cex = 0.6,
                       title.name = nhood, color.use = cell_type.detailed.colors, show.legend = FALSE)
  dev.off()

  png(paste0(file_prefix, ".png"), width = 16, height=14, units='in', res=300)
  netVisual_chord_gene(cellchat_list[[nhood]], slot.name ="netP", small.gap = 0.2, lab.cex = 0.6, title.name = nhood)
  dev.off()


  ## plot number of interactions per neighborhood by cell types
  groupSize <- as.numeric(table(cellchat_list[[nhood]]@idents))

  file_prefix <- paste0("figures/scrna/ligand_receptor/scrna.neighborhood_cellchat.", make.names(nhood), ".number_interactions")
  svglite::svglite(paste0(file_prefix, ".svg"), width = 16, height=14)
  netVisual_circle(cellchat_list[[nhood]]@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = nhood,
                   color.use = cell_type.detailed.colors)
  dev.off()


  png(paste0(file_prefix, ".png"), width = 16, height=14, units='in', res=300)
  netVisual_circle(cellchat_list[[nhood]]@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = nhood,
                   color.use = cell_type.detailed.colors)
  dev.off()
}



for (nhood in neighborhoods){
  ## plot number of interactions per neighborhood by cell types

  file_prefix <- paste0("figures/scrna/ligand_receptor/scrna.neighborhood_cellchat.", make.names(nhood), ".number_interactions")
  svglite::svglite(paste0(file_prefix, ".svg"), width = 16, height=14)
  netVisual_circle(cellchat_list[[nhood]]@net$count, weight.scale = TRUE, label.edge= FALSE, title.name = nhood,
                   color.use = cell_type.detailed.colors)
  dev.off()


  png(paste0(file_prefix, ".png"), width = 16, height=14, units='in', res=300)
  netVisual_circle(cellchat_list[[nhood]]@net$count, weight.scale = TRUE, label.edge= FALSE, title.name = nhood,
                   color.use = cell_type.detailed.colors)
  dev.off()
}

for (nhood in neighborhoods){

  cell_props <- sort(table(cellchat_list[[nhood]]@idents), decreasing = TRUE)/length(cellchat_list[[nhood]]@idents)*100
  cell_types <-  names(cell_props[cell_props>5])

  file_prefix <- paste0("figures/scrna/ligand_receptor/scrna.neighborhood_cellchat.", make.names(nhood), ".pathway_circos.abundance_filt.")
  svglite::svglite(paste0(file_prefix, ".svg"), width = 16, height=14)
  netVisual_chord_gene(cellchat_list[[nhood]], slot.name ="netP",
                       small.gap = 0.2, lab.cex = 0.6, thresh = 0.01,
                       color.use = cell_type.detailed.colors[cell_types],
                       title.name = nhood, show.legend = TRUE)
  dev.off()

  png(paste0(file_prefix, ".png"), width = 16, height=14, units='in', res=300)
  netVisual_chord_gene(cellchat_list[[nhood]], slot.name ="netP",
                       small.gap = 0.2, lab.cex = 0.6, thresh = 0.01,
                       color.use = cell_type.detailed.colors[cell_types],
                       title.name = nhood, show.legend = TRUE)
  dev.off()
}


for (nhood in neighborhoods){

  cell_props <- sort(table(cellchat_list[[nhood]]@idents), decreasing = TRUE)/length(cellchat_list[[nhood]]@idents)*100
  cell_types <-  names(cell_props[cell_props>5])

  file_prefix <- paste0("figures/scrna/ligand_receptor/scrna.neighborhood_cellchat.", make.names(nhood), ".LR_circos.abundance_filt.")
  svglite::svglite(paste0(file_prefix, ".svg"), width = 16, height=14)
  netVisual_chord_gene(cellchat_list[[nhood]],
                       small.gap = 0.2, lab.cex = 0.6, thresh = 0.01,
                       color.use = cell_type.detailed.colors[cell_types],
                       title.name = nhood, show.legend = TRUE)
  dev.off()

  png(paste0(file_prefix, ".png"), width = 16, height=14, units='in', res=300)
  netVisual_chord_gene(cellchat_list[[nhood]],
                       small.gap = 0.2, lab.cex = 0.6, thresh = 0.01,
                       color.use = cell_type.detailed.colors[cell_types],
                       title.name = nhood, show.legend = TRUE)
  dev.off()
}


library(Matrix)
library(Seurat)
setwd("D:/Dropbox/_projects/NormalSkinAtlas")
obj <- qs::qread("data/merfish/BAYSOR/seurat_objects/ns-atlas.merfish_baysor.scanvi_integrated.cellcharter.annotated.seurat_object.QS")
obj <- subset(obj,
              subset = cell_type.detailed != "Doublet" &
                cell_type.detailed != "Unknown" &
                tissue_compartment != "outside tissue")
obj <- UpdateSeuratObject(obj)
Idents(obj) <- "neighborhood.detailed"
meta <- obj@meta.data

# restrict our cellchat analysis to the cell types that were analyzed for the MERFISH neighborhood cellchat
# this will effectively simulate the neighborhoods in the scRNA-seq for LR purposes.
# use the graph layouts of the merfish cellchat to get the list of cell types
library(tidyverse)
library(DirichletReg)

cell_counts <- data.frame(table(meta$cell_type.detailed, meta$neighborhood))
colnames(cell_counts) <- c("cell_type.detailed", "neighborhood", "cell_count")

x <- table(meta$neighborhood, as.character(meta$cell_type.detailed))
class(x) <- 'matrix'

# makes cell type proportions
drmatrix <- DR_data(x)
cell_prop <- t(as.matrix(drmatrix))
class(cell_prop) <- 'matrix'
cell_prop <- data.frame(cell_prop)
cell_prop <- na.omit(cell_prop)

## only do signaling on cell types making up at least 5% of total neighborhood proportion
neighborhoods <- levels(obj@meta.data$neighborhood)
neighborhood.list <- apply(cell_prop, MARGIN = 2, function(x){names(sort(x[x>=0.01], decreasing = TRUE))})
names(neighborhood.list) <- gsub("\\.", "-", names(neighborhood.list))

library(CellChat)

## https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/FAQ_on_applying_CellChat_to_spatial_transcriptomics_data.html#seqfishmerfishstarmap
conversion.factor = 1
contact.range = 50
max.distance = 200 # based on literature
spot.size = 15
tol <- spot.size/2
spatial.factors = data.frame(ratio = conversion.factor, tol = tol)

CellChatDB <- CellChatDB.human # use CellChatDB.human if running on human data
# subset the expression data of signaling genes for saving computation cost
CellChatDB.use <- subsetDB(CellChatDB, key = "annotation")

object.list <- list()
interaction_list <- list()
pathway_list <- list()
options(future.globals.maxSize=4e09)
neighborhoods <- levels(obj@meta.data$neighborhood)
for (nhood in neighborhoods){
  print(nhood)

  Idents(obj) <- "neighborhood"
  obj_filt <- subset(obj, idents = nhood)

  Idents(obj_filt) <- "cell_type.detailed"
  obj_filt <- subset(obj_filt, idents = neighborhood.list[[nhood]])

  ## include samples with at least 100 cells in the given neighborhood
  samples <- table(obj_filt$sample_barcode)
  print(sort(samples))
  samples <- names(samples[samples >= 100])
  Idents(obj_filt) <- "sample_barcode"
  obj_filt <- subset(obj_filt, idents = samples)

  data.input <- GetAssayData(obj_filt, layer = 'data')
  spatial.locs <- data.frame(Embeddings(obj_filt, reduction = 'spatial'))

  meta <- obj_filt@meta.data[, c("cell_barcode", "sample_barcode", "cell_type.detailed")]
  colnames(meta) <- c("cell_barcode", "samples", "labels")
  meta$samples <- factor(meta$samples)
  meta$labels <- factor(meta$labels)

  cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels",
                             datatype = "spatial", coordinates = spatial.locs,
                             spatial.factors = spatial.factors)
  rm(data.input, meta, spatial.locs, obj_filt); invisible(gc())

  # subset the expression data of signaling genes for saving computation cost
  cellchat@DB <- CellChatDB.use
  cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database

  # https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/FAQ_on_applying_CellChat_to_spatial_transcriptomics_data.html#application-to-a-dataset-with-a-small-panel-of-genes
  cellchat <- identifyOverExpressedGenes(cellchat, do.DE = FALSE, min.cells = 10)
  cellchat <- identifyOverExpressedInteractions(cellchat)

  cellchat <- computeCommunProb(cellchat,
                                type = 'truncatedMean',
                                trim=0.1,
                                distance.use = FALSE,
                                scale.distance = NULL,
                                interaction.range = max.distance,
                                contact.dependent = TRUE,
                                contact.range = contact.range,
                                seed.use = 123
                                )
  cellchat <- filterCommunication(cellchat)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)

  if(length(cellchat@netP$pathways) > 0){

    cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "net")
    cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

    df.net <- data.frame(subsetCommunication(cellchat))
    df.net$neighborhood <- nhood
    interaction_list[[nhood]] <- df.net
    write.csv(df.net, file = paste0('data/merfish/BAYSOR/ligand_receptor/CellChat/neighborhood/neighborhood_', nhood, '.cellchat.interaction_results.csv'))
    qs::qsave(cellchat, file = paste0('data/merfish/BAYSOR/ligand_receptor/CellChat/neighborhood/neighborhood_', nhood, '.cellchat_object.QS'))
    object.list[[nhood]] <- cellchat

    df.net <- data.frame(subsetCommunication(cellchat, slot.name = 'netP'))
    df.net$neighborhood <- nhood
    pathway_list[[nhood]] <- df.net
    write.csv(df.net, file =
                paste0('data/merfish/BAYSOR/ligand_receptor/CellChat/neighborhood/neighborhood_', nhood,
                       '.cellchat_results.pathways.csv'))


    rm(df.net); invisible(gc())
  }
  rm(cellchat); invisible(gc())
}
qs::qsave(interaction_list, file = "data/merfish/BAYSOR/ligand_receptor/CellChat/neighborhood/neighborhood.cellchat.result_tables.list.QS")
qs::qsave(pathway_list, file = "data/merfish/BAYSOR/ligand_receptor/CellChat/neighborhood/neighborhood.cellchat.result_tables.pathways.list.QS")
qs::qsave(object.list, file = "data/merfish/BAYSOR/ligand_receptor/CellChat/neighborhood/neighborhood.cellchat_object.list.QS")
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
qs::qsave(cellchat, "data/merfish/BAYSOR/ligand_receptor/CellChat/neighborhood/neighborhood.cellchat_object.merged.QS")


df <- netVisual_bubble(cellchat, comparison = 1:length(object.list), angle.x = 45, remove.isolate = TRUE, return.data = TRUE)[[1]]
write.csv(df, file = "data/merfish/BAYSOR/ligand_receptor/CellChat/neighborhood/neighborhood.net_ineractions.table.csv")


pathway_rank <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL,
                        targets.use = NULL, stacked = TRUE, do.stat = FALSE, comparison = 1:length(object.list), return.data = TRUE)[[1]][, 1:4]
write.csv(pathway_rank, file = "data/merfish/BAYSOR/ligand_receptor/CellChat/neighborhood/neighborhood.pathway_ranking.weight.table.csv")


pathway_rank <- rankNet(cellchat, mode = "comparison", measure = "count", sources.use = NULL,
                        targets.use = NULL, stacked = TRUE, do.stat = FALSE, comparison = 1:length(object.list), return.data = TRUE)[[1]][, 1:4]
write.csv(pathway_rank, file = "data/merfish/BAYSOR/ligand_receptor/CellChat/neighborhood/neighborhood.pathway_ranking.count.table.csv")

object.list <- qs::qread("data/merfish/BAYSOR/ligand_receptor/CellChat/neighborhood/neighborhood.cellchat_object.list.QS")


cellchat <- qs::qread("data/merfish/BAYSOR/ligand_receptor/CellChat/neighborhood/neighborhood.cellchat_object.merged.QS")

neighborhoods <- names(object.list)
neighborhood.colors <- pals::glasbey(10)


interactions.barplot <- compareInteractions(cellchat, group = names(object.list), color.use = neighborhood.colors) + coord_flip() + NoLegend() |
  compareInteractions(cellchat, group = names(object.list), measure = 'weight', color.use = neighborhood.colors) + coord_flip() + NoLegend()


svglite::svglite("figures/merfish/BAYSOR/ligand_receptor/CellChat/neighborhood/number_interactions.barplot.svg", height=5, width=10)
interactions.barplot
dev.off()

png("figures/merfish/BAYSOR/ligand_receptor/CellChat/neighborhood/number_interactions.barplot.png", height=5, width=10, units='in', res=300)
interactions.barplot
dev.off()


# interactions.barplot
netp <- rankNet(cellchat, mode = "comparison", measure = "weight",color.use = neighborhood.colors, sources.use = NULL, targets.use = NULL, stacked = TRUE, do.stat = FALSE, comparison = 1:length(object.list), return.data = TRUE)[[1]]
top_pathways_per_neighborhod <- rankNet(cellchat, mode = "comparison", measure = "weight",color.use = neighborhood.colors, sources.use = NULL, targets.use = NULL, stacked = TRUE, do.stat = FALSE, comparison = 1:length(object.list))
top_pathways_per_neighborhod


svglite::svglite("figures/merfish/BAYSOR/ligand_receptor/CellChat/neighborhood/top_pathways_per_neighborhood.ranknet_propbarplot.svg", height=7, width=10)
top_pathways_per_neighborhod
dev.off()

png("figures/merfish/BAYSOR/ligand_receptor/CellChat/neighborhood/top_pathways_per_neighborhood.ranknet_propbarplot.png", height=7, width=5, units='in', res=300)
top_pathways_per_neighborhod
dev.off()

cell_type.detailed.colors <- readRDS("data/reference/cell_type.detailed.color_palette.varibow.non_shuffled.RDS")
cell_type.detailed.colors <- cell_type.detailed.colors[2:length(cell_type.detailed.colors)]
detailed.cell_types <- names(cell_type.detailed.colors)

#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
num.link <- unlist(sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)}))
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax, color.use = cell_type.detailed.colors)
}

signaling_role.panel <- patchwork::wrap_plots(plots = gg, nrow = 2)

svglite::svglite('figures/merfish/BAYSOR/ligand_receptor/CellChat/neighborhood/signaling_role.per_celltype.scatterplot_panel.svg', height=8, width=24)
signaling_role.panel
dev.off()

png('figures/merfish/BAYSOR/ligand_receptor/CellChat/neighborhood/signaling_role.per_celltype.scatterplot_panel.png', height=8, width=24, res=300, units='in')
signaling_role.panel
dev.off()


library(ComplexHeatmap)

# combining all the identified signaling pathways from different datasets
pathway.union <- unique(unlist(lapply(object.list, function(x){x@netP$pathways})))

heatmap_list <- lapply(names(object.list), function(site){
  ht <- netAnalysis_signalingRole_heatmap(object.list[[site]], pattern = "all", signaling = pathway.union, title = site, width = 10, height=8)
  ht <- grid.grabExpr(draw(ht))
  return(ht)
})

png("figures/merfish/BAYSOR/ligand_receptor/CellChat/neighborhood/neighborhood.signaling_role.all_patterns.heatmap.panel.png", height=8, width = 24, res = 150, units = 'in')
cowplot::plot_grid(plotlist = heatmap_list, nrow=2)
dev.off()


svglite::svglite("figures/merfish/BAYSOR/ligand_receptor/CellChat/neighborhood/neighborhood.signaling_role.all_patterns.heatmap.panel.svg", height=12, width = 36)
cowplot::plot_grid(plotlist = heatmap_list, nrow=2)
dev.off()

heatmap_list <- lapply(names(object.list), function(site){
  ht <- netAnalysis_signalingRole_heatmap(object.list[[site]], pattern = "outgoing", signaling = pathway.union, title = site, width = 10, height=8)
  ht <- grid.grabExpr(draw(ht))
  return(ht)
})

png("figures/merfish/BAYSOR/ligand_receptor/CellChat/neighborhood/neighborhood.signaling_role.outgoing_patterns.heatmap.panel.png", height=8, width = 24, res = 150, units = 'in')
cowplot::plot_grid(plotlist = heatmap_list, nrow=2)
dev.off()


svglite::svglite("figures/merfish/BAYSOR/ligand_receptor/CellChat/neighborhood/neighborhood.outgoing_role.all_patterns.heatmap.panel.svg", height=8, width = 24)
cowplot::plot_grid(plotlist = heatmap_list, nrow=2)
dev.off()

heatmap_list <- lapply(names(object.list), function(site){
  ht <- netAnalysis_signalingRole_heatmap(object.list[[site]], pattern = "incoming", signaling = pathway.union, title = site, width = 10, height=8)
  ht <- grid.grabExpr(draw(ht))
  return(ht)
})

png("figures/merfish/BAYSOR/ligand_receptor/CellChat/neighborhood/neighborhood.signaling_role.incoming_patterns.heatmap.panel.png", height=8, width = 24, res = 300, units = 'in')
cowplot::plot_grid(plotlist = heatmap_list, nrow=2)
dev.off()


svglite::svglite("figures/merfish/BAYSOR/ligand_receptor/CellChat/neighborhood/neighborhood.incoming_role.all_patterns.heatmap.panel.svg", height=8, width = 24)
cowplot::plot_grid(plotlist = heatmap_list, nrow=2)
dev.off()


library(Seurat)
library(Matrix)
library(tidyverse)
library(CellChat)
library(igraph)
library(tidygraph)
library(ggraph)

setwd("D:/Dropbox/_projects/NormalSkinAtlas")
cellchat_list <- object.list
cell_type.detailed.colors <- readRDS("data/reference/cell_type.detailed.color_palette.varibow.non_shuffled.RDS")
cell_type.detailed.colors <- cell_type.detailed.colors[2:length(cell_type.detailed.colors)]
detailed.cell_types <- names(cell_type.detailed.colors)
meta <- obj@meta.data
meta$neighborhood.detailed <- factor(meta$neighborhood, levels = neighborhoods)

# get neighborhood proportion data to have info for the nodes
nhood.residents <- lapply(neighborhoods, function(nhood){

  df <- meta %>%
    filter(neighborhood.detailed == nhood) %>%
    group_by(neighborhood.detailed, cell_type.detailed) %>%
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
    ungroup() %>%
    mutate(neighborhood.detailed = nhood)
  return(secondLevel)
})
names(nhood.residents) <- neighborhoods


nhood.edgelists <- lapply(neighborhoods, function(nhood){

  count_matrix <- cellchat_list[[nhood]]@net$count
  count_graph <- graph_from_adjacency_matrix(count_matrix, weighted = TRUE, mode = "directed")
  count_edgelist <- igraph::as_data_frame(count_graph)
  colnames(count_edgelist) <- c("from", "to", "interaction_count")

  weight_matrix <- cellchat_list[[nhood]]@net$weight
  weight_graph <- graph_from_adjacency_matrix(weight_matrix, weighted = TRUE, mode = "directed")
  weight_edgelist <- igraph::as_data_frame(weight_graph)
  colnames(weight_edgelist) <- c("from", "to", "interaction_weight")

  edgelist <- inner_join(count_edgelist, weight_edgelist) %>%
    mutate(
      neighborhood.detailed = nhood,
      weighted_interaction_count = interaction_count * interaction_weight,
      count_zscore = (interaction_count - mean(interaction_count))/sd(interaction_count),
      weight_zscore = (interaction_weight - mean(interaction_weight))/sd(interaction_weight),
      interaction_zscore = (weighted_interaction_count - mean(weighted_interaction_count))/sd(weighted_interaction_count),
      weighted_interaction_count.scaled = as.numeric(scale(weighted_interaction_count, center = FALSE, scale = TRUE)),
      interaction_count.scaled = as.numeric(scale(interaction_count, center = FALSE, scale = TRUE)),
      interaction_weight.scaled = as.numeric(scale(interaction_weight, center = FALSE, scale = TRUE))
    ) %>%
    filter(from != to)

  return(edgelist)
})
names(nhood.edgelists) <- neighborhoods
qs::qsave(nhood.edgelists, file = 'data/merfish/BAYSOR/ligand_receptor/CellChat/neighborhood/neighborhood.cellchat.interaction_network.edgelist_dataframe.list_object.QS')

nhood.edgelist.combined <- do.call(rbind, nhood.edgelists)
qs::qsave(nhood.edgelist.combined, file = 'data/merfish/BAYSOR/ligand_receptor/CellChat/neighborhood/neighborhood.cellchat.interaction_network.edgelist_dataframe.combined_table.QS')
library(igraph)
library(tidygraph)
library(ggraph)
nhood.graphs <- lapply(neighborhoods, function(nhood){
  print(nhood)

  group.size <- data.frame(table(cellchat_list[[nhood]]@idents))
  colnames(group.size) <- c("cell_type.detailed", "no.cellchat_cells")

  node_data <- nhood.residents[[nhood]] %>%
    left_join(group.size) %>%
    mutate(name = cell_type.detailed,
           rank = as.numeric(factor(cell_type.detailed, levels = names(cell_type.detailed.colors)))) %>%
    arrange(rank) %>%
    data.frame()


  graph <- tbl_graph(nodes = node_data, edges = nhood.edgelists[[nhood]], node_key = 'name') %>%
    mutate(
      total_degree = centrality_degree(mode = "all"),
      out_degree = centrality_degree(mode = 'out'),
      in_degree = centrality_degree(mode = "in"),
      weight.hub_score = centrality_hub(weights = interaction_weight),
      count.hub_score = centrality_hub(weights = interaction_count),
      zscore.hub_score = centrality_hub(weights = interaction_zscore),
      weighted_count.hub_score = centrality_hub(weights = weighted_interaction_count),
      weight.authority_score = centrality_authority(weights = interaction_weight),
      count.authority_score = centrality_authority(weights = interaction_weight),
      weighted_count.authority_score = centrality_authority(weights = weighted_interaction_count),
      zscore.authority_score = centrality_authority(weights = interaction_zscore)
    )
  return(graph)
})
names(nhood.graphs) <- neighborhoods
qs::qsave(nhood.graphs, file = 'data/merfish/BAYSOR/ligand_receptor/CellChat/neighborhood/neighborhood.cellchat.interaction_network.tidygraph.list_object.QS')

set.seed(123)
nhood.network_layout.list <- lapply(neighborhoods, function(nhood){
  graph <- nhood.graphs[[nhood]] %>% filter(total_degree > 0)
  layout <- create_layout(graph, layout = "igraph", algorithm = 'nicely')
  return(layout)
})
names(nhood.network_layout.list) <- neighborhoods
qs::qsave(nhood.network_layout.list, file = 'data/merfish/BAYSOR/ligand_receptor/CellChat/neighborhood/neighborhood.cellchat.interaction_network.tidygraph_network_layout.list_object.QS')


nhood.chord_layout.list <- lapply(neighborhoods, function(nhood){
  graph <- nhood.graphs[[nhood]]
  layout <- create_layout(graph, layout = "linear", circular = TRUE)
  return(layout)
})
names(nhood.chord_layout.list) <- neighborhoods
qs::qsave(nhood.chord_layout.list, file = 'data/merfish/BAYSOR/ligand_receptor/CellChat/neighborhood/neighborhood.cellchat.interaction_network.tidygraph_chord_layout.list_object.QS')



library(Seurat)
library(Matrix)
library(tidyverse)
library(CellChat)
library(igraph)
library(tidygraph)
library(ggraph)
library(patchwork)
setwd("D:/Dropbox/_projects/NormalSkinAtlas")

detailed.cell_types <- names(cell_type.detailed.colors)

## ggplot themes
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

## load layouts
nhood.network_layout.list <- qs::qread('data/merfish/BAYSOR/ligand_receptor/CellChat/neighborhood/neighborhood.cellchat.interaction_network.tidygraph_network_layout.list_object.QS')
nhood.chord_layout.list <- qs::qread('data/merfish/BAYSOR/ligand_receptor/CellChat/neighborhood/neighborhood.cellchat.interaction_network.tidygraph_chord_layout.list_object.QS')
neighborhoods <- names(nhood.chord_layout.list)


## set seeds
set.seed(123)


## 1. network plots:
library(viridis)
### hub and authority scores
hub_network_list <- lapply(neighborhoods, function(nhood){
  layout <- nhood.network_layout.list[[nhood]]
  plt <- ggraph(layout) +
    geom_edge_fan(aes(color = interaction_zscore, edge_width = weighted_interaction_count.scaled),
                  alpha = 0.7, end_cap = circle(3, 'mm'),  start_cap = circle(3, 'mm'),
                  arrow = arrow(length = unit(3, 'mm')),  lineend = 'round', linejoin = 'round') +
    geom_node_point(aes(color = zscore.hub_score, size = zscore.hub_score)) +
    geom_node_text(aes(label = name), repel = TRUE, size = 5, point.padding = 0.25) +
    theme_classic() + spatial_theme + NoAxes() +
    scale_edge_color_gradientn(colours = c("gray90", "gray20")) +
    scale_color_viridis()  +
    labs(title = nhood)

  return(plt)
})

hub_network_panel <- patchwork::wrap_plots(hub_network_list, nrow=2)

svglite::svglite("figures/merfish/BAYSOR/ligand_receptor/CellChat/neighborhood/interaction_networks.hub_score.svg", height=16, width=36)
hub_network_panel
dev.off()

png("figures/merfish/BAYSOR/ligand_receptor/CellChat/neighborhood/interaction_networks.hub_score.png", height=16, width=36, units='in', res=300)
hub_network_panel
dev.off()


### hub and authority scores
authority_network_list <- lapply(neighborhoods, function(nhood){
  layout <- nhood.network_layout.list[[nhood]]
  plt <- ggraph(layout) +
    geom_edge_fan(aes(color = interaction_zscore, edge_width = weighted_interaction_count.scaled),
                  alpha = 0.7, end_cap = circle(3, 'mm'),  start_cap = circle(3, 'mm'),
                  arrow = arrow(length = unit(3, 'mm')),  lineend = 'round', linejoin = 'round') +
    geom_node_point(aes(color = zscore.authority_score, size = zscore.authority_score)) +
    geom_node_text(aes(label = name), repel = TRUE, size = 5, point.padding = 0.25) +
    theme_classic() + spatial_theme + NoAxes() +
    scale_edge_color_gradientn(colours = c("gray90", "gray20")) +
    scale_color_viridis() +
    labs(title = nhood)

  return(plt)
})

authority_network_panel <- patchwork::wrap_plots(authority_network_list, nrow=2)


svglite::svglite("figures/merfish/BAYSOR/ligand_receptor/CellChat/neighborhood/interaction_networks.authority_score.svg", height=16, width=36)
authority_network_panel
dev.off()

png("figures/merfish/BAYSOR/ligand_receptor/CellChat/neighborhood/interaction_networks.authority_score.png", height=16, width=36, units='in', res=300)
authority_network_panel
dev.off()




receivers_network_list <- lapply(neighborhoods, function(nhood){
  layout <- nhood.network_layout.list[[nhood]]
  plt <- ggraph(layout) +
    geom_edge_fan(aes(color = interaction_zscore, edge_width = weighted_interaction_count.scaled),
                  alpha = 0.7, end_cap = circle(3, 'mm'),  start_cap = circle(3, 'mm'),
                  arrow = arrow(length = unit(3, 'mm')),  lineend = 'round', linejoin = 'round') +
    geom_node_point(aes(color = in_degree, size = zscore.hub_score)) +
    geom_node_text(aes(label = name), repel = TRUE, size = 5, point.padding = 0.25) +
    theme_classic() + spatial_theme + NoAxes() +
    scale_edge_color_gradientn(colours = c("gray90", "gray20")) +
    scale_color_viridis() +
    labs(title = nhood)

  return(plt)
})

receivers_network_panel <- patchwork::wrap_plots(receivers_network_list, nrow=2)


svglite::svglite("figures/merfish/BAYSOR/ligand_receptor/CellChat/neighborhood/interaction_networks.receivers.svg", height=16, width=36)
receivers_network_panel
dev.off()

png("figures/merfish/BAYSOR/ligand_receptor/CellChat/neighborhood/interaction_networks.receivers.png", height=16, width=36, units='in', res=300)
receivers_network_panel
dev.off()


senders_network_list <- lapply(neighborhoods, function(nhood){
  layout <- nhood.network_layout.list[[nhood]]
  plt <- ggraph(layout) +
    geom_edge_fan(aes(color = interaction_zscore, edge_width = weighted_interaction_count.scaled),
                  alpha = 0.7, end_cap = circle(3, 'mm'),  start_cap = circle(3, 'mm'),
                  arrow = arrow(length = unit(3, 'mm')),  lineend = 'round', linejoin = 'round') +
    geom_node_point(aes(color = out_degree, size = zscore.hub_score)) +
    geom_node_text(aes(label = name), repel = TRUE, size = 5, point.padding = 0.25) +
    theme_classic() + spatial_theme + NoAxes() +
    scale_edge_color_gradientn(colours = c("gray90", "gray20")) +
    scale_color_viridis() +
    labs(title = nhood)

  return(plt)
})

senders_network_panel <- patchwork::wrap_plots(senders_network_list, nrow=2)
# senders_network_panel

svglite::svglite("figures/merfish/BAYSOR/ligand_receptor/CellChat/neighborhood/interaction_networks.senders.svg", height=164, width=36)
senders_network_panel
dev.off()

png("figures/merfish/BAYSOR/ligand_receptor/CellChat/neighborhood/interaction_networks.senders.png", height=14, width=36, units='in', res=300)
senders_network_panel
dev.off()



sender_identity_network_list <- lapply(neighborhoods, function(nhood){
  layout <- nhood.network_layout.list[[nhood]]

  node.index <- as.character(1:length(layout$name))
  names(node.index) <- as.character(layout$name)
  edge.colors <- cell_type.detailed.colors[names(node.index)]
  names(edge.colors) <- node.index

  plt <- ggraph(layout) +
    geom_edge_fan(aes(color = as.factor(from), edge_width = interaction_count),
                  alpha = 0.25, end_cap = circle(3, 'mm'),  start_cap = circle(3, 'mm'),
                  lineend = 'round', linejoin = 'round') +
    geom_node_point(aes(color = name, size = weight.hub_score)) +
    geom_node_text(aes(label = name, color=name), repel = TRUE, size = 6) +
    theme_classic() + spatial_theme + theme + NoAxes() +
    scale_size_continuous(limits = c(0, 1), breaks = c(0.00, 0.25, 0.50, 0.75, 1.00)) +
    scale_edge_width_continuous(limits = c(0, 35), breaks = c(0, 10, 20, 30)) +
    scale_edge_color_manual(values = edge.colors, guide = 'none') +
    scale_color_manual(values = cell_type.detailed.colors, guide = 'none') +
    labs(title = nhood)

  return(plt)
})


sender_identity_network_panel <- patchwork::wrap_plots(sender_identity_network_list, nrow=2,  guides = 'collect')
sender_identity_network_panel


svglite::svglite("figures/merfish/BAYSOR/ligand_receptor/CellChat/neighborhood/interaction_networks.sender_identity.svg", height=14, width=36)
sender_identity_network_panel
dev.off()

png("figures/merfish/BAYSOR/ligand_receptor/CellChat/neighborhood/interaction_networks.sender_identity.png", height=14, width=36, units='in', res=300)
sender_identity_network_panel
dev.off()



sender_identity_network_list.hublab <- lapply(neighborhoods, function(nhood){
  layout <- nhood.network_layout.list[[nhood]]

  node.index <- as.character(1:length(layout$name))
  names(node.index) <- as.character(layout$name)
  edge.colors <- cell_type.detailed.colors[names(node.index)]
  names(edge.colors) <- node.index

  top_hubs <- layout %>% top_n(10, weight.hub_score) %>% dplyr::select(name) %>% deframe() %>% as.character()
  layout$label <- ifelse(layout$name %in% top_hubs, as.character(layout$name), NA)

  plt <- ggraph(layout) +
    geom_edge_fan(aes(color = as.factor(from), edge_width = interaction_count),
                  alpha = 0.25, end_cap = circle(3, 'mm'),  start_cap = circle(3, 'mm'),
                  lineend = 'round', linejoin = 'round') +
    geom_node_point(aes(color = name, size = weight.hub_score)) +
    geom_node_text(aes(label = label, color=name), repel = TRUE, size = 6) +
    theme_classic() + spatial_theme+theme + NoAxes() +
    scale_size_continuous(limits = c(0, 1), breaks = c(0.00, 0.25, 0.50, 0.75, 1.00)) +
    scale_edge_width_continuous(limits = c(0, 35), breaks = c(0, 10, 20, 30)) +
    scale_edge_color_manual(values = edge.colors, guide = 'none') +
    scale_color_manual(values = cell_type.detailed.colors, guide = 'none') +
    labs(title = nhood)

  return(plt)
})


sender_identity_network_panel.hublab <- patchwork::wrap_plots(sender_identity_network_list.hublab, nrow=2,  guides = 'collect')
sender_identity_network_panel.hublab


svglite::svglite("figures/merfish/BAYSOR/ligand_receptor/CellChat/neighborhood/interaction_networks.sender_identity.top_hubs.labeled.svg", height=14, width=36)
sender_identity_network_panel.hublab
dev.off()

png("figures/merfish/BAYSOR/ligand_receptor/CellChat/neighborhood/interaction_networks.sender_identity.top_hubs.labeled.png", height=14, width=36, units='in', res=300)
sender_identity_network_panel.hublab
dev.off()



library(ggraph)
interaction_identity_chord_list <- lapply(neighborhoods, function(nhood){
  layout <- nhood.chord_layout.list[[nhood]]
  plt <- ggraph(layout) +
    geom_edge_arc(aes(color = interaction_zscore, edge_width = weighted_interaction_count.scaled),
                  alpha = 0.7, end_cap = circle(3, 'mm'),  start_cap = circle(3, 'mm'),
                  arrow = arrow(length = unit(3, 'mm')),  lineend = 'round', linejoin = 'round') +
    geom_node_point(aes(color = name, size = weighted_count.authority_score)) +
    geom_node_text(aes(label = name), repel = FALSE, size = 3, point.padding = 0.5) +
    theme_classic() + spatial_theme + NoAxes() +
    scale_edge_color_gradientn(colours = c("gray80", "black")) +
    scale_color_manual(values = cell_type.detailed.colors, guide = 'none') +
    labs(title = nhood) + coord_fixed()
  return(plt)
})


interaction_identity_chord_panel <- patchwork::wrap_plots(interaction_identity_chord_list, nrow=2) & NoLegend()


svglite::svglite("figures/merfish/BAYSOR/ligand_receptor/CellChat/neighborhood/interaction_identity.chord_panel.svg", height=16, width=36)
interaction_identity_chord_panel
dev.off()

png("figures/merfish/BAYSOR/ligand_receptor/CellChat/neighborhood/interaction_identity.chord_panel.png", height=16, width=36, units='in', res=300)
interaction_identity_chord_panel
dev.off()



sender_identity_chord_list <- lapply(neighborhoods, function(nhood){
  layout <- nhood.chord_layout.list[[nhood]]
  node.index <- as.character(1:length(layout$name))
  names(node.index) <- as.character(layout$name)
  edge.colors <- cell_type.detailed.colors[names(node.index)]
  names(edge.colors) <- node.index

  plt <- ggraph(layout) +
    geom_edge_arc(aes(color = as.factor(from), edge_width = weighted_interaction_count.scaled), alpha = 0.5,
                  end_cap = circle(3, 'mm'),  start_cap = circle(3, 'mm'),
                  lineend = 'round', linejoin = 'round') +
    geom_node_point(aes(color = name, size = weighted_count.authority_score)) +
    geom_node_text(aes(label = name), repel = FALSE, size = 3) +
    theme_classic() + spatial_theme + NoAxes() +
    scale_edge_color_manual(values = edge.colors, guide = 'none') +
    scale_color_manual(values = cell_type.detailed.colors, guide = 'none') +
    labs(title = nhood) + coord_fixed()
  return(plt)
})


sender_identity_chord_panel <- patchwork::wrap_plots(sender_identity_chord_list, nrow=1) & NoLegend()

svglite::svglite("figures/merfish/BAYSOR/ligand_receptor/CellChat/neighborhood/neighborhood_interaction_circos_plots.sender_identity_chord_plots.svg", height=5, width=36)
sender_identity_chord_panel
dev.off()

png("figures/merfish/BAYSOR/ligand_receptor/CellChat/neighborhood/neighborhood_interaction_circos_plots.sender_identity_chord_plots.png", height=16, width=36, units='in', res=300)
sender_identity_chord_panel
dev.off()


sender_identity_chord_list.nolab <- lapply(neighborhoods, function(nhood){
  layout <- nhood.chord_layout.list[[nhood]]
  node.index <- as.character(1:length(layout$name))
  names(node.index) <- as.character(layout$name)
  edge.colors <- cell_type.detailed.colors[names(node.index)]
  names(edge.colors) <- node.index

  plt <- ggraph(layout) +
    geom_edge_arc(aes(color = as.factor(from), edge_width = weighted_interaction_count.scaled), alpha = 0.5,
                  end_cap = circle(3, 'mm'),  start_cap = circle(3, 'mm'),
                  lineend = 'round', linejoin = 'round') +
    geom_node_point(aes(color = name, size = total_degree)) +
    theme_classic() + spatial_theme + NoAxes() +
    scale_edge_color_manual(values = edge.colors, guide = 'none') +
    scale_color_manual(values = cell_type.detailed.colors, guide = 'none') +
    labs(title = nhood) + coord_fixed()
  return(plt)
})


sender_identity_chord_panel.nolab <- patchwork::wrap_plots(sender_identity_chord_list.nolab, nrow=1) & NoLegend()
svglite::svglite("figures/merfish/BAYSOR/ligand_receptor/CellChat/neighborhood/neighborhood_interaction_circos_plots.sender_identity_chord_list.nolab.svg", height=5, width=36)
sender_identity_chord_panel.nolab
dev.off()

png("figures/merfish/BAYSOR/ligand_receptor/CellChat/neighborhood/neighborhood_interaction_circos_plots.sender_identity_chord_list.nolab..png", height=16, width=36, units='in', res=300)
sender_identity_chord_panel.nolab
dev.off()



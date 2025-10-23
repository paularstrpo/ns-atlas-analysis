# scRNA: Seurat Predicted IDs from MERFISH
## RUN CELL CHAT
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
library(CellChat)
setwd("D:/Dropbox/_projects/NormalSkinAtlas")

# restrict our cellchat analysis to the cell types that were analyzed for the MERFISH neighborhood cellchat
# this will effectively simulate the neighborhoods in the scRNA-seq for LR purposes.
# use the graph layouts of the merfish cellchat to get the list of cell types
merfish <- qs::qread("data/merfish/BAYSOR/seurat_objects/ns-atlas.merfish_baysor.scanvi_integrated.cellcharter.annotated.seurat_object.QS")
merfish <- subset(merfish,
              subset = cell_type.detailed != "Doublet" &
                cell_type.detailed != "Unknown" &
                tissue_compartment != "outside tissue")
merfish <- UpdateSeuratObject(merfish)
meta <- merfish@meta.data
rm(merfish); gc()

library(tidyverse)
library(DirichletReg)

# make a table of samples by clusters/cell-types -> number of cells per cluster per sample
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
neighborhoods <- levels(meta$neighborhood)
neighborhood.list <- apply(cell_prop, MARGIN = 2, function(x){names(x[x>=0.01])})
rm(drmatrix, x, cell_counts, meta); gc()
qs::qsave(neighborhood.list, file = 'data/merfish/BAYSOR/metadata/neighborhood.resident_celltypes.list.QS')

obj <- qs::qread("data/scrna/seurat_objects/normal_skin.scRNA.harmony.integrated.seurat_object.reclustered.annotated.filtered.QS")
genes <- read.csv("data/scrna/seurat_objects/normal_skin.scrna.harmony.integrated.reclustered.annotated.genes_filtered.gene_metadata.csv.gz")[,1]

full_thickness_studies <- c("Ji_Cell_2020", "Cheng_CellReports_2018", "Takahashi_JID_2020",
                            "VorstandlechnerNatComms2021", "YuImmunity2024", "Ganier_PNAS_2024",
                            "SoleBoleBolo_CommunBiol_2020", "Alkon_JAllergyClinImmunol2023",
                            "Rohajn_JAllergyClinImmunol_2020", "Zou_DevCell2021")
obj@meta.data$is_fullthickness <- obj@meta.data$study_id %in% full_thickness_studies

obj[['RNA']] <- JoinLayers(obj[['RNA']])
Idents(obj) <- "cell_type.detailed"
obj <- subset(obj, cell_type.detailed != "Doublet")
obj <- subset(obj, study_id %in% full_thickness_studies)
meta.full <- obj@meta.data[, c("cell_barcode", "sample_barcode", "cell_type.detailed")]
expr.full <- GetAssayData(obj[['RNA']], layer = 'data')
expr.full <- Matrix(expr.full[genes,], sparse = TRUE)
rm(obj); gc()

object.list <- list()
interaction_list <- list()
pathway_list <- list()

options(future.globals.maxSize=8e09)
CellChatDB <- CellChatDB.human # use CellChatDB.human if running on human data
# subset the expression data of signaling genes for saving computation cost
# exclude ECM and non-protein signaling since they will drown out other signal
CellChatDB.use <- subsetDB(CellChatDB, key = "annotation") # use Secreted Signaling

for(nhood in names(neighborhood.list)){
  print(nhood)
  neighborhood.celltypes <- intersect(unique(neighborhood.list[[nhood]]), unique(as.character(meta.full$cell_type.detailed)))

  meta <- meta.full[meta.full$cell_type.detailed %in% neighborhood.celltypes, c("cell_barcode", "sample_barcode", "cell_type.detailed")]
  colnames(meta) <- c("cell_barcode", "samples", "labels")
  meta$samples <- factor(meta$samples)
  meta$labels <- factor(meta$labels)

  cellchat <- createCellChat(object = expr.full[, meta$cell_barcode], meta = meta, group.by = "labels")
  gc()

  cellchat@DB <- CellChatDB.use
  cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database

  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(cellchat)
  cellchat <- filterCommunication(cellchat)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)


  if(length(cellchat@netP$pathways) > 0){
    cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "net")
    cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

    df.net <- data.frame(subsetCommunication(cellchat))
    df.net$neighborhood <- nhood
    interaction_list[[nhood]] <- df.net
    write.csv(df.net, file =
                paste0('data/scrna/ligand_receptor/neighborhoods/neighborhood_', make.names(nhood), '.cellchat_results.csv'))

    saveRDS(cellchat, file =
              paste0('data/scrna/ligand_receptor/neighborhoods/neighborhood_', make.names(nhood), '.cellchat_object.RDS'))

    object.list[[nhood]] <- cellchat



    df.net <- data.frame(subsetCommunication(cellchat, slot.name = 'netP'))
    df.net$neighborhood <- nhood
    pathway_list[[nhood]] <- df.net
    write.csv(df.net, file =
                paste0('data/scrna/ligand_receptor/neighborhoods/neighborhood_', make.names(nhood), '.cellchat_results.pathways.csv'))


    rm(df.net); invisible(gc())
  }
  rm(cellchat); invisible(gc())

}

qs::qsave(interaction_list, file = "data/scrna/ligand_receptor/neighborhoods/scrna.neighborhood.cellchat.result_tables.list.QS")
qs::qsave(pathway_list, file = "data/scrna/ligand_receptor/neighborhoods/scrna.neighborhood.cellchat.result_tables.pathways.list.QS")
qs::qsave(object.list, file = "data/scrna/ligand_receptor/neighborhoods/scrna.neighborhood.cellchat_object.list.QS")
gc()

cellchat <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = TRUE)
qs::qsave(cellchat, "data/scrna/ligand_receptor/neighborhoods/neighborhood.cellchat_object.merged.QS")

neighborhoods <- names(object.list)
neighborhood.colors <- pals::glasbey(10)

cellchat <- qs::qread("data/scrna/ligand_receptor/neighborhoods/neighborhood.cellchat_object.merged.QS")


df <- netVisual_bubble(cellchat, comparison = 1:length(object.list), angle.x = 45, remove.isolate = TRUE, return.data = TRUE)[[1]]
write.csv(df, file = "data/scrna/ligand_receptor/neighborhoods/neighborhood.net_ineractions.table.csv")


pathway_rank <- rankNet(cellchat,
                        mode = "comparison",
                        thresh=0.01,
                        comparison = 1:length(object.list),
                        measure = "weight",
                        sources.use = NULL,
                        targets.use = NULL,
                        do.stat = FALSE,
                        return.data = TRUE)[[1]][, 1:4]

write.csv(pathway_rank, file = "data/scrna/ligand_receptor/neighborhoods/neighborhood.pathway_ranking_table.thresh_0.1.csv")



pathway_rank <- rankNet(cellchat,
                        mode = "comparison",
                        thresh=0.01,
                        signaling.type = 'Secreted Signaling',
                        comparison = 1:length(object.list),
                        measure = "weight",
                        sources.use = NULL,
                        targets.use = NULL,
                        stacked = TRUE,
                        do.stat = FALSE,
                        return.data = TRUE)[[1]][, 1:4]
pathway_rank$signaling_type <- "Secreted Signaling"
write.csv(pathway_rank, file = "data/scrna/ligand_receptor/neighborhoods/neighborhood.pathway_ranking_table.thresh_0.1.secreted.csv")


pathway_rank <- rankNet(cellchat,
                        mode = "comparison",
                        thresh=0.01,
                        signaling.type = 'Cell-Cell Contact',
                        comparison = 1:length(object.list),
                        measure = "weight",
                        sources.use = NULL,
                        targets.use = NULL,
                        stacked = TRUE,
                        do.stat = FALSE,
                        return.data = TRUE)[[1]][, 1:4]
pathway_rank$signaling_type <- "Cell-Cell Contact"
write.csv(pathway_rank, file = "data/scrna/ligand_receptor/neighborhoods/neighborhood.pathway_ranking_table.thresh_0.1.contact.csv")

pathway_rank <- rankNet(cellchat,
                        mode = "comparison",
                        thresh=0.01,
                        signaling.type = 'ECM-Receptor',
                        comparison = 1:length(object.list),
                        measure = "weight",
                        sources.use = NULL,
                        targets.use = NULL,
                        stacked = TRUE,
                        do.stat = FALSE,
                        return.data = TRUE)[[1]][, 1:4]
pathway_rank$signaling_type <- "ECM-Receptor"
write.csv(pathway_rank, file = "data/scrna/ligand_receptor/neighborhoods/neighborhood.pathway_ranking_table.thresh_0.1.ecm.csv")

library(ComplexHeatmap)

for (nhood in names(object.list)){
  cc <- object.list[[nhood]]
  groupSize <- as.numeric(table(cc@idents))


  ht1 <- netAnalysis_signalingRole_heatmap(cc, pattern = "outgoing")
  ht2 <- netAnalysis_signalingRole_heatmap(cc, pattern = "incoming")

  pdf(file=paste0('figures/scrna/ligand_receptor/neighborhoods/cellchat_cell-cell_communication_pathway_heatmap_interaction_strength_', nhood, '.pdf'),width = 14,height = 12)
  draw(ht1 + ht2)
  dev.off()


  mat <- cc@net$weight
  pdf(file=paste0('figures/scrna/ligand_receptor/neighborhoods/cellchat_cell-cell_communication_heatmap_interaction_strength_', nhood, '.pdf'),width = 14,height = 12)
  pheatmap::pheatmap(mat,angle_col = 90)
  dev.off()

  pdf(file=paste0('figures/scrna/ligand_receptor/neighborhoods/cellchat_cell-cell_communicationNet_1x1_interaction_strength_', nhood, '.pdf'),width = 14,height = 12)
  for (i in 1:nrow(mat)) {
    mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
    mat2[i, ] <- mat[i, ]
    print(
      netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
    )
  }
  dev.off()
}


object.list <- qs::qread("data/scrna/ligand_receptor/neighborhoods/neighborhood.cellchat_object.list.QS")
cellchat <- qs::qread("data/scrna/ligand_receptor/neighborhoods/neighborhood.cellchat_object.merged.QS")


interactions.barplot <- compareInteractions(cellchat, group = names(object.list), color.use = neighborhood.colors) + coord_flip() + NoLegend() |
  compareInteractions(cellchat, group = names(object.list), measure = 'weight', color.use = neighborhood.colors) + coord_flip() + NoLegend()


svglite::svglite("figures/scrna/ligand_receptor/neighborhoods/number_interactions.barplot.svg", height=5, width=10)
interactions.barplot
dev.off()

png("figures/scrna/ligand_receptor/neighborhoods/number_interactions.barplot.png", height=5, width=10, units='in', res=300)
interactions.barplot
dev.off()


# interactions.barplot
top_pathways_per_neighborhod <- rankNet(cellchat, thresh = 0.01,
                                        mode = "comparison", measure = "weight",
                                        color.use = neighborhood.colors,
                                        sources.use = NULL, targets.use = NULL,
                                        stacked = TRUE, do.stat = FALSE,
                                        comparison = 1:length(object.list))


svglite::svglite("figures/scrna/ligand_receptor/neighborhoods/top_pathways_per_neighborhood.ranknet_propbarplot.svg", height=15, width=5)
top_pathways_per_neighborhod
dev.off()

png("figures/scrna/ligand_receptor/neighborhoods/top_pathways_per_neighborhood.ranknet_propbarplot.png", height=15, width=5, units='in', res=300)
top_pathways_per_neighborhod
dev.off()

# cell_type.detailed.colors <- readRDS("data/reference/cell_type.detailed.color_palette.final.RDS")
# cell_type.detailed.colors <- cell_type.detailed.colors[2:length(cell_type.detailed.colors)]
# detailed.cell_types <- names(cell_type.detailed.colors)

#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
num.link <- unlist(sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)}))
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax) # , color.use = cell_type.detailed.colors)
}

signaling_role.panel <- patchwork::wrap_plots(plots = gg, nrow=2)

svglite::svglite('figures/scrna/ligand_receptor/neighborhoods/signaling_role.per_celltype.scatterplot_panel.svg', height=10, width=24)
signaling_role.panel
dev.off()

png('figures/scrna/ligand_receptor/neighborhoods/signaling_role.per_celltype.scatterplot_panel.png', height=10, width=24, res=300, units='in')
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

png("figures/scrna/ligand_receptor/neighborhoods/neighborhood.signaling_role.all_patterns.heatmap.panel.png", height=8, width = 24, res = 150, units = 'in')
cowplot::plot_grid(plotlist = heatmap_list, nrow=2)
dev.off()


svglite::svglite("figures/scrna/ligand_receptor/neighborhoods/neighborhood.signaling_role.all_patterns.heatmap.panel.svg", height=12, width = 36)
cowplot::plot_grid(plotlist = heatmap_list, nrow=2)
dev.off()

heatmap_list <- lapply(names(object.list), function(site){
  ht <- netAnalysis_signalingRole_heatmap(object.list[[site]], pattern = "outgoing", signaling = pathway.union, title = site, width = 10, height=8)
  ht <- grid.grabExpr(draw(ht))
  return(ht)
})

png("figures/scrna/ligand_receptor/neighborhoods/neighborhood.signaling_role.outgoing_patterns.heatmap.panel.png", height=8, width = 24, res = 150, units = 'in')
cowplot::plot_grid(plotlist = heatmap_list, nrow=2)
dev.off()


svglite::svglite("figures/scrna/ligand_receptor/neighborhoods/neighborhood.outgoing_role.all_patterns.heatmap.panel.svg", height=8, width = 24)
cowplot::plot_grid(plotlist = heatmap_list, nrow=2)
dev.off()

heatmap_list <- lapply(names(object.list), function(site){
  ht <- netAnalysis_signalingRole_heatmap(object.list[[site]], pattern = "incoming", signaling = pathway.union, title = site, width = 10, height=8)
  ht <- grid.grabExpr(draw(ht))
  return(ht)
})

png("figures/scrna/ligand_receptor/neighborhoods/neighborhood.signaling_role.incoming_patterns.heatmap.panel.png", height=8, width = 24, res = 300, units = 'in')
cowplot::plot_grid(plotlist = heatmap_list, nrow=2)
dev.off()


svglite::svglite("figures/scrna/ligand_receptor/neighborhoods/neighborhood.incoming_role.all_patterns.heatmap.panel.svg", height=8, width = 24)
cowplot::plot_grid(plotlist = heatmap_list, nrow=2)
dev.off()


library(ComplexHeatmap)
#
merfish.cellchat.list <- qs::qread("data/merfish/BAYSOR/ligand_receptor/CellChat/neighborhood/neighborhood.cellchat_object.list.QS")
merfish.interactions.list <- lapply(names(merfish.cellchat.list), function(x){
  data.frame(subsetCommunication(merfish.cellchat.list[[x]])) %>%
    filter(pval < 0.01) %>%
    mutate(
      interacting_celltypes = paste0(source, "→", target),
      interaction_id = paste0(source, "→", target, ": ", interaction_name_2),
      neighborhood=x)
})
names(merfish.interactions.list) <- names(merfish.cellchat.list)
merfish.interactions.df <- bind_rows(merfish.interactions.list)
write.csv(merfish.interactions.df, file = "data/merfish/BAYSOR/ligand_receptor/CellChat/neighborhood/merfish.neighborhood.cellchat_interactions.concatenated.table.csv")

object.list <- qs::qread("data/scrna/ligand_receptor/neighborhoods/scrna.neighborhood.cellchat_object.list.QS")
#
## get LR hits that are unique to each neighborhood
names(object.list) <- names(merfish.cellchat.list)
scrna.interactions <- lapply(names(object.list), function(x){
  data.frame(subsetCommunication(object.list[[x]])) %>%
    filter(pval < 0.01) %>%
    mutate(
      interacting_celltypes = paste0(source, "→", target),
      interaction_id = paste0(source, "→", target, ": ", interaction_name_2),
      neighborhood=x,
      interaction.in_merfish =  paste0(source, "→", target, ": ", interaction_name_2) %in% unique(merfish.interactions.list[[x]]$interaction_id),
      cell_pair.in_merfish = paste0(source, "→", target) %in% unique(merfish.interactions.list[[x]]$interacting_celltypes)
      )
})
names(scrna.interactions) <- names(object.list)
# dfnet.filtered <- dfnet.combined %>% inner_join(merfish.interacting_celltypes)
neighborhoods <- names(object.list)
neighborhood_interactions <- lapply(scrna.interactions, function(x){
  x %>%
  dplyr::select(interaction_id) %>%
  unique() %>%
  deframe()
})
names(neighborhood_interactions) <- neighborhoods

library(ComplexHeatmap)
overlap_matrix <- list_to_matrix(neighborhood_interactions)
m <- make_comb_mat(overlap_matrix)
UpSet(m) # , comb_order = order(comb_size(combination_matrix)))

mfilt = m[comb_degree(m) == 1]
unique_interaction_codes =  comb_name(mfilt)
names(unique_interaction_codes) = rownames(mfilt)[rowSums(mfilt) != 0]
overlap_matrix <- data.frame(overlap_matrix)
codes <- apply(overlap_matrix, MARGIN = 1, paste0, collapse="")

interactiondf <- data.frame(interaction_id = names(codes[codes %in% unique_interaction_codes]),
                            code = codes[codes %in% unique_interaction_codes])

interactiondf$neighborhood <- factor(interactiondf$code, levels = unique_interaction_codes, labels = names(unique_interaction_codes))
scrna.interactions.df <- bind_rows(scrna.interactions) %>% mutate(unique_to_nhood=interaction_id %in% interactiondf$interaction_id)
write.csv(scrna.interactions.df, file = "data/scrna/ligand_receptor/neighborhoods/scrna.neighborhood.cellchat_interactions.concatenated.table.csv")
qs::qsave(scrna.interactions.df, "data/scrna/ligand_receptor/neighborhoods/neighborhood.cellchat_interactions.concatenated.table.QS")


scrna.interactions.filt <- scrna.interactions.df %>% filter(unique_to_nhood & cell_pair.in_merfish)
top_neighborhood_per_pathway <- pathway_rank %>% group_by(name) %>% top_n(1, wt =contribution) %>% ungroup()

table(scrna.interactions.filt$neighborhood)

pathways.show <- top_neighborhood_per_pathway %>%
  group_by(group) %>%
  top_n(15, contribution) %>%
  arrange(group) %>%
  ungroup() %>%
  dplyr::select(name) %>%
  deframe()


cellchat <- qs::qread("data/scrna/ligand_receptor/neighborhoods/neighborhood.cellchat_object.merged.QS")

neighborhood.colors <- pals::glasbey(n=10)
names(neighborhood.colors) <- make.names(neighborhoods)
top_pathways_per_neighborhod <- rankNet(cellchat, mode = "comparison", signaling = pathways.show, measure = "weight",color.use = neighborhood.colors,
        sources.use = NULL, targets.use = NULL, stacked = TRUE, do.stat = FALSE, comparison = 1:length(object.list))


svglite::svglite("figures/scrna/ligand_receptor/neighborhoods/top_pathways_per_neighborhood.ranknet_propbarplot.filtered.svg", height=15, width=5)
top_pathways_per_neighborhod
dev.off()

png("figures/scrna/ligand_receptor/neighborhoods/top_pathways_per_neighborhood.ranknet_propbarplot.filtered.png", height=15, width=5, units='in', res=300)
top_pathways_per_neighborhod
dev.off()


pathway_rank_mat <- pathway_rank %>%
  dplyr::select(group, name, contribution) %>%
  pivot_wider(names_from = "name", values_from="contribution") %>% data.frame()
rownames(pathway_rank_mat) <- pathway_rank_mat$group
pathway_rank_mat$group <- NULL

mat <- as.matrix(pathway_rank_mat)
mat <- apply(mat, MARGIN=2, function(x){ x / sum(x)})



cluster_pathways <- hclust(dist(t(mat)))
pathwayrank2 <- reshape2::melt(mat)
pathway_order = cluster_pathways$labels[cluster_pathways$order]
pathwayrank2$Var2 <- factor(pathwayrank2$Var2, levels=pathway_order)
library(ggdendro)
barplots <- ggplot(pathwayrank2) +
  aes(x=Var2, y = value, fill = Var1) +
  geom_bar(stat='identity', color='black') +
  scale_y_continuous(labels = scales::percent) +
  theme_classic() + rotate_x_text() +
  labs(x=NULL, y=NULL) +
  paletteer::scale_fill_paletteer_d("pals::glasbey")

hcdata <- dendro_data(cluster_pathways)
dend <- ggdendrogram(cluster_pathways)

svglite::svglite("figures/scrna/ligand_receptor/neighborhoods/pathways_per_neighborhood.propbarplot.clustered.svg", height=6, width=18)
(dend / barplots) + plot_layout(heights = c(0.5, 2), axes = 'collect_x')
dev.off()

panel.list <- list()
## per-neighborhood pathway rank plots
for (nhood in names(object.list)){
  cc <- object.list[[nhood]]
  pathways.show <- cc@netP$pathways

  # interaction counts
  sec <- rankNet(cc, mode = "single", measure = "count", thresh = 0.01, signaling.type = 'Secreted Signaling',
                 sources.use = NULL, targets.use = NULL) +
    labs(x='Secreted Signaling', y = "No. Interactions") +
    theme_classic()
  con <- rankNet(cc, mode = "single", measure = "count", thresh = 0.01, signaling.type = 'Cell-Cell Contact',
                 sources.use = NULL, targets.use = NULL) +
    labs(x='Cell-Cell Contact', y = "No. Interactions") +
    theme_classic()
  ecm <- rankNet(cc, mode = "single", measure = "count", thresh = 0.01, signaling.type = 'ECM-Receptor',
                 sources.use = NULL, targets.use = NULL) +
    labs(x='ECM-Receptor', y = "No. Interactions") +
    theme_classic()

  counts <- (sec / con / ecm) + plot_layout(guides='collect', axes='collect', heights = c(2,2,0.5), ncol=1)

  # interaction weights
  sec <- rankNet(cc, mode = "single", measure = "weight", thresh = 0.01, signaling.type = 'Secreted Signaling',
                        sources.use = NULL, targets.use = NULL) +
    labs(x='Secreted Signaling', y = "Interaction Weight") +
    theme_classic()
  con <- rankNet(cc, mode = "single", measure = "weight", thresh = 0.01, signaling.type = 'Cell-Cell Contact',
                 sources.use = NULL, targets.use = NULL) +
    labs(x='Cell-Cell Contact', y = "Interaction Weight") +
    theme_classic()
  ecm <- rankNet(cc, mode = "single", measure = "weight", thresh = 0.01, signaling.type = 'ECM-Receptor',
                 sources.use = NULL, targets.use = NULL) +
    labs(x='ECM-Receptor', y = "Interaction Weight") +
    theme_classic()

  weights <- (sec / con / ecm) + plot_layout(guides='collect', axes='collect', heights = c(2,2,0.5), ncol=1)


  panel <- (counts | weights) +
    plot_layout(guides="collect", axes='collect') +
    plot_annotation(title=nhood)

  svglite::svglite(paste0("figures/scrna/ligand_receptor/neighborhoods/", nhood, ".top_pathways.ranked_barplot.individual.svg"), height=10, width=5)
  print(panel)
  dev.off()

  png(paste0("figures/scrna/ligand_receptor/neighborhoods/", nhood, ".top_pathways.ranked_barplot.individual.png"), height=10, width=5, units='in', res=300)
  print(panel)
  dev.off()

  panel.list[[nhood]] <- panel
}

pdf("figures/scrna/ligand_receptor/neighborhoods/top_pathways.ranked_barplot.by_signalingtype.per_neighborhood.pdf", height=12, width=5)
print(panel.list)
dev.off()


panel.list <- list()
## per-neighborhood pathway rank plots
for (nhood in names(object.list)){
  cc <- object.list[[nhood]]
  pathways.show <- cc@netP$pathways

  # interaction counts
  counts <- rankNet(cc, mode = "single", measure = "count", thresh = 0.01,
                 sources.use = NULL, targets.use = NULL) +
    labs( y = "No. Interactions") +
    theme_classic()

  # interaction weights
  weights <- rankNet(cc, mode = "single", measure = "weight", thresh = 0.01,
                 sources.use = NULL, targets.use = NULL) +
    labs(y = "Interaction Weight") +
    theme_classic()

  panel <- (counts | weights) +
    plot_layout(guides="collect", axes='collect') +
    plot_annotation(title=nhood)

  svglite::svglite(paste0("figures/scrna/ligand_receptor/neighborhoods/", nhood, ".top_pathways.ranked_barplot.individual.summary.svg"), height=10, width=5)
  print(panel)
  dev.off()

  png(paste0("figures/scrna/ligand_receptor/neighborhoods/", nhood, ".top_pathways.ranked_barplot.individual.summary.png"), height=10, width=5, units='in', res=300)
  print(panel)
  dev.off()

  panel.list[[nhood]] <- panel
}

pdf("figures/scrna/ligand_receptor/neighborhoods/top_pathways.ranked_barplot.per_neighborhood.pdf", height=12, width=5)
print(panel.list)
dev.off()


pathway_summary <- data.frame()
## per-neighborhood pathway rank plots
for (nhood in names(object.list)){
  cc <- object.list[[nhood]]
  pathways.show <- cc@netP$pathways
  count_res <- rankNet(cc, mode = "single", measure = "count", thresh = 0.01,
                       sources.use = NULL, targets.use = NULL, return.data = TRUE)[[1]][, 1:2]
  count_res <- data.frame(neighborhood = nhood, pathway = count_res$name, interaction_count = count_res$contribution)

  weight_res <- rankNet(cc, mode = "single", measure = "weight", thresh = 0.01,
                        sources.use = NULL, targets.use = NULL, return.data = TRUE)[[1]][, 1:2]
  weight_res <- data.frame(neighborhood = nhood, pathway = weight_res$name, interaction_weight = weight_res$contribution)

  pathway_summary <- rbind(pathway_summary, data.frame(inner_join(count_res, weight_res)))
}


pathway_summary <- pathway_summary %>%
  group_by(pathway) %>%
  mutate(weighted_count = interaction_count * interaction_weight,
         relative_count = interaction_count / sum(interaction_count),
         relative_wcount = weighted_count / sum(weighted_count),
         relative_weight = interaction_weight / sum(interaction_weight)
         ) %>%
  ungroup()

pathway_variance <- pathway_summary %>%
  group_by(pathway) %>%
  summarise(var_wcount=var(weighted_count),
            var_weight=var(interaction_weight),
            var_count=var(interaction_count)) %>%
  ungroup() %>%
  arrange(-var_wcount, -var_weight, -var_count)

pathways.show <-pathway_summary %>%
  group_by(neighborhood) %>%
  top_n(15, wt=interaction_count) %>%
  ungroup() %>%
  dplyr::select(pathway) %>%
  deframe() %>%
  unique() %>%
  as.character() %>%
  sort()
pathways.show

ggplot(pathway_summary[pathway_summary$pathway %in% pathways.show,]) +
  aes(x=pathway, y = relative_wcount, fill = neighborhood) +
  geom_bar(stat='identity', color='black') +
  scale_y_continuous(labels = scales::percent) +
  theme_classic() + rotate_x_text() +
  labs(x=NULL, y=NULL) # +
  # scale_fill_manual(values = neighborhood.colors)

pathways.show <- pathway_summary %>%
  group_by(neighborhood) %>%
  top_n(15, wt=relative_wcount) %>%
  ungroup() %>%
  dplyr::select(pathway) %>%
  deframe() %>%
  unique() %>%
  as.character() %>%
  sort()
names(neighborhood.colors) <- names(object.list)
ggplot(pathway_summary[pathway_summary$pathway %in% pathways.show,]) +
  aes(x=pathway, y = relative_wcount, fill = neighborhood) +
  geom_bar(stat='identity', color='black') +
  scale_y_continuous(labels = scales::percent) +
  theme_classic() + rotate_x_text() +
  labs(x=NULL, y=NULL) +
  scale_fill_manual(values = neighborhood.colors)


ggplot(pathway_summary) +
  aes(x=pathway, y = neighborhood, color = relative_weight, size = relative_count) +
  geom_point(stat='identity') +
  theme_classic() + rotate_x_text() +
  scale_color_gradientn(colors = c('lightcyan', 'dodgerblue', "midnightblue"))

rel_wcount_matrix <- pathway_summary %>%
  dplyr::select(neighborhood, pathway, relative_wcount) %>%
  pivot_wider(names_from=neighborhood, values_from=relative_wcount, values_fill=0) %>%
  data.frame()
rownames(rel_wcount_matrix) <- rel_wcount_matrix$pathway
rel_wcount_matrix$pathway <- NULL
rel_wcount_matrix <- as.matrix(rel_wcount_matrix)

library(ComplexHeatmap)
col <- circlize::colorRamp2(c(0, 1), c('white', 'firebrick1'))
Heatmap(t(rel_wcount_matrix), col = col,
        width = unit(nrow(rel_wcount_matrix) * 3.5, "mm"),
        height = unit(ncol(rel_wcount_matrix) * 3.5, "mm"))



rel_count_matrix <- pathway_summary %>%
  dplyr::select(neighborhood, pathway, relative_count) %>%
  pivot_wider(names_from=neighborhood, values_from=relative_count, values_fill=0) %>%
  data.frame()
rownames(rel_count_matrix) <- rel_count_matrix$pathway
rel_count_matrix$pathway <- NULL
rel_count_matrix <- as.matrix(rel_count_matrix)

library(ComplexHeatmap)
col <- circlize::colorRamp2(c(0, 1), c('white', 'dodgerblue'))
Heatmap(t(rel_count_matrix), col = col,
        width = unit(nrow(rel_count_matrix) * 3.5, "mm"),
        height = unit(ncol(rel_count_matrix) * 3.5, "mm"))



secreted_signaling <- data.frame()
## per-neighborhood pathway rank plots
for (nhood in names(object.list)){
  cc <- object.list[[nhood]]
  pathways.show <- cc@netP$pathways
  count_res <- rankNet(cc, mode = "single", measure = "count", thresh = 0.01, signaling.type = 'Secreted Signaling',
                       sources.use = NULL, targets.use = NULL, return.data = TRUE)[[1]][, 1:2]
  count_res <- data.frame(neighborhood = nhood, pathway = count_res$name, interaction_count = count_res$contribution)

  weight_res <- rankNet(cc, mode = "single", measure = "weight", thresh = 0.01, signaling.type = 'Secreted Signaling',
                        sources.use = NULL, targets.use = NULL, return.data = TRUE)[[1]][, 1:2]
  weight_res <- data.frame(neighborhood = nhood, pathway = weight_res$name, interaction_weight = weight_res$contribution)

  secreted_signaling <- rbind(secreted_signaling, data.frame(inner_join(count_res, weight_res)))
}



secreted_signaling <- secreted_signaling %>%
  group_by(pathway) %>%
  mutate(weighted_count = interaction_count * interaction_weight,
         relative_count = interaction_count / sum(interaction_count),
         relative_wcount = weighted_count / sum(weighted_count),
         relative_weight = interaction_weight / sum(interaction_weight)
  ) %>%
  ungroup()

ggplot(secreted_signaling) +
  aes(x=pathway, y = relative_count, fill = neighborhood) +
  geom_bar(stat='identity', color='black') +
  scale_y_continuous(labels = scales::percent) +
  theme_classic() +
  labs(x=NULL, y=NULL) +
  scale_fill_manual(values = neighborhood.colors) +
  coord_flip()

ggplot(secreted_signaling) +
  aes(x=pathway, y = relative_wcount, fill = neighborhood) +
  geom_bar(stat='identity', color='black') +
  scale_y_continuous(labels = scales::percent) +
  theme_classic() +
  labs(x=NULL, y=NULL) +
  scale_fill_manual(values = neighborhood.colors) +
  coord_flip()

ggplot(secreted_signaling) +
  aes(x=reorder(pathway, interaction_count), y = interaction_count) +
  geom_bar(stat='identity', color='black') +
  facet_wrap(~neighborhood, scales='free', nrow=2) +
  theme_classic() +
  labs(x=NULL, y=NULL) +
  coord_flip()

normalize <- function(x) {
  return((x- min(x)) /(max(x)-min(x)))
}

panel.list <- list()
## per-neighborhood pathway rank plots
for (nhood in names(object.list)){
  cc <- object.list[[nhood]]
  pathways.show <- secreted_signaling %>%
    filter(neighborhood==nhood) %>%
    top_n(15, wt=interaction_count)

  panel.list[[nhood]] <- ggplot(pathways.show) +
    aes(x = interaction_count, y = reorder(pathway, interaction_count)) +
    geom_bar(stat='identity', fill=neighborhood.colors[nhood]) +
    theme_classic() +
    labs(title=nhood, x="No. Interactions", y = NULL)
}

svglite::svglite("figures/scrna/ligand_receptor/neighborhoods/secreted_signaling.interaction_count.barplot.perneighborhood.panel.svg", height=3, width=18)
wrap_plots(panel.list, nrow=1)
dev.off()


png("figures/scrna/ligand_receptor/neighborhoods/secreted_signaling.interaction_count.barplot.perneighborhood.panel.png", height=3, width=18, units='in', res=300)
wrap_plots(panel.list, nrow=1)
dev.off()

nhood='PERIVASC-I'
pathways.show <- secreted_signaling %>%
  filter(neighborhood==nhood) %>%
  top_n(15, wt=interaction_count) %>%
  dplyr::select(pathway) %>%
  deframe() %>%
  as.character()

scrna.interactions.df <- qs::qread("data/scrna/ligand_receptor/neighborhoods/neighborhood.cellchat_interactions.concatenated.table.QS")
# top_senders
# top_receivers
pv1 <- scrna.interactions.df %>%
  filter(neighborhood=='PERIVASC-I') %>%
  filter(cell_pair.in_merfish) %>%
  filter(pathway_name %in% pathways.show) %>%
  filter(source %in% top_senders & target %in% top_receivers)
library(igraph)
library(ggraph)
library(tidygraph)
merfish.graph_layouts <- qs::qread("data/merfish/BAYSOR/ligand_receptor/CellChat/neighborhood/neighborhood.cellchat.interaction_network.tidygraph.list_object.QS")

graph <- merfish.graph_layouts[['PERIVASC-I']]
hub_scores = colSums(object.list[['PERIVASC-I']]@net$count) + rowSums(object.list[['PERIVASC-I']]@net$count)

netVisual_bubble(object.list[['PERIVASC-I']],
                 pairLR.use = data.frame(interaction_name=unique(pv1$interaction_name)),
                 sources.use = unique(as.character(pv1$source)),
                 targets.use = unique(as.character(pv1$target)),
                 thresh = 0.01,
                 remove.isolate = TRUE)

pdf("figures/scrna/ligand_receptor/neighborhoods/pathway_genes.perivasc.I.circos.pdf")
for (p in pathways.show){
  print(
    netVisual_chord_gene(object.list[['PERIVASC-I']], signaling = p, thresh = 0.01, legend.pos.x = 100, title.name = p)
  )
}
dev.off()

netVisual_chord_cell(object.list[['PERIVASC-I']],
                     signaling = 'TNF',
                     thresh = 0.01)

netVisual_chord_gene(object.list[['PERIVASC-I']], signaling = pathways.show,
                     sources=c("CD8+ Tc", "CD4+ Th", "VEC", "HEC"),
                     targets = c("VEC", "Peri", "Perivasc Fib II", "Perivasc Fib I"),
                     thresh = 0.01, legend.pos.x = 100)


nhood <- "PERIVASC-I"
cc <- object.list[[nhood]]
senders <- sort(colSums(cc2@net$weight), decreasing=TRUE)
receivers <- sort(rowSums(cc2@net$weight), decreasing=TRUE)

pathways.show <- secreted_signaling %>%
  filter(neighborhood==nhood) %>%
  top_n(15, wt=interaction_count) %>%
  dplyr::select(pathway) %>%
  deframe() %>%
  as.character()

netVisual_chord_gene(cc,
                     signaling = pathways.show,
                     sources = names(senders[1:3]),
                     targets = names(receivers[1:3]),
                     thresh = 0.001, legend.pos.x = 100)

#
#
# netVisual_chord_gene(object.list[['PERIVASC-I']], signaling = pathways.show,
#                      targets = c("CD8+ Tc", "CD4+ Th", "VEC", "HEC", "Naïve T"),
#                      sources = c("VEC", "Peri", "Perivasc Fib II", "Perivasc Fib I", "Retic Fib I"),
#                      thresh = 0.001, legend.pos.x = 100)

## tangram - LR scores - neighborhood LR variationa cross sites

library(tidyverse)
library(Matrix)
library(Seurat)
library(reticulate)
library(rhdf5)
library(anndata)
library(SingleCellExperiment)
library(scater)
library(edgeR)
library(limma)
library(variancePartition)
library(dreamlet)
library(variancePartition)
library(sva)


setwd("D:/Dropbox/_projects/NormalSkinAtlas")

merfish.obj <- qs::qread("data/merfish/BAYSOR/seurat_objects/ns-atlas.merfish_baysor.scanvi_integrated.cellcharter.annotated.filtered.seurat_object.QS")


# Add observation metadata
meta <- merfish.obj@meta.data
h5ad <- "data/merfish/BAYSOR/ligand_receptor/Tangram/merfish.tangram_projected.lr_scores.full.anndata_object.h5ad"
adata <- anndata::read_h5ad(h5ad)
expr <- adata$layers['gmean']
expr <- t(expr)
gc()


barcodes <- intersect(rownames(meta), colnames(expr))
full_info <- meta[barcodes,]
full_expr <- expr[, barcodes]
#rm(adata, expr, meta); gc()


full_info$neighborhood.cell_type <- paste0(make.names(full_info$neighborhood), "_", make.names(full_info$cell_type.detailed))
full_info$observation_id <- make.names(paste0(full_info$sample_barcode, ".", full_info$neighborhood,  ".", full_info$cell_type.detailed))



neighborhood.list <- qs::qread('data/merfish/BAYSOR/metadata/neighborhood.resident_celltypes.list.QS')
neighborhood_celltypes <- unlist(
  sapply(names(neighborhood.list), function(x){
    paste0(make.names(x),'_', make.names(neighborhood.list[[x]]))
  })
)


full_info <- full_info[full_info$neighborhood.cell_type %in% neighborhood_celltypes, ]
full_info <- full_info[full_info$donor_id %in% c("D165","D077", "D082", "D107", "D145", "D149", "D151"), ]
barcodes <- intersect(rownames(full_info), colnames(full_expr))
full_info <- full_info[barcodes,]
full_expr <- full_expr[, barcodes]

full_obs <- unique(full_info$observation_id)
full_lrmat <- sapply(full_obs, function(sn){
  rowMeans(full_expr[, full_info$observation_id == sn, drop = FALSE], na.rm = TRUE)
})

full_info <- unique(full_info[, c("observation_id", "sample_barcode", "cell_type.detailed", "neighborhood", 'anatomic_site', 'donor_id', 'donor_sex', 'batch', 'gene_panel')])
rownames(full_info) <- full_info$observation_id
full_info <- full_info[colnames(full_lrmat),]
full_lrmat <- as.matrix(full_lrmat[,rownames(full_info)])

write.csv(full_lrmat, file = "data/merfish/BAYSOR/ligand_receptor/Tangram/knn.lr_scores.tangram_imputed.sample-nhood-celltype.pseudobulked.csv")
## adjust for donor ID
full_lrmat_adj <- sva::ComBat(full_lrmat, batch = as.factor(as.character(full_info$batch)))
write.csv(full_lrmat_adj, file = "data/merfish/BAYSOR/ligand_receptor/Tangram/knn.lr_scores.tangram_imputed.sample-nhood-celltype.pseudobulked.batch_corrected.csv")

vpform <- ~ (1|neighborhood) + (1|cell_type.detailed)+(1|anatomic_site)+(1|donor_id)+(1|donor_sex)+(1|batch)

varPart <- fitExtractVarPartModel(full_lrmat_adj, formula = vpform, full_info)
vpres <- sortCols(varPart)
plotVarPart(vpres)

pca.obj <- prcomp(t(as.matrix(full_lrmat_adj)), center=TRUE, scale=TRUE)
pcdf <- pca.obj$x
full_info$PC1 <- pcdf[, 1]
full_info$PC2 <- pcdf[, 2]


(ggplot(full_info) + aes(x = PC1, y = PC2, color=batch) + geom_point() + theme_bw()) +
  (ggplot(full_info) + aes(x = PC1, y = PC2, color=anatomic_site) + geom_point() + theme_bw()) +
  (ggplot(full_info) + aes(x = PC1, y = PC2, color=gene_panel) + geom_point() + theme_bw()) +
  (ggplot(full_info) + aes(x = PC1, y = PC2, color=neighborhood) + geom_point() + theme_bw())

library(ggrepel)

umap <- uwot::umap(t(as.matrix(full_lrmat_adj)))
full_info$umap_1 <- umap[, 1]
full_info$umap_2 <- umap[, 2]
(ggplot(full_info) + aes(x = umap_1, y = umap_2, color=batch) + geom_point() + theme_bw()) +
  (ggplot(full_info) + aes(x = umap_1, y = umap_2, color=anatomic_site) + geom_point() + theme_bw()) +
  (ggplot(full_info) + aes(x = umap_1, y = umap_2, color=cell_type.detailed) + geom_point() + theme_bw()) +
  (ggplot(full_info) + aes(x = umap_1, y = umap_2, color=neighborhood) + geom_point() + theme_bw())

result.list <- list()
neighborhoods <- levels(unique(full_info$neighborhood))
for (nhood in neighborhoods) {
  print(nhood)
  info <- full_info[full_info$neighborhood == nhood,]
  lrmat <- full_lrmat_adj[, rownames(info)]

  info$anatomic_site <- factor(make.names(as.character(info$anatomic_site)))
  info$batch <- factor(make.names(info$batch))

  sites <- unique(info$anatomic_site)

  ### across anatomic sites only
  form <-  ~ 0 + anatomic_site + (1|cell_type.detailed) + (1|donor_sex) + (1|batch)
  contrast_list <- lapply(sites, function(site){

    side1 <- paste0("anatomic_site", site)
    side2 <- sites[sites != site]
    len <- length(side2)

    side2 <- paste0("(", paste(paste0("anatomic_site", side2), collapse=" + "), ") / ", len)

    contr <- paste0(side1, " - ", side2)

  })
  names(contrast_list) <- sites
  contrasts <- unlist(contrast_list)

  L <- makeContrastsDream(form, data = info, contrasts = contrasts)
  fit <- dream(exprObj = lrmat, form, info, L)
  # attr(., 'errors')

  fit <- eBayes(fit)

  site.df <- lapply(1:length(sites), function(x){
    df <- data.frame(topTable(fit, number = Inf, coef = x, confint = TRUE))
    df$neighborhood <- nhood
    df$anatomic_site <- sites[x]
    df$LR <- rownames(df)
    df$is.signif <- df$adj.P.Val < 0.05
    df$direction <- ifelse(df$logFC < 0, "Down", "Up")
    return(df)
  }) %>%
    bind_rows()

  result.list[[nhood]] <- site.df
  result.list
  rm(lrmat, info, contrasts, contrast_list, site.df, fit, L); gc()

}

gc()

resultdf <- bind_rows(result.list)
write.csv(resultdf, file = 'data/merfish/BAYSOR/ligand_receptor/Tangram/tangram.imputed_lr.knn_scores.per-neighborhood.anatomic_site-comparison.dLR_results.csv')




library(ggpubr)
library(ggrepel)
library(tidyverse)
library(ComplexHeatmap)
resultdf_full <- read.csv('data/merfish/BAYSOR/ligand_receptor/Tangram/tangram.imputed_lr.knn_scores.per-neighborhood.anatomic_site-comparison.dLR_results.csv')
neighborhoods <- unique(resultdf_full$neighborhood)

for (nhood in neighborhoods){

  resultdf <- resultdf_full %>% filter(neighborhood == nhood)

  top_lrs <- resultdf %>%
    filter(is.signif & direction == "Up") %>%
    group_by(anatomic_site) %>%
    top_n(10, wt=logFC)

  pmat <- resultdf[resultdf$LR %in% top_lrs$LR, ] %>%
    mutate(pstar = gsub(" ", "", gsub('\\.', '', gtools::stars.pval(adj.P.Val)))) %>%
    dplyr::select(LR, anatomic_site, pstar) %>%
    unique() %>%
    pivot_wider(names_from='anatomic_site', values_from='pstar') %>%
    data.frame()
  rownames(pmat) <- pmat$LR
  pmat$LR <- NULL

  pmat <- as.matrix(pmat)


  lfc <- resultdf[resultdf$LR %in% top_lrs$LR, ] %>%
    dplyr::select(LR, anatomic_site, z.std) %>%
    unique() %>%
    pivot_wider(names_from='anatomic_site', values_from='z.std') %>%
    data.frame()
  rownames(lfc) <- lfc$LR
  lfc$LR <- NULL


  svglite::svglite(paste0("figures/merfish/BAYSOR/ligand_receptor/",nhood,".tangram_knn.site_lr_heatmap.svg"), height=25, width=25)
  print({
    Heatmap(as.matrix(lfc), name='Z-Score', width = ncol(lfc)*unit(5, "mm"), height = nrow(lfc)*unit(5, "mm"),
          cell_fun = function(j, i, x, y, width, height, fill) {
            grid.text(sprintf("%s", pmat[i, j]), x, y, gp = gpar(fontsize = 10))
          })
  })
  dev.off()
}



library(tidyverse)
library(Matrix)
library(Seurat)
library(SeuratObject)

run_prefix <- "ns-atlas.merfish"
data_path <- "data/merfish/"

setwd("D:/Dropbox/_projects/NormalSkinAtlas/")
# setwd("~/Dropbox/_projects/NormalSkinAtlas")


## ggplot themes
theme <- theme(
  legend.position = 'right',
  text = element_text(color = 'black', size = 10),
  axis.text = element_text(color = 'black', size = 10),
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



common_genes <- read.csv('data/reference/NormalSkin.CrossPanelOverlap.GeneList.txt')[,1]
obj <- qs::qread("data/merfish/BAYSOR/seurat_objects/ns-atlas.merfish_baysor.scanvi_integrated.cellcharter.annotated.seurat_object.QS")

count_dist_by_batch <- VlnPlot(obj, features = c("nCount_RNA", 'nFeature_RNA'),
        pt.size = 0, raster=FALSE, group.by = 'batch', stack = TRUE, log = TRUE) +
  geom_boxplot(width=0.2, aes(fill=NULL), outliers = FALSE, show.legend=FALSE) +
  theme_bw() + theme + labs(x=NULL, y = "MERFISH Imaging Batch")+ NoLegend()

count_dist_by_batch

svglite::svglite("figures/merfish/BAYSOR/quality_control/post_filtering.qc_distributions.violin_plots.by_batch.svg", height=4, width=6)
print(count_dist_by_batch)
dev.off()


median(obj$nCount_RNA)
mean(obj$nCount_RNA)
sd(obj$nCount_RNA)

median(obj$nFeature_RNA)
mean(obj$nFeature_RNA)
sd(obj$nFeature_RNA)

median(obj$volume)
mean(obj$volume)

median(obj$transcript_density)
mean(obj$transcript_density)


exprs <- t(as.matrix(GetAssayData(obj, layer='counts')))

cts <- data.frame(t(rowsum(exprs, group = obj@meta.data$gene_panel)))
cts$gene <- rownames(cts)

cts <- cts[common_genes,]

cor.test(cts$NSv1, cts$NSv2, method = 'pearson', exact = TRUE)
cor.test(cts$NSv1, cts$NSv2, method = 'spearman', exact = TRUE)


panel.corr.counts_in_cells <- ggplot(cts) + aes(x = NSv1, y = NSv2) +
  geom_point() + geom_smooth(method='lm', color='red', se=FALSE) +
  theme_bw() + theme + scale_x_log10() + scale_y_log10() + coord_fixed() +
  labs(title = "Gene Panel Correlation: Transcript counts within Segmented Cells")
panel.corr.counts_in_cells

counts_per_gene.in_cells = colSums(exprs)


library(tidyverse)
library(jsonlite)
library(vroom)
setwd("D:/Dropbox/_projects/NormalSkinAtlas/")






tx_files <- paste0("D:/MERFISH/baysor/", summary_df$run_name, "/", summary_df$region_name, "/baysor_detected_transcripts.csv")
names(tx_files) <- summary_df$sample_id


tx_list <- lapply(names(tx_files), function(sample){
  vroom(tx_files[[sample]], num_threads = 12, progress = TRUE) %>%
    filter(transcript_id != "-1") %>%
    mutate(sample_id = sample) %>%
    data.frame()
})

transcript_df <- tx_list %>% bind_rows()
write.csv(transcript_df, file='data/merfish/BAYSOR/qc/detected_transcripts.all_samples.csv')


fovdf <- transcript_df %>%
  group_by(sample_id, fov) %>%
  summarise(transcript_count = n()) %>%
  ungroup() %>%
  inner_join(summary_df)
write.csv(fovdf, file = 'data/merfish/BAYSOR/qc/fov_metrics.csv')

pseudobulk <- transcript_df %>%
  inner_join(summary_df[, c("sample_id", "sample_barcode")]) %>%
  group_by(sample_barcode, gene) %>%
  summarise(transcript_count = n()) %>%
  ungroup() %>%
  dplyr::select(sample_barcode, gene, transcript_count) %>%
  pivot_wider(names_from = gene, values_from = transcript_count, values_fill = 0) %>%
  data.frame()
rownames(pseudobulk) <- pseudobulk$sample_barcode
pseudobulk$sample_barcode <- NULL
gc()
write.csv(pseudobulk, file = "data/merfish/BAYSOR/qc/pseudobulk_from_txdf.counts_matrix.csv")

gene_info <- read.csv('data/merfish/BAYSOR/seurat_objects/ns-atlas.merfish_baysor.scvi_integrated.gene_data.csv', row.names = 1)
pseudobulk <- read.csv("data/merfish/BAYSOR/qc/pseudobulk_from_txdf.counts_matrix.csv", row.names = 1)

pbmat <- t(pseudobulk)
summarydf <- read.csv("data/merfish/BAYSOR/qc/baysor_sample_qc_metrics_annotated.csv")
annodf <- summarydf[, c('sample_barcode', 'anatomic_site', 'collection_source', 'donor_id', 'gene_panel')]
rownames(annodf) <- annodf$sample_barcode
annodf$sample_barcode <- NULL
library(ComplexHeatmap)

site_colors <- setNames(Polychrome::light.colors(n=length(unique(annodf$anatomic_site))), nm=unique(annodf$anatomic_site))
donor_colors <- setNames(Polychrome::glasbey.colors(n=length(unique(annodf$donor_id))), nm=unique(annodf$donor_id))
collection_colors <- c("autopsy"='hotpink', 'mohs'='goldenrod1', 'panniculectomy'='dodgerblue')
panel_colors <- c("NSv1"="blue", "NSv2"='orange')
cmat.spearman <- cor(pbmat, method='spearman')

hm_anno <- HeatmapAnnotation(df = annodf[rownames(cmat.spearman),],
                             col = list(
                               anatomic_site=site_colors,
                               collection_source = collection_colors,
                               donor_id = donor_colors,
                               gene_panel=panel_colors
                               )
                             )

row_anno <- rowAnnotation(df = annodf[rownames(cmat.spearman),],
                          col = list(
                            anatomic_site=site_colors,
                            collection_source = collection_colors,
                            donor_id = donor_colors,
                            gene_panel=panel_colors
                          )
)



svglite::svglite('figures/merfish/BAYSOR/quality_control/sample_tx_correlation.svg', width = 30, height = 30)
ComplexHeatmap::Heatmap(cmat.spearman, name = 'Spearman',
                        col = c('white', 'firebrick1'),
                        bottom_annotation = hm_anno,
                        left_annotation = row_anno,
                        heatmap_width = unit(5*ncol(cmat.spearman), "mm"),
                        heatmap_height = unit(5*nrow(cmat.spearman), "mm"))
dev.off()


svglite::svglite('figures/merfish/BAYSOR/quality_control/sample_tx_correlation.blue2red.svg', width = 30, height=30)
ComplexHeatmap::Heatmap(cmat.spearman, name = 'Spearman',
                        bottom_annotation = hm_anno,
                        left_annotation = row_anno,
                        heatmap_width = unit(5*ncol(cmat.spearman), "mm"),
                        heatmap_height = unit(5*nrow(cmat.spearman), "mm"))
dev.off()

pseudobulk$sample_barcode <- rownames(pseudobulk)

cts <- transcript_df %>%
  filter(gene %in% common_genes) %>%
  inner_join(summary_df[, c("sample_id", "sample_barcode")]) %>%
  group_by(sample_barcode, gene) %>%
  summarise(transcript_count = n()) %>%
  ungroup() %>%
  inner_join(summary_df[, c("sample_barcode", 'gene_panel')]) %>%
  group_by(gene, gene_panel) %>%
  summarise(sum_transcripts = sum(transcript_count),
            mean_transcripts = mean(transcript_count)) %>%
  ungroup()


sumdf <- cts %>%
  dplyr::select(gene, gene_panel, sum_transcripts) %>%
  pivot_wider(names_from = 'gene_panel', values_from = 'sum_transcripts')

cor.test(log1p(sumdf$NSv1), log1p(sumdf$NSv2), method = 'pearson', exact = TRUE)
cor.test(sumdf$NSv1, sumdf$NSv2, method = 'spearman', exact = TRUE)

panel.corr.sum <- ggplot(sumdf) +
  aes(x = log1p(NSv1), y = log1p(NSv2)) +
  geom_point() +
  geom_smooth(method='lm', color='red', se=TRUE, alpha=0.25) +
  theme_bw() + theme +
  labs(title = "MERFISH Gene Panel Correlation",
       x = "Gene Panel 1: Log1p(Counts)",
       y='Gene Panel 2: Log1p(Counts)')
panel.corr.sum

svglite::svglite("figures/merfish/BAYSOR/quality_control/gene_panel_txdf_correlation.svg", height=5, width=5)
panel.corr.sum
dev.off()

svglite::svglite("figures/merfish/BAYSOR/quality_control/transcripts_per_fov.per_sample.svg", height=15, width=5)
ggplot(fovdf) +
  aes(y = sample_barcode, x = transcript_count) +
  geom_boxplot() + theme_bw() + theme +
  scale_x_log10() +
  labs(x='Tx/FOV')
dev.off()

svglite::svglite("figures/merfish/BAYSOR/quality_control/transcripts_per_fov.per_sample.log1p.svg", height=15, width=5)
ggplot(fovdf) +
  aes(y = sample_barcode, x = log1p(transcript_count)) +
  geom_boxplot() + theme_bw() + theme +
  labs(x='Log1p(Tx/FOV)')
dev.off()

svglite::svglite("figures/merfish/BAYSOR/quality_control/transcripts_per_fov.global_histogram.svg", height=5, width=7)
ggplot(fovdf) +
  aes(x = transcript_count) +
  geom_histogram(bins=100, color='black', fill='gray70') +
  theme_bw() + theme +
  labs(x = "Trancripts per FOV", y = "Count") +
  scale_x_log10()
dev.off()


svglite::svglite("figures/merfish/BAYSOR/quality_control/transcripts_per_fov.global_histogram.log1p.svg", height=5, width=7)
ggplot(fovdf) +
  aes(x = log1p(transcript_count)) +
  geom_histogram(bins=100, color='black', fill='gray70') +
  theme_bw() + theme +
  labs(x = "log1p(Tx/FOV)", y = "Count")
dev.off()

mean(fovdf$transcript_count)
sd(fovdf$transcript_count)
median(fovdf$transcript_count)

fovdf <- read.csv('data/merfish/BAYSOR/qc/fov_metrics.csv')
tx_per_fov <- fovdf %>%
  group_by(batch, run_name) %>%
  summarise(tx_per_fov = mean(transcript_count),
            total_tx_count = sum(transcript_count)) %>%
  ungroup()

summary_df <- readxl::read_xlsx('data/metadata/sample_metadata.xlsx') %>%
  dplyr::select(run_name, DV200, pmi_hrs, gene_panel) %>%
  unique() %>%
  inner_join(tx_per_fov)

ffpe <- summary_df[!is.na(summary_df$DV200),]

cor.test(ffpe$DV200, ffpe$tx_per_fov, method = "pearson", exact = TRUE)
cor.test(ffpe$DV200, ffpe$tx_per_fov, method = 'spearman', exact = TRUE)

dv200_vs_txperfov <- ggplot(summary_df[!is.na(summary_df$DV200),]) +
  aes(x = log1p(tx_per_fov), y = DV200) +
  geom_smooth(method='lm', color='red') +
  geom_point() +
  theme_bw() + theme +
  labs(x= "Log1p(Tx/FOV)", y = "DV200",
       title=paste0("Spearman: ", round(cor.test(ffpe$DV200, ffpe$tx_per_fov, method = 'spearman', exact = TRUE)$estimate, digits=2)))
dv200_vs_txperfov
svglite::svglite("figures/annotated/dv200_vs_txperfov_corr.svg", height = 5, width = 5)
dv200_vs_txperfov
dev.off()

cor.test(summary_df$DV200, summary_df$tx_per_fov, na.action = 'omit', method='pearson')
cor.test(summary_df$DV200, summary_df$tx_per_fov, na.action = 'omit', method='spearman')

library(ggpubr)
#
# svglite::svglite("figures/Apr2024/MERFISH/quality_control/anatomicsite_vs_txperfov_box_raw.svg", height = 5, width = 7)
# ggplot(summary_df) +
#   aes(y = anatomic_site, x = log1p(tx_per_fov), fill = anatomic_site) +
#   geom_boxplot(show.legend=FALSE) + theme_bw()
# dev.off()

autopsy <- summary_df[!is.na(summary_df$pmi_hrs),]

cor.test(autopsy$pmi_hrs, autopsy$tx_per_fov, method='pearson')
cor.test(autopsy$pmi_hrs, autopsy$tx_per_fov, method='spearman')
txperfov_pmi <- ggplot(summary_df[!is.na(summary_df$pmi_hrs),]) +
  aes( x = log1p(tx_per_fov), y = pmi_hrs) +
  geom_smooth(method='lm', color='red') +
  geom_point() +
  theme_bw() + theme +
  labs(x= "Log1p(Tx/FOV)", y = "PMI (hrs)",
       title=paste0("Spearman: ", round(cor.test(autopsy$pmi_hrs, autopsy$tx_per_fov, method='spearman')$estimate, digits=2)))

svglite::svglite("figures/annotated/pmi_vs_txperfov_corr.svg", height = 5, width = 5)
txperfov_pmi
dev.off()

library(patchwork)
(dv200_vs_txperfov / txperfov_pmi)

obj <- qs::qread("data/merfish/seurat_objects/ns-atlas.projected_integration.processed.seurat_object.annotated.QS")
obj <- UpdateSeuratObject(obj)
metadf <- obj@meta.data

donor_metadata <- metadf %>% dplyr::select(donor_id, donor_age, donor_sex, collection_source, collection_type) %>% unique()


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



## -
svglite::svglite("figures/merfish/donor_sex_distribution.svg", height = 5, width = 7)
ggplot(donor_metadata) + aes(x = donor_sex, fill = donor_sex) + facet_wrap(~collection_type) +
  theme_classic() + geom_bar(color='black', show.legend = FALSE) +
  scale_fill_manual(values = c("Male" = "skyblue2", "Female" = "hotpink")) + theme
dev.off()

ggplot(donor_metadata) + aes(x = donor_age, fill = collection_source) + geom_density() + geom_rug() + theme_classic() +  facet_wrap(~collection_source)


svglite::svglite("figures/merfish/donor_age_distribution.svg", height = 5, width = 7)
ggplot(donor_metadata) + aes(x =collection_type, y = donor_age, fill = collection_type) +
  geom_violin(show.legend=FALSE) + geom_boxplot(aes(fill=NULL),width=0.2, show.legend = FALSE) + theme_classic() +
  facet_wrap(~collection_type, scales ='free_x') + theme +
  scale_fill_manual(values = c("autopsy"= "royalblue", "surgical discards" = "darkorchid3"))
dev.off()



abd <- run_summary_df %>% filter(anatomic_site == "abdomen")



## AH21 technical replicate corrleation
v1.genes <- read.csv("data/reference/FinalGeneList_NormSkinNew_VA00227.csv")$Vizgen.Gene

repdf <- data.frame(t(pseudobulk[c("D018_SKIN_NS_S02_R01","D018_SKIN_NS_S02_R02"), intersect(v1.genes, colnames(pseudobulk))]))
colnames(repdf) <- c("R01", "R02")
repdf$Gene <- rownames(repdf)
repdf <- repdf %>% filter(R01 > 0 & R02>0)
cor.test(repdf$R01, repdf$R02, method='pearson')
cor.test(repdf$R01, repdf$R02, method='spearman')
d018_replicates <- ggplot(repdf) +
  aes(x=log1p(R01), y = log1p(R02)) +
  labs(title = "D018 Technical Replicates",
       x='D018 Replicate 1: Log1p(Counts)',
       y='D018 Replicate 2: Log1p(Counts)') +
  geom_point() +
  geom_smooth(method='lm', color='red', alpha=0.25) +
  theme_bw() + theme

svglite::svglite("figures/merfish/BAYSOR/quality_control/D018_AH21_TechnicalReplicateCorrelation.svg", height = 5, width = 6)
d018_replicates
dev.off()


panel.corr.sum/d018_replicates

## GTex Correlation
lowerleg_tpm <- read.table(gzfile("data/bulkRNA/GTEx/gene_tpm_2017-06-05_v8_skin_sun_exposed_lower_leg.gct.gz"),
                           fill = TRUE, skip=2, header = TRUE) %>%
  dplyr::select(-id)

rownames(lowerleg_tpm) <- lowerleg_tpm$Name
genes <- str_split(lowerleg_tpm$Name, pattern = "\\.", simplify = TRUE)[,1]
lowerleg_tpm$Description <- NULL
lowerleg_tpm$Name <- NULL
lowerleg_tpm <- rowsum(lowerleg_tpm, group = genes, reorder = TRUE)
lowerleg_tpm.means <- rowMeans(lowerleg_tpm)
lowerleg_tpm.medians <- matrixStats::rowMedians(as.matrix(lowerleg_tpm))

lowerleg <- data.frame(Ensemble.ID = names(lowerleg_tpm.means),
                       GTEx_LowerLeg_TPM_Mean = log1p(lowerleg_tpm.means),
                       GTEx_LowerLeg_TPM_Median = log1p(lowerleg_tpm.medians))

suprapubic_tpm <- read.table(gzfile("data/bulkRNA/GTEx/gene_tpm_2017-06-05_v8_skin_not_sun_exposed_suprapubic.gct.gz"),
                           fill = TRUE, skip=2, header = TRUE) %>%
  dplyr::select(-id)

rownames(suprapubic_tpm) <- suprapubic_tpm$Name
genes <- str_split(suprapubic_tpm$Name, pattern = "\\.", simplify = TRUE)[,1]
suprapubic_tpm$Name <- NULL
suprapubic_tpm$Description <- NULL
suprapubic_tpm <- rowsum(suprapubic_tpm, group = genes, reorder = TRUE)
suprapubic_tpm.means <- rowMeans(suprapubic_tpm)
suprapubic_tpm.medians <- matrixStats::rowMedians(as.matrix(suprapubic_tpm))

suprapubic <- data.frame(Ensemble.ID = names(suprapubic_tpm.means),
                         GTEx_Suprapubic_TPM_Mean = log1p(suprapubic_tpm.means),
                         GTEx_Suprapubic_TPM_Median= log1p(suprapubic_tpm.medians))

gtex.tpm <- inner_join(lowerleg, suprapubic)

all.genes <- colnames(merfish)
overlapping.genes <- intersect(suprapubic$gene, all.genes)


merfish <- read.csv("data/merfish/BAYSOR/qc/pseudobulk_from_txdf.counts_matrix.csv", row.names = 1)
sample_meta <- read.csv("data/merfish/BAYSOR/qc/baysor_sample_qc_metrics_annotated.csv", row.names = 1)
rownames(sample_meta) <- sample_meta$sample_barcode
merfish <- merfish[rownames(sample_meta),]
sample_barcodes <- rownames(sample_meta)


v1_panel <- read.csv("data/reference/FinalGeneList_NormSkinNew_VA00227.csv")
v1_panel$logFPKM <- log1p(v1_panel$Abundance)
v2_panel <- read.csv("data/reference/FinalGeneList_NormalSkinV2_CP1251.csv")
v2_panel$logFPKM <- log1p(v2_panel$Abundance)

v1_counts <- data.frame(t(merfish[sample_meta$gene_panel == "NSv1",]))
v1_mean <- rowMeans(v1_counts)
v1_median <-  matrixStats::rowMedians(as.matrix(v1_counts))
v1_sum <- rowSums(v1_counts)
v1_samples <- colnames(v1_counts)
v1_counts$Vizgen.Gene <- rownames(v1_counts)
v1_counts$Median_MERFISH_Expr <- log1p(v1_median)
v1_counts$Mean_MERFISH_Expr <- log1p(v1_mean)
v1_counts$Total_MERFISH_Counts <- v1_sum

v1df <- v1_panel %>% inner_join(gtex.tpm) %>% inner_join(v1_counts)

v2_counts <- data.frame(t(merfish[sample_meta$gene_panel == "NSv2",]))
v2_mean <- rowMeans(v2_counts)
v2_median <- matrixStats::rowMedians(as.matrix(v2_counts))
v2_sum <- rowSums(v2_counts)

v2_samples <- colnames(v2_counts)
v2_counts$Vizgen.Gene <- rownames(v2_counts)

v2_counts$Median_MERFISH_Expr <- log1p(v2_median)
v2_counts$Mean_MERFISH_Expr <- log1p(v2_mean)
v2_counts$Total_MERFISH_Counts <- v2_sum

v2df <- v2_panel %>% inner_join(gtex.tpm) %>% inner_join(v2_counts)

v1cor <- lapply(v1_samples, function(x){

  df <- v1df
  df[, x] <- log1p(df[,x])
  se.mean.pearson <- cor.test(df[,x], df$GTEx_LowerLeg_TPM_Mean, method = 'pearson') %>% broom::tidy()
  se.median.pearson <- cor.test(df[,x], df$GTEx_LowerLeg_TPM_Median, method = 'pearson') %>% broom::tidy()

  se.mean.spearman <- cor.test(df[,x], df$GTEx_LowerLeg_TPM_Mean, method = 'spearman') %>% broom::tidy()
  se.median.spearman <- cor.test(df[,x], df$GTEx_LowerLeg_TPM_Median, method = 'spearman') %>% broom::tidy()


  se.fpkm.spearman <- cor.test(df[,x], df$logFPKM, method = 'spearman') %>% broom::tidy()
  se.fpkm.pearson <- cor.test(df[,x], df$logFPKM, method = 'pearson') %>% broom::tidy()


  ne.mean.pearson <- cor.test(df[,x], df$GTEx_Suprapubic_TPM_Mean, method = 'pearson') %>% broom::tidy()
  ne.median.pearson <- cor.test(df[,x], df$GTEx_Suprapubic_TPM_Median, method = 'pearson') %>% broom::tidy()

  ne.mean.spearman <- cor.test(df[,x], df$GTEx_Suprapubic_TPM_Mean, method = 'spearman') %>% broom::tidy()
  ne.median.spearman <- cor.test(df[,x], df$GTEx_Suprapubic_TPM_Median, method = 'spearman') %>% broom::tidy()


  res <- data.frame(sample_barcode = x,
                    gene_panel = "NSv1",

                    GTex.SunExposed_LowerLeg.FPKM.Pearson.corr = se.fpkm.pearson$estimate,
                    GTex.SunExposed_LowerLeg.FPKM.Pearson.pval = se.fpkm.pearson$p.value,

                    GTex.SunExposed_LowerLeg.FPKM.Spearman.corr = se.fpkm.spearman$estimate,
                    GTex.SunExposed_LowerLeg.FPKM.Spearman.pval = se.fpkm.spearman$p.value,

                    GTex.SunExposed_LowerLeg.Mean_TPM.Pearson.corr = se.mean.pearson$estimate,
                    GTex.SunExposed_LowerLeg.Mean_TPM.Pearson.pval = se.mean.pearson$p.value,

                    GTex.SunExposed_LowerLeg.Mean_TPM.Spearman.corr = se.mean.spearman$estimate,
                    GTex.SunExposed_LowerLeg.Mean_TPM.Spearman.pval = se.mean.spearman$p.value,

                    GTex.SunExposed_LowerLeg.Median_TPM.Pearson.corr = se.median.pearson$estimate,
                    GTex.SunExposed_LowerLeg.Median_TPM.Pearson.pval = se.median.pearson$p.value,

                    GTex.SunExposed_LowerLeg.Median_TPM.Spearman.corr = se.median.spearman$estimate,
                    GTex.SunExposed_LowerLeg.Median_TPM.Spearman.pval = se.median.spearman$p.value,

                    GTex.NonExposed_Suprapubic.Mean_TPM.Pearson.corr = ne.mean.pearson$estimate,
                    GTex.NonExposed_Suprapubic.Mean_TPM.Pearson.pval = ne.mean.pearson$p.value,

                    GTex.NonExposed_Suprapubic.Mean_TPM.Spearman.corr = ne.mean.spearman$estimate,
                    GTex.NonExposed_Suprapubic.Mean_TPM.Spearman.pval = ne.mean.spearman$p.value,

                    GTex.NonExposed_Suprapubic.Median_TPM.Pearson.corr = ne.median.pearson$estimate,
                    GTex.NonExposed_Suprapubic.Median_TPM.Pearson.pval = ne.median.pearson$p.value,

                    GTex.NonExposed_Suprapubic.Median_TPM.Spearman.corr = ne.median.spearman$estimate,
                    GTex.NonExposed_Suprapubic.Median_TPM.Spearman.pval = ne.median.spearman$p.value
  )

  return(res)

}) %>% bind_rows()



v2cor <- lapply(v2_samples, function(x){

  df <- v2df
  df[, x] <- log1p(df[,x])
  se.mean.pearson <- cor.test(df[,x], df$GTEx_LowerLeg_TPM_Mean, method = 'pearson') %>% broom::tidy()
  se.median.pearson <- cor.test(df[,x], df$GTEx_LowerLeg_TPM_Median, method = 'pearson') %>% broom::tidy()

  se.mean.spearman <- cor.test(df[,x], df$GTEx_LowerLeg_TPM_Mean, method = 'spearman') %>% broom::tidy()
  se.median.spearman <- cor.test(df[,x], df$GTEx_LowerLeg_TPM_Median, method = 'spearman') %>% broom::tidy()


  se.fpkm.spearman <- cor.test(df[,x], df$logFPKM, method = 'spearman') %>% broom::tidy()
  se.fpkm.pearson <- cor.test(df[,x], df$logFPKM, method = 'pearson') %>% broom::tidy()


  ne.mean.pearson <- cor.test(df[,x], df$GTEx_Suprapubic_TPM_Mean, method = 'pearson') %>% broom::tidy()
  ne.median.pearson <- cor.test(df[,x], df$GTEx_Suprapubic_TPM_Median, method = 'pearson') %>% broom::tidy()

  ne.mean.spearman <- cor.test(df[,x], df$GTEx_Suprapubic_TPM_Mean, method = 'spearman') %>% broom::tidy()
  ne.median.spearman <- cor.test(df[,x], df$GTEx_Suprapubic_TPM_Median, method = 'spearman') %>% broom::tidy()


  res <- data.frame(sample_barcode = x,
                    gene_panel = "NSv2",

                    GTex.SunExposed_LowerLeg.FPKM.Pearson.corr = se.fpkm.pearson$estimate,
                    GTex.SunExposed_LowerLeg.FPKM.Pearson.pval = se.fpkm.pearson$p.value,

                    GTex.SunExposed_LowerLeg.FPKM.Spearman.corr = se.fpkm.spearman$estimate,
                    GTex.SunExposed_LowerLeg.FPKM.Spearman.pval = se.fpkm.spearman$p.value,


                    GTex.SunExposed_LowerLeg.Mean_TPM.Pearson.corr = se.mean.pearson$estimate,
                    GTex.SunExposed_LowerLeg.Mean_TPM.Pearson.pval = se.mean.pearson$p.value,

                    GTex.SunExposed_LowerLeg.Mean_TPM.Spearman.corr = se.mean.spearman$estimate,
                    GTex.SunExposed_LowerLeg.Mean_TPM.Spearman.pval = se.mean.spearman$p.value,

                    GTex.SunExposed_LowerLeg.Median_TPM.Pearson.corr = se.median.pearson$estimate,
                    GTex.SunExposed_LowerLeg.Median_TPM.Pearson.pval = se.median.pearson$p.value,

                    GTex.SunExposed_LowerLeg.Median_TPM.Spearman.corr = se.median.spearman$estimate,
                    GTex.SunExposed_LowerLeg.Median_TPM.Spearman.pval = se.median.spearman$p.value,

                    GTex.NonExposed_Suprapubic.Mean_TPM.Pearson.corr = ne.mean.pearson$estimate,
                    GTex.NonExposed_Suprapubic.Mean_TPM.Pearson.pval = ne.mean.pearson$p.value,

                    GTex.NonExposed_Suprapubic.Mean_TPM.Spearman.corr = ne.mean.spearman$estimate,
                    GTex.NonExposed_Suprapubic.Mean_TPM.Spearman.pval = ne.mean.spearman$p.value,

                    GTex.NonExposed_Suprapubic.Median_TPM.Pearson.corr = ne.median.pearson$estimate,
                    GTex.NonExposed_Suprapubic.Median_TPM.Pearson.pval = ne.median.pearson$p.value,

                    GTex.NonExposed_Suprapubic.Median_TPM.Spearman.corr = ne.median.spearman$estimate,
                    GTex.NonExposed_Suprapubic.Median_TPM.Spearman.pval = ne.median.spearman$p.value
  )

  return(res)

}) %>% bind_rows()



cordf <- bind_rows(v1cor, v2cor)
cor.cols <-  colnames(cordf)[grepl("corr", colnames(cordf))]
corlong <- cordf[, c("sample_barcode", "gene_panel", cor.cols)] %>% pivot_longer(cor.cols, names_to = 'name', values_to = "correlation")


test_info <- str_split(corlong$name, "\\.", simplify=TRUE)
corlong$GTex.Tissue_Type <- test_info[,2]
corlong$GTex.Expr_Metric <- test_info[,3]
corlong$Correlation_Test <- test_info[,4]


gtex_cor <- ggplot(corlong) +
  aes(x = GTex.Tissue_Type, y = correlation, fill = GTex.Expr_Metric) +
  facet_wrap(~Correlation_Test) + geom_boxplot() + theme_bw()


final.cordf <- corlong %>% filter(GTex.Expr_Metric == "Median_TPM")

svglite::svglite("figures/merfish/BAYSOR/quality_control/gtex_median_tpm_correlation.boxplot_per_sample.svg", height=5, width=7)
ggplot(final.cordf) + aes(x =  GTex.Tissue_Type, y = correlation, color = GTex.Tissue_Type) +
  theme_bw() + facet_wrap(~Correlation_Test) +
  ggbeeswarm::geom_quasirandom(show.legend=FALSE) +
  geom_boxplot(aes(color=NULL), alpha=0, show.legend=FALSE)
dev.off()



ggplot(v2df) + aes(x = Median_MERFISH_Expr, y = GTEx_LowerLeg_TPM_Median) + geom_point() + theme_bw() +
  geom_smooth(method='lm')

ggplot(v2df) + aes(x = Median_MERFISH_Expr, y = GTEx_Suprapubic_TPM_Median) + geom_point() + theme_bw() +
  geom_smooth(method='lm')


ggplot(v1df) + aes(x = Median_MERFISH_Expr, y = GTEx_LowerLeg_TPM_Median) + geom_point() + theme_bw() +
  geom_smooth(method='lm')


ggplot(v1df) + aes(x = Mean_MERFISH_Expr, y = GTEx_LowerLeg_TPM_Median) + geom_point() + theme_bw() +
  geom_smooth(method='lm')

cols <- c("Vizgen.Gene", "Gene_Panel", "GTEx_LowerLeg_TPM_Median", "GTEx_Suprapubic_TPM_Median", "Median_MERFISH_Expr", "logFPKM")
v1df$Gene_Panel <- "NSv1"
v2df$Gene_Panel <- "NSv2"
global.corvalues <- rbind(v1df[,cols], v2df[, cols]) %>% filter(Median_MERFISH_Expr > 0)


lowerleg_global_scatter <- ggplot(global.corvalues) +
  aes(x = Median_MERFISH_Expr, y = GTEx_LowerLeg_TPM_Median) +
  geom_point() + geom_smooth(method='lm', se=FALSE, color='red') +
  facet_wrap(~Gene_Panel) + theme_bw() + theme +
  labs(x='MERFISH Median Log1p(Counts)', y= 'Lower Leg', title='Correlation with Bulk RNA-seq')

suprapubic_global_scatter <- ggplot(global.corvalues) +
  aes(x = Median_MERFISH_Expr, y = GTEx_Suprapubic_TPM_Median) +
  geom_point() + geom_smooth(method='lm', se=FALSE, color='red') +
  facet_wrap(~Gene_Panel) + theme_bw() + theme +
  labs(x='MERFISH Median Log1p(Counts)', y= 'Suprapubic')

library(patchwork)

svglite::svglite("figures/merfish/BAYSOR/quality_control/gtex_median_tpm_correlation.scatter_global.svg", height=10, width=10)
(lowerleg_global_scatter / suprapubic_global_scatter) + plot_layout(axes = 'collect')
dev.off()



library(ggbeeswarm)

seg_benchmark <- ggplot(segdf) +
  aes(x = segmentation, y = tx_pct.within_cells, color = segmentation) +
  geom_quasirandom(show.legend = FALSE) +
  geom_boxplot(width=0.2, aes(color=NULL),alpha=0, show.legend = FALSE) +
  theme_bw() + theme + NoLegend() +
  scale_color_manual(values = c("dodgerblue", "goldenrod1", "orchid2")) +
  labs(x = "Segmentation Method", y = "Transcripts within a Cell (%)")

SF1_SCATTER <- ((dv200_vs_txperfov / txperfov_pmi) | (panel.corr.sum/d018_replicates) | (lowerleg_global_scatter / suprapubic_global_scatter) | (seg_benchmark)) +
  plot_layout(axes = 'collect', widths = c(1,1,2, 1))

svglite::svglite("figures/merfish/BAYSOR/quality_control/SF1_PANEL.svg", height=8, width=20)
print(SF1_SCATTER)
dev.off()

### ---- ###
library(ComplexHeatmap)
# Jaccard similarity index function
jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}


library(Matrix)
## scRNA vs MERFISH marker qc
merfish <- qs::qread("data/merfish/BAYSOR/seurat_objects/ns-atlas.merfish_baysor.scanvi_intergated.cellcharter.seurat_object.QS") %>%
  subset(subset = cell_type.detailed != "Doublet") %>%
  UpdateSeuratObject()
merfish <- UpdateSeuratObject(merfish)

merfish[['RNA']]$data <- as.matrix(merfish[['RNA']]$data)
scrna <- qs::qread("data/scrna/seurat_objects/normal_skin.scRNA.harmony.integrated.seurat_object.reclustered.annotated.filtered.QS") %>%
  UpdateSeuratObject()

genes <- intersect(common_genes, rownames(scrna))

## Only test on features in both MERFISH and scRNA
Idents(merfish) <- "cell_type.broad"
merfish.markerdf <- FindAllMarkers(merfish, assay = 'RNA',
                                   only.pos=TRUE,
                                   logfc.threshold = 0,
                                   return.thresh = 1,
                                   features = genes)

Idents(scrna) <- "cell_type.broad"
scrna.markerdf <- FindAllMarkers(scrna, assay = 'RNA',
                                 only.pos=TRUE,
                                 return.thresh = 1,
                                 logfc.threshold = 0,
                                 features = genes)

# rm(scrna, merfish); gc()
merfish.celltypes <- c("Bas KC", "Spn KC", "Grn KC",
                       'Ecc Duct', "Ecc Gland","Seb", "HFE",
                       "Fibro","Peri", "EC", "LEC", "Schwann",
                       "Plasma", "Lym", "LC", "Mac/DC", "Mast")
scrna.celltypes <- c("Cyc KC", "Bas KC", "Spn KC", "Grn KC",
                     "Ecc Gland","Melano", "Bas HFE", "Diff HFE",
                     "Fib", "Peri", "EC", "LEC",
                     "Plasma", "CD4+ T Lym", "CD8+ T Lym",
                     "LC", "DC", "Mono", "Mac", "Mast")
# make cell type marker lists for each assay
merfish.markerlist <- split(merfish.markerdf, f = merfish.markerdf$cluster) %>%
  lapply(FUN = function(x){
    x %>% filter(p_val_adj < 0.05) %>% dplyr::select(gene) %>% unique() %>% deframe()
  })
merfish.celltypes <- names(merfish.markerlist)

scrna.markerlist <- split(scrna.markerdf, f = scrna.markerdf$cluster) %>%
  lapply(FUN = function(x){
    x %>% filter(p_val_adj < 0.05) %>% dplyr::select(gene) %>% unique() %>% deframe()
  })


scrna.markerlist <- scrna.markerlist[scrna.celltypes]
merfish.markerlist <- merfish.markerlist[merfish.celltypes]
jaccard.matrix <- sapply(merfish.markerlist, function(m){
  sapply(scrna.markerlist, function(s){jaccard(m, s)})
})

col_fun = circlize::colorRamp2(c(0,0.2, 0.4, 0.6, 0.8, 1), c("white",'salmon', 'tomato1', "red2", 'red3', "red4"))
svglite::svglite("figures/annotated/merfish2scrna.global.marker_jaccard.heatmap.svg", height=5, width=6)
Heatmap(jaccard.matrix,
        col = col_fun,
        name='Jaccard Index',
        column_title = "MERFISH",
        row_title = "scRNA",
        cluster_rows = FALSE,
        cluster_columns = FALSE)
dev.off()

genes <- intersect(common_genes, rownames(scrna))
scrna.counts <- AggregateExpression(scrna, features = genes, group.by = "cell_type.broad")[[1]]
scrna.counts <- as.matrix(scrna.counts[genes, ])

merfish.counts <- AggregateExpression(merfish, features = genes, group.by = "cell_type.broad")[[1]]
merfish.counts <- as.matrix(merfish.counts[genes,])


spearman.df <- lapply(merfish.celltypes, function(m){
  mdat <- log1p(merfish.counts[, m])
  lapply(names(scrna.markerlist), function(s){
    gg <- union(scrna.markerlist[[s]], merfish.markerlist[[m]])
    sdat <- log1p(scrna.counts[, s])

    cor.test(mdat[gg], sdat[gg], method = 'spearman') %>%
      broom::tidy() %>%
      mutate(scrna.cell_type = s,
             merfish.cell_type = m,
             num_genes = length(gg),
             jaccard_index = jaccard(scrna.markerlist[[s]], merfish.markerlist[[m]]))
  }) %>% bind_rows()
}) %>%
  bind_rows() %>%
  mutate(p.adj = p.adjust(p.value, method = "BH"))

cormat <- spearman.df %>% dplyr::select(scrna.cell_type, merfish.cell_type, estimate) %>% pivot_wider(names_from = "merfish.cell_type", values_from = "estimate") %>% data.frame()
rownames(cormat) <- cormat$scrna.cell_type
cormat$scrna.cell_type <- NULL
colnames(cormat) <- merfish.celltypes
smat = cormat
pmat <- spearman.df %>% dplyr::select(scrna.cell_type, merfish.cell_type, p.adj) %>% pivot_wider(names_from = "merfish.cell_type", values_from = "p.adj") %>% data.frame()
rownames(pmat) <- pmat$scrna.cell_type
pmat$scrna.cell_type <- NULL
colnames(pmat) <- merfish.celltypes
pmat <- as.matrix(pmat)
pstar <- gtools::stars.pval(pmat)
pstar[pstar == "."] <- " "
pstar <- as.matrix(pstar)
colnames(pstar) <- merfish.celltypes

#
# col_fun = circlize::colorRamp2(c(0, max(smat)), c("white", "red"))
# Heatmap(as.matrix(smat),
#         name = 'Spearman', col = col_fun,
#         cluster_rows = FALSE,
#         cluster_columns = FALSE,
#         column_title = "MERFISH",
#         row_title = "scRNA")
#
#
#

col_fun = circlize::colorRamp2(breaks = c(-1, -0.5, 0, 0.5, 1),
                               colors = c("blue", 'dodgerblue', "white", "tomato","red3"))
svglite::svglite("figures/annotated/merfish2scrna.global.marker_correlation.heatmap.svg", height=5, width=6)
Heatmap(as.matrix(smat),
        name = 'Spearman',
        col = col_fun,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        column_title = "MERFISH",
        row_title = "scRNA")
dev.off()



#### cell type detailed


### ---- ###
library(ComplexHeatmap)
# Jaccard similarity index function
jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}


## scRNA vs MERFISH marker qc
merfish <- qs::qread("data/merfish/BAYSOR/seurat_objects/ns-atlas.merfish_baysor.scanvi_intergated.cellcharter.seurat_object.QS") %>%
  subset(subset = cell_type.detailed != "Doublet") %>%
  UpdateSeuratObject()

scrna <- qs::qread("data/scrna/seurat_objects/normal_skin.scRNA.harmony.integrated.seurat_object.reclustered.annotated.filtered.QS") %>%
  UpdateSeuratObject()

genes <- intersect(common_genes, rownames(scrna))
merfish.counts <- as.matrix(merfish.counts[genes,])
merfish.celltypes <- c("Bas KC", "Spn KC", "Grn KC",
                       'Ecc Duct', "Ecc Gland","Seb", "HFE",
                       "Fibro","Peri", "EC", "LEC", "Schwann",
                       "Plasma", "Lym", "LC", "Mac/DC", "Mast")
scrna.celltypes <- c("Cyc KC", "Bas KC", "Spn KC", "Grn KC",
                     "Ecc Gland","Melano", "Bas HFE", "Diff HFE",
                     "Fib", "Peri", "EC", "LEC",
                     "Plasma", "CD4+ T Lym", "CD8+ T Lym",
                     "LC", "DC", "Mono", "Mac", "Mast")


scrna.counts <- AggregateExpression(scrna, features = genes, group.by = "cell_type.broad")[[1]]
scrna.counts <- as.matrix(scrna.counts[genes, scrna.celltypes])

merfish.counts <- AggregateExpression(merfish, features = genes, group.by = "cell_type.broad")[[1]]
merfish.counts <- as.matrix(merfish.counts[genes, merfish.celltypes])

spearman.df <- lapply(merfish.celltypes, function(m){
  lapply(scrna.celltypes, function(s){
    cor.test(merfish.counts[, m], scrna.counts[, s], method = 'spearman') %>%
      broom::tidy() %>%
      mutate(scrna.cell_type = s,
             merfish.cell_type = m)
  }) %>% bind_rows()
}) %>%
  bind_rows() %>%
  mutate(p.adj = p.adjust(p.value, method = "BH"))

spearman_cormat <- spearman.df %>%
  dplyr::select(scrna.cell_type, merfish.cell_type, estimate) %>%
  pivot_wider(names_from = "merfish.cell_type", values_from = "estimate") %>%
  data.frame()


rownames(spearman_cormat) <- spearman_cormat$scrna.cell_type
spearman_cormat$scrna.cell_type <- NULL
colnames(spearman_cormat) <- merfish.celltypes

library(ComplexHeatmap)
svglite::svglite("figures/annotated/merfish2scrna.broad.all_genes.spearman_correlation.redblue.heatmap.svg", height=10, width=10)
Heatmap(as.matrix(spearman_cormat),
        name = 'Spearman', # col =  col_fun,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        column_title = "MERFISH",
        row_title = "scRNA",
        width = ncol(spearman_cormat)*unit(3, "mm"),
        height = nrow(spearman_cormat)*unit(3, "mm")
)
dev.off()



### DETAILED
genes <- intersect(common_genes, rownames(scrna))
merfish.celltypes <- levels(merfish$cell_type.detailed)
scrna.celltypes <- intersect(merfish.celltypes, unique(scrna$cell_type.detailed))


scrna.counts <- AggregateExpression(scrna, features = genes, group.by = "cell_type.detailed")[[1]]
scrna.counts <- as.matrix(scrna.counts[genes, scrna.celltypes])

merfish.counts <- AggregateExpression(merfish, features = genes, group.by = "cell_type.detailed")[[1]]
merfish.counts <- as.matrix(merfish.counts[genes, merfish.celltypes])

spearman.df <- lapply(merfish.celltypes, function(m){
  lapply(scrna.celltypes, function(s){
    cor.test(merfish.counts[, m], scrna.counts[, s], method = 'spearman') %>%
      broom::tidy() %>%
      mutate(scrna.cell_type = s,
             merfish.cell_type = m)
  }) %>% bind_rows()
}) %>%
  bind_rows() %>%
  mutate(p.adj = p.adjust(p.value, method = "BH"))

spearman_cormat <- spearman.df %>%
  dplyr::select(scrna.cell_type, merfish.cell_type, estimate) %>%
  pivot_wider(names_from = "merfish.cell_type", values_from = "estimate") %>%
  data.frame()


rownames(spearman_cormat) <- spearman_cormat$scrna.cell_type
spearman_cormat$scrna.cell_type <- NULL
colnames(spearman_cormat) <- merfish.celltypes

library(ComplexHeatmap)
svglite::svglite("figures/annotated/merfish2scrna.detailed.all_genes.spearman_correlation.redblue.heatmap.svg", height=10, width=10)
Heatmap(as.matrix(spearman_cormat),
        name = 'Spearman', # col =  col_fun,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        column_title = "MERFISH",
        row_title = "scRNA",
        width = ncol(spearman_cormat)*unit(3, "mm"),
        height = nrow(spearman_cormat)*unit(3, "mm")
)
dev.off()

## prepare public scRNAseq data

# 1. Libraries
library(tidyverse)
library(Matrix)
library(Seurat)
library(SeuratObject)
library(rhdf5)
library(reticulate)
library(bit64)
library(loomR)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(biomaRt)
library(DoubletFinder)
setwd("D:/Dropbox/_projects/NormalSkinAtlas")

sc <- import("scanpy")
ad <- import("anndata")

## Ganier - PNAS - 2024
adata <- sc$read_h5ad("D:/Dropbox/_data/scRNA/GanierPNAS2024/ganier_pnas_2024.healthy.anndata.h5ad")
counts <- Matrix(t(adata$layers['counts']), sparse = TRUE)
gc()
colnames(counts) <- adata$obs_names$to_list()
rownames(counts) <- adata$var_names$to_list()

meta <- adata$obs
rm(adata)

meta <- meta[, c("01_sample", '03_location', '05_subcelltypes', 'percent.mt')]
colnames(meta) <- c("sample_barcode",'anatomic_site', 'reported.cell_type', "pct.mito")
meta$orig.cell_barcode <- rownames(meta)
meta$study_id <- "Ganier_PNAS_2024"
meta$donor_id <- meta$sample_barcode
meta$sample <- meta$sample_barcode
meta$sample_barcode <- paste0("CG_", meta$sample)
meta$donor_sex <- "male"
meta$donor_sex[meta$sample == 'face_cheek3b'] <- 'female'


samples <- levels(meta$sample)
age <- c(85, 62, 81, 90, 77, 86, 85, 70, 80, 90, 73, 56, 59, 78)
names(age) <- samples
meta$donor_age <- as.numeric(as.character(factor(as.character(meta$sample), labels = age)))
rm(age); gc()

meta <- meta[! (meta$reported.cell_type %in% c("Chondrocyte", "Skeletal muscle cells", "Chondrocytes")),]

counts <- counts[, rownames(meta)]
levels(meta$reported.cell_type) <- c("T Helper", "NK", "Fib", "T Cell", "T Reg", "Mac", "DC", "Fib", "Peri", "Fib", "Peri", "EC",
                                     "T Cell", "ILC/NK", "Spn/Grn KC", "B Cell", "Basal KC", "Mono", "Mast", "Melano", "Chondrocyte", "DC", "LEC", "Plasma", "Fib",
                                     "Mig DC", "DC", "Schwann", "SM", "Skeletal Muscle")
meta$reported.cell_type <- as.character(meta$reported.cell_type)

ganier <- CreateSeuratObject(counts, meta.data = meta, project = 'Ganier_PNAS_2024', min.cells = 3, min.features = 200)
ganier <- UpdateSeuratObject(ganier)

summary(ganier@meta.data$pct.mito)
summary(ganier@meta.data$nCount_RNA)
ganier <- subset(ganier, subset =  pct.mito < 10 & nCount_RNA > 100)

saveRDS(ganier, "data/published/scRNA/Ganier_PNAS_2024.normal_skin.seurat_object.RDS")
rm(ganier, counts, meta); gc()
gc()


## Weidemann - Cell Reports - 2023
# GEO files are truncated, cellxgene seems to be more comprehensive
## use data from cellxgene
weidemann <- readRDS("D:/Dropbox/_data/scRNA/Weidemann2023CellReports/weidemann_cellreports_2023.cellxgene.seurat_object.rds")
weidemann@meta.data$anatomic_site <- weidemann@meta.data$location
weidemann@meta.data$donor_sex <- weidemann@meta.data$sex
weidemann@meta.data$donor_age <- as.numeric(str_split(weidemann@meta.data$development_stage, "-", simplify = TRUE)[,1])
weidemann@meta.data$donor_id <- gsub("subject_", "JW_Subj0", weidemann@meta.data$donor_id)
levels(weidemann@meta.data$tissue) <- c("Epi", 'Derm', "WS")
weidemann@meta.data$sample_barcode <- paste0(weidemann@meta.data$donor_id, "_", stringr::str_to_title(weidemann@meta.data$anatomic_site),'-', weidemann@meta.data$tissue)
weidemann@meta.data$orig.cell_barcode <- rownames(weidemann@meta.data)
weidemann@meta.data$reported.cell_type <- weidemann@meta.data$author_cell_type
weidemann@meta.data$study_id <- "Wiedemann_CellReports_2023"
weidemann@meta.data$pct.mito <- weidemann@meta.data$percent.mt
levels(weidemann@meta.data$reported.cell_type) <- c("Spn KC"," Spn KC", "EC", "Fib", "Fib", "Immune", "Melano", "Peri", "Grn KC", "Basal KC", "Basal KC")


counts <- GetAssayData(weidemann[['RNA']], layer='counts')
meta <- weidemann@meta.data[, c("orig.cell_barcode", "study_id", 'sample_barcode', "donor_id",'donor_age', 'donor_sex', "anatomic_site", 'reported.cell_type', "pct.mito")]
rownames(counts) <- as.character(weidemann@assays[["RNA"]]@meta.features$feature_name)

weidemann <- CreateSeuratObject(counts = counts, meta.data = meta, project = 'Weidemann_CellReports_2023', min.cells = 3, min.features = 200)

summary(weidemann@meta.data$pct.mito)
summary(weidemann@meta.data$nCount_RNA)
weidemann <- subset(weidemann, subset =  pct.mito < 10 & nCount_RNA > 100)
weidemann <- UpdateSeuratObject(weidemann)
saveRDS(weidemann, file = "data/published/scRNA/Wiedemann_CellReports_2023.seurat_object.rds")
rm(weidemann, counts, meta); gc()

## - load data from Raman
int <- readRDS("D:/Dropbox/_data/scRNA/normal_skin_multidata/merged_human_AD_int_by_dataset_v3.rds")
int <- subset(int, subset = mitoPercent < 10 & nCount_RNA > 100 & condition == "Healthy")
gc()
int@meta.data$pct.mito <- int@meta.data$mitoPercent
int@meta.data$study_id <- factor(int@meta.data$dataset,
                                 levels = c("EXCLUDE", "Brunner23", "Brunner20", "Guttman", "Tabib"),
                                 labels = c("EXCLUDE", "Alkon_JAllergyClinImmunol_2023",
                                            "Rohajn_JAllergyClinImmunol_2020", "EXCLUDE", "EXCLUDE"))

int <- subset(int, subset = study_id != "EXCLUDE")
gc()
brunner20 <- subset(int, study_id == "Rohajn_JAllergyClinImmunol_2020")
brunner23 <- subset(int, study_id == "Alkon_JAllergyClinImmunol_2023")
rm(int); gc()

## Rohajn
meta <- data.frame(donor_id = c("HC3", "HC4", 'HC6', 'HC7'),
                   donor_sex = rep("female", times=4),
                   donor_age = c(47, 49, 27, 42),
                   anatomic_site = rep('antecubital fossa', times=4)
                   )

meta$dataset_sample <- paste0("Brunner20_", meta$donor_id)
meta$study_id <- "Rohajn_JAllergyClinImmunol_2020"

brunner20@meta.data$orig.cell_barcode <- rownames(brunner20@meta.data)
metadf <- brunner20@meta.data %>% inner_join(meta)
rownames(metadf) <- metadf$orig.cell_barcode

brunner20@meta.data <- metadf
brunner20@meta.data$reported.cell_type <- brunner20@meta.data$celltype
brunner20@meta.data$donor_id <- paste0("TR_", brunner20@meta.data$donor_id)
brunner20@meta.data$sample_barcode <- paste0(brunner20@meta.data$donor_id, "_AF-WS")
counts <- GetAssayData(brunner20[['RNA']], layer='counts')
meta <- brunner20@meta.data[, c("orig.cell_barcode", "study_id", 'sample_barcode', "donor_id",'donor_age', 'donor_sex', "anatomic_site", 'reported.cell_type', "pct.mito")]

brunner20 <- CreateSeuratObject(counts = counts, meta.data = meta, project = 'Rohajn_JAllergyClinImmunol_2020', min.cells = 3, min.features = 200)
brunner20 <- subset(brunner20, subset = nCount_RNA > 100 & pct.mito < 10)
saveRDS(brunner20, file = 'data/published/scRNA/Rohajn_JAllergyClinImmunol_2020.normal_skin.seurat_object.RDS')
rm(counts, meta, metadf, brunner20); gc()

## Alkon
brunner23@meta.data$anatomic_site <- "trunk"
brunner23@meta.data$orig.cell_barcode <- rownames(brunner23@meta.data)
brunner23@meta.data$reported.cell_type <- brunner23@meta.data$celltype
meta <- data.frame(dataset_sample = c("Brunner23_HC1", "Brunner23_HC2", "Brunner23_HC3", "Brunner23_HC4"),
                   donor_age = c(51, 48, 56, 44),
                   donor_sex = c("female", "male", "female", 'female')
)

metadf <- brunner23@meta.data %>% inner_join(meta)
rownames(metadf) <- metadf$orig.cell_barcode
brunner23@meta.data <- metadf
rm(meta, metadf)
brunner23@meta.data$donor_id <- gsub("Brunner23_", "NA_", brunner23@meta.data$dataset_sample)
brunner23@meta.data$sample_barcode <- paste0(brunner23@meta.data$donor_id, "_Trunk-WS")
counts <- GetAssayData(brunner23[['RNA']], layer='counts')
meta <- brunner23@meta.data[, c("orig.cell_barcode", "study_id", 'sample_barcode', "donor_id",'donor_age', 'donor_sex', "anatomic_site", 'reported.cell_type', "pct.mito")]
brunner23 <- CreateSeuratObject(counts = counts, meta.data = meta, project = 'Alkon_JAllergyClinImmunol_2023', min.cells = 3, min.features = 200)
brunner23 <- subset(brunner23, subset = nCount_RNA > 100 & pct.mito < 10)
saveRDS(brunner23, file = 'data/published/scRNA/Alkon_JAllergyClinImmunol_2023.normal_skin.seurat_object.RDS')
rm(meta, counts, brunner23); gc()


## Sole-Boldo
soleboldo <- readRDS("D:/Dropbox/_data/scRNA/SoleBoldoMolSystBiol2022/sole-boldo_communbiol_2020.cellxgene.seurat_object.rds")
soleboldo@meta.data$study_id <- "SoleBoldo_CommunBiol_2020"
soleboldo@meta.data$anatomic_site <- "groin"
soleboldo@meta.data$donor_sex <- tolower(soleboldo@meta.data$sex)
soleboldo@meta.data$donor_age <- as.numeric(str_split(soleboldo@meta.data$development_stage, "-", simplify = TRUE)[,1])
soleboldo@meta.data$donor_id <- paste0("LS_", soleboldo@meta.data$donor_id)
soleboldo@meta.data$sample_barcode <- paste0(soleboldo@meta.data$donor_id, "_Groin-WS")
soleboldo@meta.data$orig.cell_barcode <- rownames(soleboldo@meta.data)
soleboldo@meta.data$reported.cell_type <- soleboldo@meta.data$Celltype
levels(soleboldo@meta.data$reported.cell_type) <- c("Spn/Grn KC", "Basal KC", "RBC", "LEC", "Mac/DC", "Melano", "Fib", "Peri", "Fib", "Fib","Fib", "T Cell", "EC")
soleboldo@meta.data$reported.cell_type <- as.character(soleboldo@meta.data$reported.cell_type)

meta <- soleboldo@meta.data[, c("orig.cell_barcode", "study_id", 'sample_barcode', "donor_id",'donor_age', 'donor_sex', "anatomic_site", 'reported.cell_type')]
meta <- meta[meta$reported.cell_type != "RBC", ]

counts <- GetAssayData(soleboldo[['RNA']], layer='counts')
rownames(counts) <- as.character(soleboldo@assays[["RNA"]]@meta.features$feature_name)
counts <- counts[, rownames(meta)]

soleboldo <- CreateSeuratObject(counts = counts, meta.data = meta, project = 'SoleBoldo_CommunBiol_2020', min.cells = 3, min.features = 200)
soleboldo[['pct.mito']] <- PercentageFeatureSet(soleboldo, pattern = '^MT-')
soleboldo <- subset(soleboldo, subset = nCount_RNA > 100 & pct.mito < 10)
saveRDS(soleboldo, "data/published/scRNA/Sole-Boldo_CommunBiol_2020.seurat_object.rds")
rm(soleboldo, counts, meta); gc()

## NS-Int from Andrew
nsint <- readRDS("D:/Dropbox/_data/scRNA/normal_skin_multidata/ns_meta.Rds")

nsint@meta.data$pct.mito <- nsint@meta.data$percent.mt
nsint@meta.data$reported.cell_type <- factor(nsint@meta.data$subcell)
nsint@meta.data$donor_id <- paste0(nsint@meta.data$dataset, "_", nsint@meta.data$sample)
nsint@meta.data$sample_barcode <- nsint@meta.data$donor_id
nsint@meta.data$orig.cell_barcode <- rownames(nsint@meta.data)
nsint@meta.data <- nsint@meta.data[, c("dataset", "sample", "donor_id", "sample_barcode", "reported.cell_type", "orig.cell_barcode", "pct.mito")]
levels(nsint@meta.data$reported.cell_type) <- c("ASDC", "B Cell", "Basal KC", "Basal KC", "Basal KC", "CD1C+ DC", "CLEC9A+ DC", "Cyc KC", "Cyc KC", "Cyc KC", "Eccrine", "Eccrine", "EC",
                                                "Fib", "Grn KC", "HF", "IB", "INF", "INF", "INF", "INF", "IRS", "KC", "LC", "Mac", "MDSC", "Melano", "Spn KC", "Spn KC", "Multiplet", "Multiplet", "Basal KC", "Cyc KC", "NK", "OB",
                                                "pDC", "Seb", "SB", "Seb", "Seb", "Spn KC", "Spn KC", "Spn KC", "Spn KC", "Spn KC", "T Cell")
nsint@meta.data$reported.cell_type <- as.character(nsint@meta.data$reported.cell_type)
nsint@meta.data$reported.cell_type[is.na(nsint@meta.data$reported.cell_type)] <- "Unknown"

nsint <- subset(nsint, subset = pct.mito < 10)
nsint <- subset(nsint, subset = reported.cell_type != "Multiplet")

tabib <- subset(nsint, dataset == 'TT')
zou <- subset(nsint, dataset == 'ZZ')
ji <- subset(nsint, dataset == 'AJ')
takahashi <- subset(nsint, dataset == 'RT')
cheng <- subset(nsint, dataset == 'JC')
rm(nsint); gc()

### Tabib
counts <- GetAssayData(tabib[['RNA']], layer="counts")
meta <- tabib@meta.data
meta$study_id <- "Tabib_JID_2018"
meta$anatomic_site <- 'forearm'
samples <- unique(meta$donor_id)
sex <- c('male', 'male', 'female', 'female','female','male')
names(sex) <- samples
meta$donor_sex <- factor(as.character(meta$sample), labels = sex)

rm(sex); gc()
age <- c(63, 54, 66, 23, 62, 24)
names(age) <- samples
meta$donor_age <- as.numeric(as.character(factor(as.character(meta$sample), labels = age)))
rm(age); gc()

meta <- meta[, c("orig.cell_barcode", "study_id", 'sample_barcode', "donor_id",'donor_sex', "anatomic_site", 'reported.cell_type', 'pct.mito')]

tabib <- CreateSeuratObject(counts = counts, meta.data = meta, project = 'Tabib_JID_2018', min.cells=3, min.features = 200)
tabib <- subset(tabib, subset = nCount_RNA > 100 & pct.mito < 10)
saveRDS(tabib, file = "data/published/scRNA/Tabib_JID_2018.seurat_object.RDS")
rm(meta, counts, samples, tabib); gc()

## Zou Dev Cell 2020
counts <- GetAssayData(zou[['RNA']], layer="counts")
meta <- zou@meta.data
meta$study_id <- "Zou_DevCell_2021"
meta$donor_age <- as.numeric(as.character(substr(meta$sample, 2, nchar(meta$sample))))
meta$anatomic_site <- 'eyelid'
meta$donor_sex <- "female"
meta <- meta[, c("orig.cell_barcode", "study_id", 'sample_barcode', "donor_id",'donor_sex', "anatomic_site", 'reported.cell_type', 'pct.mito')]

zou <- CreateSeuratObject(counts = counts, meta.data = meta, project = 'Zou_DevCell_2021', min.cells=3, min.features = 200)
zou <- subset(zou, subset = nCount_RNA > 100 & pct.mito < 10)
saveRDS(zou, file = "data/published/scRNA/Zou_DevCell_2021.seurat_object.RDS")
rm(counts, meta, zou); gc()


## Ji Cell 2020
counts <- GetAssayData(ji[['RNA']], layer="counts")
meta <- ji@meta.data
meta$study_id <- "Ji_Cell_2020"

demo <- read.csv('D:/Dropbox/_data/scRNA/normal_skin_multidata/ji_cell_2020.sample_metadata.csv')
meta <- inner_join(meta, demo)
meta <- data.frame(meta)
rownames(meta) <- meta$orig.cell_barcode
meta <- meta[, c("orig.cell_barcode", "study_id", 'sample_barcode', "donor_id",'donor_sex', "anatomic_site", 'reported.cell_type', 'pct.mito')]

ji <- CreateSeuratObject(counts = counts, meta.data = meta, project = 'Ji_Cell_2020', min.cells=3, min.features = 200)
ji <- subset(ji, nCount_RNA > 100 & pct.mito < 10)
saveRDS(ji, file = "data/published/scRNA/Ji_Cell_2020.normal_skin.seurat_object.RDS")
rm(counts, meta, demo, ji); gc()

## Takahashi
counts <- GetAssayData(takahashi[['RNA']], layer="counts")
meta <- takahashi@meta.data
meta$study_id <- "Takahashi_JID_2020"
meta$anatomic_site <- 'scalp'

meta <- meta[, c("orig.cell_barcode", "study_id", "sample_barcode", "donor_id", 'anatomic_site', "reported.cell_type", "pct.mito")]

takahashi <- CreateSeuratObject(counts, meta.data = meta, project = "Takahashi_JID_2020", min.cells=3, min.features = 200)
takahashi <- subset(takahashi, nCount_RNA > 100 & pct.mito < 10)
saveRDS(takahashi, "data/published/scRNA/Takahashi_JID_2020.seurat_object.RDS")


## Cheng
meta <- readRDS("D:/Dropbox/_data/scRNA/ChengCellReports2018/Cheng_CellReports_2018.cellxgene.seurat_object.rds")@meta.data
donor_sex <- as.character(deframe(unique(meta[, c('donor_id', "sex")]))[unique(cheng@meta.data$sample)])
anatomic_sites <- c('abd4'="abdomen", "br41epi"='breast', 'br53epi'='breast', 's11'='scalp', 'scalp26'='scalp', 'scalp32'='scalp')

cheng@meta.data$anatomic_site <- as.character(factor(as.character(cheng@meta.data$sample), labels = anatomic_sites))
cheng@meta.data$donor_sex <- as.character(factor(as.character(cheng@meta.data$sample), labels = donor_sex))
rm(meta, donor_sex, anatomic_sites)

cheng@meta.data$study_id <- "Cheng_CellReports_2018"
meta <- cheng@meta.data[, c("orig.cell_barcode", "study_id", 'sample_barcode', "donor_id", 'donor_sex', "anatomic_site", 'reported.cell_type', 'pct.mito')]
counts <- GetAssayData(cheng[['RNA']], layer='counts')

cheng <- CreateSeuratObject(counts = counts, meta.data = meta, project = 'Cheng_CellReports_2018', min.cells = 3, min.features = 200)
cheng <- subset(cheng, subset = nCount_RNA > 100 & pct.mito < 10)
saveRDS(cheng, "data/published/scRNA/Cheng_CellReports_2018.seurat_object.RDS")
rm(counts, meta, cheng); gc()


## GTEX
adata <- sc$read_h5ad("D:/Dropbox/_data/scRNA/GTEx_scRNA/gtex.lowerleg.normal_skin.raw_counts.h5ad")
counts <- Matrix(t(adata$layers['counts']), sparse = TRUE)
gc()
colnames(counts) <- adata$obs_names$to_list()
rownames(counts) <- adata$var_names$to_list()

meta <- adata$obs


meta$sample_barcode <- meta$`Sample ID short`
meta$donor_sex <- tolower(meta$donor_sex)
levels(meta$reported.cell_type) <- c("Adipo", "LEC", "EC", "EC", "Basal KC", "Grn KC", "Spn KC", "Spn KC", "Fib", "DC/Mac", "LC", "T Cell", "Mast", "Melano", "Peri", "Peri", "Seb", "Seb", "Seb", "Eccrine", "Unknown", "Unknown")

meta$reported.cell_type <- as.character(meta$reported.cell_type)
meta <- meta[, c("orig.cell_barcode", "study_id", 'sample_barcode', "donor_id", 'donor_sex', "anatomic_site", 'reported.cell_type')]

gtex <- CreateSeuratObject(counts, meta.data = meta, project = "GTEx_Eraslan_Science_2022", min.cells = 3, min.features = 200)
gtex[['pct.mito']] <- PercentageFeatureSet(gtex, pattern = '^MT-')

gtex <- subset(gtex, subset = nCount_RNA > 100 & pct.mito < 10)
saveRDS(gtex, "data/published/scRNA/GTEx_Eraslan_Science_2022.normal_skin.seurat_object.RDS")
rm(adata, meta, counts, gtex); gc()


## Belote Nat Comms 2021
counts <- read.csv("D:/Dropbox/_data/scRNA/BeloteNatComms2021/GSE151091_raw_matrix.csv", row.names = 1)
ptdf <- read.csv(gzfile("D:/Dropbox/_data/scRNA/BeloteNatComms2021/GSE151091_Metadata.csv.gz"))
ptdf <- ptdf %>% dplyr::select(patient, age, sex, anatomical_location, plate) %>% unique() %>% data.frame()
colnames(ptdf) <- c("donor_id", "donor_age", "donor_sex", "anatomic_site", "plate_id")
cell_names <- colnames(counts)
mt <- data.frame(str_split(cell_names, pattern = "_", simplify = TRUE))
mt$cell_barcode <- cell_names
mt <- mt[nchar(mt$X4) == 0, ]
colnames(mt) <- c("donor_id", "plate_index", "plate_id", "X4", "cell_barcode")
mt$plate_index <- paste0(mt$plate_index, "_", mt$plate_id)
mt$X4 <- NULL
meta <- inner_join(mt, ptdf) %>% data.frame()
rownames(meta) <- meta$cell_barcode
meta$study_id <- "Belote_NatcellBio_2021"
meta$orig.cell_barcode <- meta$cell_barcode
meta$sample_barcode <- meta$donor_id
counts <- counts[, meta$cell_barcode]

belote <- CreateSeuratObject(counts, meta.data = meta, project = "Belote_NatCellBio_2021", min.cells = 3, min.features = 200)

belote[['pct.ribo']] <- PercentageFeatureSet(belote, pattern = "^RPL|^RPS")
belote[['pct.mito']] <- PercentageFeatureSet(belote, pattern = "^MT-")
belote[['pct.hemo']] <- PercentageFeatureSet(belote, pattern = "^HB")
belote <- subset(belote, subset = nCount_RNA > 100 & pct.mito < 10)
saveRDS(belote, "data/published/scRNA/Belote_NatCellBio_2021.seurat_object.RDS")


## VorstandlechnerNatComms2021
seurat_list <- lapply(1:3, function(x){
  obj <- Read10X(paste0("D:/Dropbox/_data/scRNA/VorstandlechnerNatComms2021/outs/human_skin_", x), unique.features = TRUE)
  obj <- CreateSeuratObject(obj, project = paste0("VorstandlechnerNatComms2021_NS-0", x), min.cells = 3, min.features = 200)
  obj@meta.data$sample_barcode <- paste0("VV_NS0", x, "-Abd")
  obj@meta.data$anatomic_site <- "abdomen"
  obj@meta.data$donor_id <- paste0("VV_NS0", x, "-Abd")
  obj@meta.data$study_id <- "VorstandlechnerNatComms2021"
  obj@meta.data$donor_sex <- "Female"
  obj[['pct.ribo']] <- PercentageFeatureSet(obj, pattern = "^RPL|^RPS")
  obj[['pct.mito']] <- PercentageFeatureSet(obj, pattern = "^MT-")
  obj[['pct.hemo']] <- PercentageFeatureSet(obj, pattern = "^HB")
  obj@meta.data$cell_barcode <-  rownames(obj@meta.data)
  obj <- subset(obj, subset = nCount_RNA > 100 & pct.mito < 10)

  return(obj)
})

vorstandlechner <- merge(seurat_list[[1]], seurat_list[2:3], add.cell.ids = paste0("VorstandlechnerNatComms2021_NS-0", 1:3))
vorstandlechner[['RNA']] <- JoinLayers(vorstandlechner[['RNA']])
saveRDS(vorstandlechner, "data/published/scRNA/Vorstandlechner_NatComms_2021.seurat_object.RDS")
rm(seurat_list); gc()


seurat_list <- lapply(c("axillar", "groin"), function(x){
  obj <- Read10X(paste0("D:/Dropbox/_data/scRNA/YuImmunity2024/GSE158955_RAW/NS_", x), unique.features = TRUE)
  obj <- CreateSeuratObject(obj, project = paste0("YuImmunity2024_NS-", x), min.cells = 3, min.features = 200)
  obj@meta.data$sample_barcode <- paste0("YWW_NS-", x)
  obj@meta.data$anatomic_site <- x
  obj@meta.data$donor_id <- paste0("YWW_NS-", x)
  obj@meta.data$study_id <- "YuImmunity2024"
  obj@meta.data$donor_sex <- "Female"
  obj[['pct.ribo']] <- PercentageFeatureSet(obj, pattern = "^RPL|^RPS")
  obj[['pct.mito']] <- PercentageFeatureSet(obj, pattern = "^MT-")
  obj[['pct.hemo']] <- PercentageFeatureSet(obj, pattern = "^HB")
  obj@meta.data$cell_barcode <-  rownames(obj@meta.data)
  obj <- subset(obj, subset = nCount_RNA > 100 & pct.mito < 10)

  return(obj)
})

yu <- merge(seurat_list[[1]], seurat_list[[2]], add.cell.ids = paste0("YuImmunity2024-NS", c("ax", "gr")))
yu[['RNA']] <- JoinLayers(yu[['RNA']])
yu <- JoinLayers(yu)

saveRDS(yu, "data/published/scRNA/Yu_Immunity_2024.seurat_object.RDS")
rm(seurat_list); gc()



### put it all together
seurat_list <- list(
  Ji_Cell_2020 = readRDS("data/published/scRNA/Ji_Cell_2020.normal_skin.seurat_object.RDS"),
  Wiedemann_CellReports_2023 = readRDS("data/published/scRNA/Wiedemann_CellReports_2023.seurat_object.rds"),
  Zou_DevCell_2021 = readRDS("data/published/scRNA/Zou_DevCell_2021.seurat_object.RDS"),
  Tabib_JID_2018 = readRDS("data/published/scRNA/Tabib_JID_2018.seurat_object.RDS"),
  GTEx_Eraslan_Science_2022 = readRDS("data/published/scRNA/GTEx_Eraslan_Science_2022.normal_skin.seurat_object.RDS"),
  Ganier_PNAS_2024 = readRDS("data/published/scRNA/Ganier_PNAS_2024.normal_skin.seurat_object.RDS"),
  SoleBoldo_CommunBiol_2020 = readRDS("data/published/scRNA/Sole-Boldo_CommunBiol_2020.seurat_object.rds"),
  Rohajn_JAllergyClinImmunol_2020 = readRDS("data/published/scRNA/Rohajn_JAllergyClinImmunol_2020.normal_skin.seurat_object.RDS"),
  Alkon_JAllergyClinImmunol_2023 = readRDS("data/published/scRNA/Alkon_JAllergyClinImmunol_2023.normal_skin.seurat_object.RDS"),
  Cheng_CellReports_2018 = readRDS("data/published/scRNA/Cheng_CellReports_2018.seurat_object.RDS"),
  Takahashi_JID_2020 = readRDS("data/published/scRNA/Takahashi_JID_2020.seurat_object.RDS"),
  Belote_NatCellBio_2021 = readRDS("data/published/scRNA/Belote_NatCellBio_2021.seurat_object.RDS"),
  Vorstandlechner_NatComms_2021 = readRDS("data/published/scRNA/Vorstandlechner_NatComms_2021.seurat_object.RDS"),
  Yu_Immunity_2024 = readRDS("data/published/scRNA/Yu_Immunity_2024.seurat_object.RDS")
)
saveRDS(seurat_list, file = "data/published/scRNA/public_scrna.normal_skin.seurat_list.RDS")

nscomb <- merge(seurat_list[[1]], y = seurat_list[2:length(seurat_list)], merge.data = FALSE)
rm(seurat_list); gc()
nscomb <- UpdateSeuratObject(nscomb)
nscomb <- JoinLayers(nscomb)
gc()

nscomb[['pct.ribo']] <- PercentageFeatureSet(nscomb, pattern = "^RPL|^RPS")
nscomb[['pct.hemo']] <- PercentageFeatureSet(nscomb, pattern = "^HB")
nscomb <- subset(nscomb, subset = pct.hemo < 5 & pct.ribo < 60)

nscomb@meta.data <- nscomb@meta.data[, c("orig.ident", 'nCount_RNA', 'nFeature_RNA', "pct.mito", "pct.ribo", "pct.hemo", "study_id", 'sample_barcode', "donor_id", "donor_sex", "donor_age", "anatomic_site", "reported.cell_type")]
nscomb@meta.data$anatomic_site.detailed <- nscomb@meta.data$anatomic_site
nscomb@meta.data$anatomic_site[nscomb@meta.data$anatomic_site == 'groin'] <- "inguinal fold"
nscomb@meta.data$anatomic_site[nscomb@meta.data$anatomic_site %in% c("nose", "cheek", "eyelid", 'temple', 'forehead', "ear")] <- 'face'
nscomb@meta.data$anatomic_site[nscomb@meta.data$anatomic_site %in% c("trunk", "chest")] <- "trunk"
nscomb@meta.data$anatomic_site[nscomb@meta.data$anatomic_site %in% c('shin', "calf", "medial calf")] <- "lower leg"
nscomb@meta.data$anatomic_site[nscomb@meta.data$anatomic_site %in% c('right thigh', "left thigh", "thigh")] <- "upper leg"
nscomb@meta.data$anatomic_site[nscomb@meta.data$anatomic_site %in% c('forearm', "anterior forearm", "arm")] <- "forearm"
nscomb@meta.data$anatomic_site[nscomb@meta.data$anatomic_site %in% c( "heel", "arch")] <- "sole"
nscomb@meta.data$anatomic_site[nscomb@meta.data$anatomic_site %in% c( "ankle or upper foot")] <- "upper foot"


saveRDS(nscomb, "data/scrna/seurat_objects/normal_skin.scRNA.merged.seurat_object.RDS")

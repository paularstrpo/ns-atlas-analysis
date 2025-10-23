library(tidyverse)
library(Seurat)
library(Matrix)
library(rhdf5)
library(reticulate)
setwd("D:/Dropbox/_projects/NormalSkinAtlas")

## Ma Nat Comms 2024
folder.list <- list.dirs("data/published/Visium/ma-natcomms-2024/GSE225475_RAW/", full.names = TRUE, recursive = FALSE)

study_id <- "MaNatComms2024"
sample_ids <- str_split(folder.list, '/', simplify=TRUE)[,6]
sample_ids <- str_split(sample_ids, "_", simplify = TRUE)[,2]
names(folder.list) <- paste0(study_id, "_", sample_ids)
seurat_list <- lapply(folder.list, function(x){
  sample <- str_split(x, '/', simplify=TRUE)[,6]
  sample <- str_split(sample, "_", simplify = TRUE)[,2]

  obj <- Load10X_Spatial(x, slice = paste0(study_id,"_", sample))
  obj <- subset(obj, nCount_Spatial > 0)
  obj$pct.mito <- PercentageFeatureSet(obj, pattern = "^MT-")
  Project(obj) <- study_id
  obj$study_id <- study_id
  obj$sample_id <- paste0(study_id, "_", sample)
  obj$donor_id <- paste0(study_id, sample)
  obj$tissue_type <- substr(sample, start = 1, stop=2)
  obj$disease <- substr(sample, start = 1, stop=2)
  obj$condition <- ifelse(grepl("NS", sample), "NL", "LS")

  return(obj)
})
names(seurat_list) <- sample_ids
MaNatComms2024 <- merge(seurat_list[[1]], seurat_list[2:length(seurat_list)], add.cell.ids = names(seurat_list))
MaNatComms2024 <- JoinLayers(MaNatComms2024)
rm(seurat_list); gc()

qs::qsave(MaNatComms2024, file = 'data/visium/seurat_objects/MaNatComms2024.visium_data.raw.seurat_object.QS')
# rm(MaNatComms2024); gc()



## Mitamura Allergy 2023
folder.list <- list.dirs("data/published/Visium/mitamura-allergy-2023/GSE197023_RAW/", full.names = TRUE, recursive = FALSE)
study_id <- "MitamuraAllergy2023"
sample_ids <- paste0(study_id, "_", str_split(folder.list, '/', simplify=TRUE)[,6])
names(folder.list) <- sample_ids
seurat_list <- lapply(folder.list, function(x){
  sample <- str_split(x, '/', simplify=TRUE)[,6]
  info <- str_split(sample, "_", simplify = TRUE)
  obj <- Load10X_Spatial(x, slice = paste0(study_id,'_', sample))
  Project(obj) <- study_id
  obj <- subset(obj, nCount_Spatial > 0)
  obj$pct.mito <- PercentageFeatureSet(obj, pattern = "^MT-")
  obj$study_id <- study_id
  obj$sample_id <- paste0(study_id, "_", sample)
  obj$donor_id <- paste0(study_id, "_", info[,1], info[,2])
  obj$tissue_type <- ifelse(grepl("LS", sample), info[,1], "NS")
  obj$condition <- ifelse(grepl("LS", sample), "LS", "NL")
  obj$disease <- ifelse(info[,1] == "HE", "NS", info[,1])
  return(obj)
})
names(seurat_list) <- sample_ids
MitamuraAllergy2023 <- merge(seurat_list[[1]], seurat_list[2:length(seurat_list)], add.cell.ids = names(seurat_list))
MitamuraAllergy2023 <- JoinLayers(MitamuraAllergy2023)
rm(seurat_list); gc()
MitamuraAllergy2023$disease <- gsub("HE", "NS", MitamuraAllergy2023$disease)
qs::qsave(MitamuraAllergy2023, file = 'data/visium/seurat_objects/MitamuraAllergy2023.visium_data.raw.seurat_object.QS')
# rm(MitamuraAllergy2023); gc()


## Ganier PNAS 2024
folder.list <- list.dirs("data/published/Visium/ganier-pnas-2024/", full.names = TRUE, recursive = FALSE)
study_id <- "GanierPNAS2024"
sample_ids <- paste0(study_id, "_", str_split(folder.list, '/', simplify=TRUE)[,5])
seurat_list <- lapply(folder.list, function(x){
  sample <- str_split(x, '/', simplify=TRUE)[,5]
  obj <- Load10X_Spatial(x, slice = paste0(study_id, "_", sample))
  obj <- subset(obj, nCount_Spatial > 0)
  obj$pct.mito <- PercentageFeatureSet(obj, pattern = "^MT-")
  Project(obj) <- study_id
  obj$study_id <- study_id
  obj$sample_barcode <- paste0(study_id, "_", sample)
  return(obj)
})
names(seurat_list) <- sample_ids
GanierPNAS2024 <- merge(seurat_list[[1]], seurat_list[2:length(seurat_list)], add.cell.ids = names(seurat_list))
GanierPNAS2024 <- JoinLayers(GanierPNAS2024)
rm(seurat_list, sample_ids); gc()

sample_info <- read.delim("data/published/Visium/ganier-pnas-2024/E-MTAB-13084.sdrf.txt")
sample_info <- unique(sample_info[, c("Source.Name", "Characteristics.age.", "Characteristics.sex.",
                                      "Characteristics.sampling.site.", 'Characteristics.disease.', "Characteristics.sample.id.")])
sample_info$sample_id <- make.names(gsub("\\.", "", make.names(paste0(study_id, "_", trimws(gsub(" ", "", sample_info$Characteristics.sample.id.))))))
sample_info$sample_barcode <- paste0(study_id, "_", sample_info$Source.Name)
sample_info$disease <- ifelse(grepl("normal", sample_info$Characteristics.disease.), "NS", "BCC")
sample_info$tissue_type <- sample_info$disease
sample_info$condition <- ifelse(grepl("normal", sample_info$Characteristics.disease.), "NL", "LS")
sample_info$anatomic_site.detailed <- sample_info$Characteristics.sampling.site.
sample_info$anatomic_site.detailed[sample_info$anatomic_site.detailed == "inguinal part of abdomen"] <- "inguinal fold"
sample_info$anatomic_site <- sample_info$anatomic_site.detailed
sample_info$anatomic_site[grepl("face", sample_info$sample_description)] <- "face"
sample_info$donor_age <- sample_info$Characteristics.age.
sample_info$donor_sex <- sample_info$Characteristics.sex.
sample_info$donor_id <- sample_info$sample_id

## this one is a typo in the array express metadata file, fix it
sample_info$sample_barcode[sample_info$sample_barcode == "GanierPNAS2024_WSSKNKCLsp104466211"] <- "GanierPNAS2024_WSSKNKCLsp10446621"

sample_info <- unique(sample_info[, c("sample_barcode",'sample_id', "donor_id", "donor_age", "donor_sex", "tissue_type", "disease", "condition", "anatomic_site", "anatomic_site.detailed")])
names(GanierPNAS2024@images) <- sample_info$sample_id

GanierPNAS2024@meta.data$cell_barcode <- rownames(GanierPNAS2024@meta.data)
meta <- inner_join(sample_info, GanierPNAS2024@meta.data[, c("sample_barcode", "cell_barcode")]) %>% data.frame()
rownames(meta) <- meta$cell_barcode
meta <- meta[, c("donor_id", "donor_age", "donor_sex", "sample_id",'tissue_type', "disease", "condition", "anatomic_site", "anatomic_site.detailed")]
GanierPNAS2024 <- AddMetaData(GanierPNAS2024, metadata = meta)
rm(meta, sample_info)

GanierPNAS2024 <- UpdateSeuratObject(GanierPNAS2024)
qs::qsave(GanierPNAS2024, file = 'data/visium/seurat_objects/GanierPNAS2024.visium_data.raw.seurat_object.QS')
# rm(GanierPNAS2024, meta, sample_info); gc()

## Bergenstrahle Nat Biotech 2022
folder.list <- list.dirs("data/published/Visium/bergenstrahle-natbiotech-2022/", full.names = TRUE, recursive = FALSE)
sample_ids <- paste0("BergenstrahleNatBiotech2022_", str_split(folder.list, '/', simplify=TRUE)[,5])
seurat_list <- lapply(folder.list, function(x){
  sample <- str_split(x, '/', simplify=TRUE)[,5]
  img <- Read10X_Image(paste(x,'spatial', sep = '/'),
                       image.name = 'tissue_hires_image.png',
                       slice = paste0("BergenstrahleNatBiotech2022_", sample))
  orig.lowres <- img@scale.factors$lowres
  orig.hires <- img@scale.factors$hires
  img@scale.factors$lowres <- orig.hires
  img@scale.factors$hires <- orig.lowres
  obj <- Load10X_Spatial(x, image = img, slice = make.names(paste0("BergenstrahleNatBiotech2022_", sample)))

  Project(obj) <- "BergenstrahleNatBiotech2022"
  obj <- subset(obj, nCount_Spatial > 0)
  obj$pct.mito <- PercentageFeatureSet(obj, pattern = "^MT-")

  obj$study_id <- "BergenstrahleNatBiotech2022"
  obj$sample_id <- make.names(paste0("BergenstrahleNatBiotech2022_", sample))
  obj$donor_id <- "P03"
  obj$donor_sex <- "Female"
  obj$donor_age <- 71
  obj$anatomic_site <- "scalp"
  obj$condition <- "LS"
  obj$disease <- "SCC"
  obj$tissue_type <- "SCC"
  return(obj)
})
names(seurat_list) <- sample_ids
BergenstrahleNatBiotech2022 <- merge(seurat_list[[1]], seurat_list[2:length(seurat_list)], add.cell.ids = names(seurat_list))
rm(seurat_list, sample_ids, folder.list); gc()
BergenstrahleNatBiotech2022 <- JoinLayers(BergenstrahleNatBiotech2022)
qs::qsave(BergenstrahleNatBiotech2022, file = 'data/visium/seurat_objects/BergenstrahleNatBiotech2022.visium_data.raw.seurat_object.QS')
# rm(BergenstrahleNatBiotech2022); gc()

## Ji Cell 2020
folder.list <- list.dirs("data/published/Visium/ji-cell-2020/", full.names = TRUE, recursive = FALSE)
sample_ids <- paste0("JiCell2020_", str_split(folder.list, '/', simplify=TRUE)[,5])
seurat_list <- lapply(folder.list, function(x){
  sample <- str_split(x, '/', simplify=TRUE)[,5]
  obj <- Load10X_Spatial(x, slice = paste0("JiCell2020_", sample))

  Project(obj) <- "JiCell2020"
  obj <- subset(obj, nCount_Spatial > 0)
  obj$pct.mito <- PercentageFeatureSet(obj, pattern = "^MT-")

  obj$study_id <- "JiCell2020"
  obj$sample_id <- paste0("JiCell2020_", sample)
  obj$condition <- "LS"
  obj$disease <- "SCC"
  obj$tissue_type <- "SCC"

  return(obj)
})
names(seurat_list) <- sample_ids
JiCell2020 <- merge(seurat_list[[1]], seurat_list[2:length(seurat_list)], add.cell.ids = names(seurat_list))
rm(seurat_list, sample_ids, folder.list); gc()
JiCell2020@meta.data$donor_sex <- "Female"
JiCell2020@meta.data$donor_id <- ifelse(grepl("T28", JiCell2020@meta.data$sample_id), "P06", "P04")
JiCell2020@meta.data$anatomic_site <- ifelse(grepl("T28", JiCell2020@meta.data$sample_id), "chest", "hand")
JiCell2020@meta.data$donor_age <- ifelse(grepl("T28", JiCell2020@meta.data$sample_id), 96, 75)
JiCell2020 <- JoinLayers(JiCell2020)
qs::qsave(JiCell2020, file = 'data/visium/seurat_objects/JiCell2020.visium_data.raw.seurat_object.QS')
# rm(JiCell2020); gc()


## Thrane JID 2023
folder.list <- list.dirs("data/published/Visium/thrane-jid-2023/", full.names = TRUE, recursive = FALSE)
sample_ids <- paste0("ThraneJID2023_", str_split(folder.list, '/', simplify=TRUE)[,5])
seurat_list <- lapply(folder.list, function(x){
  sample <- str_split(x, '/', simplify=TRUE)[,5]
  obj <- Load10X_Spatial(x, slice = paste0("ThraneJID2023_", sample))

  Project(obj) <- "ThraneJID2023"
  obj <- subset(obj, nCount_Spatial > 0)
  obj$pct.mito <- PercentageFeatureSet(obj, pattern = "^MT-")

  obj$study_id <- "ThraneJID2023"
  obj$donor_id <- str_split(sample, "_", simplify = TRUE)[,1]
  obj$sample_id <- paste0("ThraneJID2023_", sample)
  obj$condition <- "NL"
  obj$tissue_type <- "NS"
  obj$disease <- "NS"
  return(obj)
})
names(seurat_list) <- sample_ids
ThraneJID2023 <- merge(seurat_list[[1]], seurat_list[2:length(seurat_list)], add.cell.ids = names(seurat_list))
rm(seurat_list, sample_ids, folder.list); gc()

ThraneJID2023@meta.data$anatomic_site <- ifelse(grepl("P03", ThraneJID2023@meta.data$sample_id), "scalp", "ear")
ThraneJID2023@meta.data$donor_age <- ifelse(grepl("P03", ThraneJID2023@meta.data$sample_id), 75, 71)
ThraneJID2023@meta.data$donor_sex <- ifelse(grepl("P03", ThraneJID2023@meta.data$sample_id), "Female", "Male")
ThraneJID2023 <- JoinLayers(ThraneJID2023)
qs::qsave(ThraneJID2023, file = 'data/visium/seurat_objects/ThraneJID2023.visium_data.raw.seurat_object.QS')


## Yu Immunity 2024
folder.list <- list.dirs("data/published/Visium/yu-immunity-2024/", full.names = TRUE, recursive = FALSE)
sample_ids <- paste0("YuImmunity2024_", str_split(folder.list, '/', simplify=TRUE)[,5])
# sample_ids <- sample_ids[c(1:8, 10)] # only visium, no cytassist
seurat_list <- lapply(folder.list, function(x){
  sample <- str_split(x, '/', simplify=TRUE)[,5]
  obj <- Load10X_Spatial(x, slice = paste0("YuImmunity2024_", sample))

  Project(obj) <- "YuImmunity2024"
  obj <- subset(obj, nCount_Spatial > 0)
  obj$pct.mito <- PercentageFeatureSet(obj, pattern = "^MT-")

  obj$study_id <- "YuImmunity2024"
  obj$sample_id <- paste0("YuImmunity2024_", sample)
  obj$image_id <- obj$sample_id

  return(obj)
})
names(seurat_list) <- sample_ids

library(loupeR)

s9 <- seurat_list[[9]]
spatial_coords <- s9

coords <- s9@images$YuImmunity2024_S09$centroids@coords
rownames(coords) <- s9@images$YuImmunity2024_S09$centroids@cells
colnames(coords) <- c("coords_1", "coords_2")
s9[['coords']] <- CreateDimReducObject(embeddings = coords, assay="Spatial")
loupeR::create_loupe_from_seurat(s9, output_dir = "data/published/Visium/yu-immunity-2024/S09/", output_name = 'S09')


s9.meta <- read.csv("data/published/Visium/yu-immunity-2024/S09/sample_id.assignments.csv")
rownames(s9.meta) <- s9.meta$Barcode
s9.meta$Barcode <- NULL
s9@meta.data$image_id <- s9@meta.data$sample_id
s9@meta.data$sample_id <- NULL
s9 <- AddMetaData(s9, metadata = s9.meta)
s9@meta.data$sample_id[nchar(s9@meta.data$sample_id)==0] <- "YuImmunity2024_S09A"
seurat_list[[9]] <- s9

YuImmunity2024 <- merge(seurat_list[[1]], seurat_list[2:length(seurat_list)], add.cell.ids = names(seurat_list))


sample.meta <- read.csv('data/published/Visium/yu-immunity-2024/sample_metadata.csv')
cell.meta <- YuImmunity2024@meta.data
cell.meta$spot_barcode <- rownames(YuImmunity2024@meta.data)
cell.meta <- cell.meta %>% left_join(sample.meta) %>% data.frame()
rownames(cell.meta) <- cell.meta$spot_barcode
cols <- setdiff(colnames(cell.meta), colnames(YuImmunity2024@meta.data))
YuImmunity2024 <- AddMetaData(YuImmunity2024, metadata = cell.meta[, cols])

qs::qsave(YuImmunity2024, file = 'data/visium/seurat_objects/YuImmunity2024.visium_data.raw.seurat_object.QS')


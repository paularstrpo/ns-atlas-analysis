library(Seurat)
library(SeuratObject)
library(Matrix)
library(tidyverse)
library(qs)
setwd("D:/Dropbox/_projects/NormalSkinAtlas")
## load data
obj <- qs::qread("data/merfish/BAYSOR/seurat_objects/ns-atlas.merfish_baysor.scanvi_integrated.cellcharter.annotated.seurat_object.QS")

obj$orig.x <- obj$center_x
obj$orig.y <- obj$center_y

meta.list <- split(obj@meta.data, f = obj@meta.data$sample_barcode)

rotation_angles <- read.csv("data/metadata/tissue_rotation.csv")
sample_barcodes <- names(meta.list)

new_coords <- lapply(names(meta.list), function(x){
    df <- meta.list[[x]][, c("center_x", "center_y")]
    angle <- rotation_angles$rotation_angle[rotation_angles$sample_barcode == x] * pi/180
    rot <- spdep::Rotation(df, angle)
    rot <- data.frame(rot)
    colnames(rot) <- c("x.rotated", "y.rotated")
    rot$cell_barcode <- rownames(rot)
    rot$sample_barcode <- x
    rot$tissue_compartment <- meta.list[[x]]$tissue_compartment
    return(rot)
}) %>%
    bind_rows() %>%
    data.frame()
rownames(new_coords) <- new_coords$cell_barcode

png("figures/merfish/BAYSOR/annotated/new_coords.png", height = 30, width=40, units="in", res=150)
ggplot(new_coords) +
    aes(x = x.rotated, y = y.rotated, color = tissue_compartment) +
    geom_point(size = 0.25, alpha=0.5) +
    facet_wrap(~ sample_barcode, scales = 'free', nrow=9) +
    scale_color_manual(values = c(epidermis="orchid1", dermis='goldenrod1', subcutis='dodgerblue')) +
    NoLegend() + NoAxes()
dev.off()


new_coords <- new_coords[, c("x.rotated", "y.rotated")]

obj <- AddMetaData(obj, new_coords)
obj$spatial_rotated_1 <- obj$x.rotated
obj$spatial_rotated_2 <- obj$y.rotated

spatial_coords <- as.matrix(obj@meta.data[, c("spatial_rotated_1", "spatial_rotated_2")])
obj[["spatial_rotated"]] <- CreateDimReducObject(spatial_coords, key = "spatialrotated_", assay = 'RNA')
levels(obj@meta.data$cell_type.detailed)[32] <- "Perineural Fib"

spatial_rotated <- data.frame(Embeddings(obj, reduction='spatial_rotated'))
write.csv(spatial_rotated, file='data/merfish/BAYSOR/metadata/rotated_tissue_coordinates.csv')

qs::qsave(obj, file = "data/merfish/BAYSOR/seurat_objects/ns-atlas.merfish_baysor.scanvi_integrated.cellcharter.annotated.seurat_object.QS")

obj <- subset(obj, cell_type.detailed != "Doublet" & cell_type.detailed != "Unknown" & tissue_compartment != "outside tissue")
obj <- UpdateSeuratObject(obj)

qs::qsave(obj, file = "data/merfish/BAYSOR/seurat_objects/ns-atlas.merfish_baysor.scanvi_integrated.cellcharter.annotated.filtered.seurat_object.QS")


cols <- c("cell_barcode", # unique cell barcode
          "center_x", "center_y", # x/y spatial coordinates
          "nFeature_RNA", "nCount_RNA", "volume", "gene_panel", # technical variables
          "batch", # batch id
          "sample_barcode", # unique tissue level ID
          "anatomic_site", # anatomic site
          "donor_id", # individual donor ID
          "donor_sex", # donor sex
          "donor_age", # donor age
          "scvi.leiden_1", # unsupervised cluster id
          "collection_source", # is autopsy, mohs discards or panniculectomy discards?
          "tissue_compartment", # manual tissue layer annotation
          "neighborhood", # neighborhood
          "cell_category", # overall population category
          "cell_type.broad", # l1 cell type
          "cell_type.detailed" # l2 cell type
          )
obj@meta.data <- obj@meta.data[, cols]
obj@meta.data$l1.cell_type <- obj@meta.data$cell_category
obj@meta.data$l2.cell_type <- obj@meta.data$cell_type.broad
obj@meta.data$l3.cell_type <- obj@meta.data$cell_type.detailed

obj@meta.data$UMAP_1 <- Embeddings(obj, reduction='scvi_umap')[,1]
obj@meta.data$UMAP_1 <- Embeddings(obj, reduction='scvi_umap')[,2]
obj@meta.data$spatial_1 <- obj@meta.data$center_x
obj@meta.data$spatial_2 <- obj@meta.data$center_y

obj2 <- DietSeurat(obj, assays = "RNA", layers = c("counts", "data"), dimreducs = c("scvi_umap", "scvi", "spatial", 'spatial_rotated'))
obj2 <- UpdateSeuratObject(obj2)

saveRDS(obj2, file = "data/website/ns-atlas.merfish.annotated.website.07-28-2025.RDS")
qs::qsave(obj2, "data/website/ns-atlas.merfish.annotated.website.07-28-2025.QS")

obj2 <- qs::qread("data/website/ns-atlas.merfish.annotated.website.07-28-2025.QS")


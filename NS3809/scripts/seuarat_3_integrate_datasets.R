# Load Required Libraries
library(Seurat)
library(future)
library(future.apply)
library(SoupX)
library(DropletUtils)
library(ggplot2)
library(DoubletFinder)
library(knitr)
library(dplyr)
library(patchwork)
library(reshape2)
library(rlang)
source("NS3809/scripts/utils/seurat_mod_integration.R")
source("NS3809/scripts/utils/seurat_mod_labels.R")

plan("multisession", workers = 4)  # Works with RStudio
options(future.globals.maxSize = 60 * 1024^3)

seurat_list <- readRDS( file="NS3809/results/seurat/objects/3a_soupx_seurats.Rds")

# Extract Seurat objects and additional outputs
seurat_list <- lapply(seurat_list, `[[`, "seurat_obj")
# Apply SCTransform and RunPCA to each Seurat object in the list using lapply
seurat_list <- lapply(seurat_list, function(x) {
  x <- SCTransform(x, vst.flavor = "v2", verbose = FALSE)
  return(x)
})

features <- SelectIntegrationFeatures(object.list = seurat_list, nfeatures = 3000)
seurat_list <- PrepSCTIntegration(object.list = seurat_list, anchor.features = features)
anchors <- FindIntegrationAnchors(object.list = seurat_list, normalization.method = "SCT", anchor.features = features)
seurat_integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
#saveRDS(seurat_integrated,file="NS3809/results/seurat/integrated_seurat.Rds")


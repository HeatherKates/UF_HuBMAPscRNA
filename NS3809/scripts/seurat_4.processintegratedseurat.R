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

plan("multisession", workers = 4)  # Works with RStudio
options(future.globals.maxSize = 60 * 1024^3)

seurat <- readRDS(file="NS3809/results/seurat/integrated_seurat.Rds")

seurat <- RunPCA(seurat, verbose = FALSE)
seurat <- RunUMAP(seurat, reduction = "pca", dims = 1:30, verbose = FALSE)
seurat <- FindNeighbors(seurat, reduction = "pca", dims = 1:30)
seurat <- FindClusters(seurat, resolution = 0.1)
#Really I should go back and add meta.data to initial datasets, but for now:
seurat@meta.data$dataset <- sapply(strsplit(colnames(seurat@assays[["integrated"]]@data), "_"), function(x) x[2])

DimPlot(seurat)
#Cluster 5 = Macrophage, Mast, Stellate1, Beta
#Cluster 7 = Endothelial
#Cluster 0,2,3,4 = Acinar
#Cluster 1 = Ductal
#Cluster 6 = Stellate 2

library(TissueEnrich)
library(clusterProfiler)
library(org.Hs.eg.db)
source("NS3809/scripts/utils/perform_enrichment.R")
# Example usage with a Seurat object
enrichment_results <- perform_enrichment_with_outputs(seurat)
# Combine top tissue enrichment (sorted by fold.change) and add Ident column
combined_tissue_enrichment <- bind_rows(lapply(names(enrichment_results), function(ident) {
  df <- enrichment_results[[ident]]$teEnrichmentOutput
  if (!is.null(df) && nrow(df) > 0) {
    df <- df[order(-df$Log10PValue), ][1, ]  # Take the top row based on fold.change
    df$Ident <- ident  # Add Ident column
    return(df)
  }
}))
combined_tissue_enrichment$pvalue <- 10^(-combined_tissue_enrichment$Log10PValue)
# Combine all cell type enrichment results and add Ident column
combined_celltype_enrichment <- bind_rows(lapply(names(enrichment_results), function(ident) {
  df <- enrichment_results[[ident]]$celltype_res_df
  if (!is.null(df) && nrow(df) > 0) {
    df$Ident <- ident  # Add Ident column
    return(df)
  }
}))
FeaturePlot(seurat,features = c("ELMO1","ROBO2","CADM1","CAPN13","GNAS","SUSD4","ITGA1","STAT3","ATP2A3","TERF2IP"))

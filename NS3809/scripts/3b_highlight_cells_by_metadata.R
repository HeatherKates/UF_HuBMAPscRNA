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
source("utils/assess_clusters.R")

plan("multisession", workers = 2)  # Works with RStudio
options(future.globals.maxSize = 60 * 1024^3)
seurat_merged <- readRDS(file="../results/seurat/objects/3a_soupx_seurat_merged.Rds")

########################################################
####### Plot marker genes from Elgamal et al. 2024 #####
########################################################

# Define the output directory
output_dir <- "../results/seurat/plots"

library(patchwork)
library(Seurat)

DimPlot(seurat_merged,split.by = "dataset")
##Cluster the cells
seurat_merged <- FindClusters(seurat_merged, resolution = 0.1)
#Save the current clusters
seurat_merged$seurat_clusters_0.1res <- Idents(seurat_merged)
#Add cluster for plotting
seurat_merged@meta.data$seurat_clusters <- paste0("cluster_",seurat_merged@meta.data$seurat_clusters)

# Example usage with seurat_merged where Idents() has been changed to 'dataset'
Idents(seurat_merged) <- seurat_merged@meta.data$seurat_clusters
#Generate the plots
lowrescluster_plots <- plot_highlight_idents(seurat_merged)
# Save the generated ident plots to a multi-page PDF
save_idents_plots_pdf(lowrescluster_plots, file_name = paste0(output_dir,"/lowres_cluster_highlight_plots.pdf"),2)
#Plot a heatmap of dataset x cluster enrichment
enrichment_heatmap(seurat_merged, "seurat_clusters", "dataset")



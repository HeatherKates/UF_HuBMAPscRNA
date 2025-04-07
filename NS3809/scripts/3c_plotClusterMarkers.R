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

seurat_merged <- PrepSCTFindMarkers(seurat_merged)
ClusterMarkers <- FindAllMarkers(seurat_merged,only.pos = TRUE,logfc.threshold = 0.2,min.diff.pct = .40)
# Filter ClusterMarkers for unique markers with strong specificity
unique_markers <- ClusterMarkers %>%
  dplyr::filter(avg_log2FC > 1,  # High fold change
                pct.1 > 0.25,    # Expressed in at least 25% of the target cluster
                pct.2 < 0.1,     # Expressed in less than 10% of other clusters
                p_val_adj < 0.05)  # Significant markers
# Load required libraries
library(org.Hs.eg.db)
library(clusterProfiler)
library(AnnotationDbi)
library(gridExtra)
library(rPanglaoDB)

relevant_tissues <- c("Pancreas", "Islets", "Endocrine pancreas", "Exocrine pancreas")
relevant_cell_types <- c(
  "Alpha cells", "Beta cells", "Delta cells", "Gamma cells", "Epsilon cells",
  "Acinar cells", "Ductal cells", "Stellate cells", "Endothelial cells",
  "Macrophages", "T cells", "B cells",
  "acinar",
  "ductal", 
  "endocrine", 
  "endothelial", 
  "erythroblast", 
  "immune", 
  "mesenchymal",
  "Schwann cells",
  "mast"
)

library(purrr)
library(dplyr)
library(stringr)
library(rPanglaoDB)
# Filter genes with high pct.1 in the cluster of interest and low pct.2 in other clusters
filter_markers <- function(markers_df, pct_threshold = 0.7, pct_diff = 0.4) {
  markers_df %>%
    filter(pct.1 > pct_threshold, (pct.1 - pct.2) > pct_diff)
}

filtered_markers <- filter_markers(ClusterMarkers)

# Query PanglaoDB for multiple genes with filtering on Tissue
result <- lapply(filtered_markers$gene[c(1:13,15:72)], function(gene) {
  res <- rPanglaoDB::getMarkers(include = gene)
  # Subset rows where "Tissue" contains "panc" (case insensitive)
  res <- res[grep("panc", res$Tissue, ignore.case = TRUE), ]
  res
})

# Combine the filtered results into a single dataframe
result_df <- do.call(rbind, result)

cluster_markers_with_anno <- merge(result_df,filtered_markers,by.x="Markers",by.y="gene")


library(gridExtra)

# Combine cell types for each Markers x cluster
cluster_markers_with_anno <- cluster_markers_with_anno %>%
  group_by(cluster, Markers) %>%
  summarise(
    Cell_Type = paste(unique(`Cell-Type`), collapse = "; "),
    avg_log2FC = first(avg_log2FC),  # Retain other metadata (e.g., avg_log2FC) if needed
    .groups = "drop"
  )

library(ggplot2)
library(gridExtra)
library(grid)

# Function to create FeaturePlots for a cluster
plot_cluster_markers <- function(seurat_obj, cluster, markers_df) {
  # Filter markers for the given cluster
  cluster_markers <- markers_df %>% filter(cluster == !!cluster)
  markers <- cluster_markers$Markers
  
  # Paginate markers if there are more than 9
  plots <- list()
  for (page_start in seq(1, length(markers), by = 9)) {
    page_markers <- markers[page_start:min(page_start + 8, length(markers))]
    page_data <- cluster_markers %>% filter(Markers %in% page_markers)
    
    # Generate plots for the page
    marker_plots <- lapply(seq_len(nrow(page_data)), function(i) {
      marker <- page_data$Markers[i]
      cell_type <- page_data$Cell_Type[i]
      # Break the title into lines (max 3 cell types per line)
      title <- paste(cluster, marker, paste(strwrap(cell_type, width = 30), collapse = "\n"))
      
      FeaturePlot(seurat_obj, features = marker, reduction = "umap", pt.size = 0.2) +  # Reduce point size
        ggtitle(title) +
        theme_minimal(base_size = 8)  # Adjust font size
    })
    
    # Combine plots into a grid
    plots[[length(plots) + 1]] <- grid.arrange(grobs = marker_plots, nrow = 3, ncol = 3,
                                               top = textGrob(paste(cluster, "Markers"),
                                                              gp = gpar(fontsize = 14, fontface = "bold")))
  }
  return(plots)
}

# Define the output directory
output_dir <- "cluster_png_plots/"
if (!dir.exists(output_dir)) dir.create(output_dir)  # Create the directory if it doesn't exist

# Loop over all unique clusters
for (cluster in unique(cluster_markers_with_anno$cluster)) {
  message("Processing cluster: ", cluster)
  
  # Generate plots for the current cluster
  cluster_plots <- plot_cluster_markers(seurat_merged, cluster, cluster_markers_with_anno)
  
  # Save each page of cluster_plots as a PNG
  for (i in seq_along(cluster_plots)) {
    # Define the output PNG file name
    output_png_file <- file.path(output_dir, paste0("cluster_", cluster, "_page_", i, ".png"))
    
    # Open a PNG graphics device
    png(filename = output_png_file, width = 1200, height = 1200, res = 150)  # Adjust resolution and size
    
    # Draw the plot
    grid.draw(cluster_plots[[i]])
    
    # Close the PNG graphics device
    dev.off()
    
    message("Saved: ", output_png_file)
  }
}

message("PNG files for all clusters saved to: ", output_dir)

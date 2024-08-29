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

# Define the base path where your data is located
base_path <- "/home/hkates/blue_garrett/Campbell-Thompson/scRNA/NS3809/results/cellranger/"

# Define the output path for SoupX-corrected data
soupx_output_path <- "/blue/timgarrett/hkates/Campbell-Thompson/scRNA/NS3809/results/soupx/"

# List all directories matching the pattern
sample_dirs <- Sys.glob(paste0(base_path, "NS*"))
options(future.globals.maxSize = 40 * 1024^3)

# Function to load, correct with SoupX, filter, and normalize each dataset
load_correct_and_preprocess_sample <- function(dir) {
  # Define paths for raw and filtered matrices
  raw_data_path <- file.path(dir, "outs", "raw_feature_bc_matrix.h5")
  filtered_data_path <- file.path(dir, "outs", "filtered_feature_bc_matrix.h5")
  
  # Load raw and filtered data
  raw.matrix  <- Read10X_h5(raw_data_path, use.names = TRUE)
  filt.matrix <- Read10X_h5(filtered_data_path, use.names = TRUE)
  sample_name <- strsplit(dir, "/")[[1]][[10]] # Extract sample name from directory path
  # Create Seurat object from filtered data for clustering
  srat  <- CreateSeuratObject(counts = filt.matrix)
  srat <- SCTransform(srat, verbose = FALSE)
  srat <- RunPCA(srat, verbose = TRUE, npcs = 30)
  srat <- RunUMAP(srat, dims = 1:30, verbose = FALSE)
  srat <- FindNeighbors(srat, dims = 1:30, verbose = FALSE)
  srat <- FindClusters(srat, verbose = TRUE, resolution = 1)
  
  # Prepare SoupChannel object using raw and filtered data
  meta    <- srat@meta.data
  umap    <- srat@reductions$umap@cell.embeddings
  soup.channel  <- SoupChannel(raw.matrix, filt.matrix)
  soup.channel  <- setClusters(soup.channel, setNames(meta$seurat_clusters, rownames(meta)))
  soup.channel  <- setDR(soup.channel, umap)
  
  # Try to estimate contamination and adjust counts
  tryCatch({
    # Redirect the plot output to a file before running autoEstCont
    pdf(file = paste0(soupx_output_path,"plots/", sample_name, "_contamination_plot.v2.pdf"))
    soup.channel  <- autoEstCont(soup.channel)
    dev.off()  # Close the PDF device
    
    # Save the plot from plotMarkerDistribution
    pdf(file = paste0(soupx_output_path,"plots/", sample_name, "_marker_distribution_plot.v2.pdf"))
    plotMarkerDistribution(soup.channel)
    dev.off()  # Close the PDF device
    
    # Identify genes highly expressed in the background
    # Extract the soup profile and sort by estimated contamination
    high_background_genes <- head(soup.channel[["soupProfile"]][order(-soup.channel[["soupProfile"]][, "est"]), ], 100)
    
    # Adjust the raw counts to remove the estimated contamination
    corrected_data <- adjustCounts(soup.channel, roundToInt = TRUE)
    
    # Save the corrected data
    saveRDS(corrected_data, file = paste0(soupx_output_path,"objects/", sample_name, "_soupx_corrected.v2.Rds"))
    
    # Create Seurat object from corrected data
    seurat_obj <- CreateSeuratObject(counts = corrected_data, project = sample_name, min.cells = 3, min.features = 200)
    
    # Calculate the percentage of mitochondrial genes
    seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
    seurat_obj[["sample"]] <- sample_name
    
    # Filter cells based on quality metrics
    seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
    
    # Process with SCTransform and clustering
    seurat_obj <- SCTransform(seurat_obj, verbose = FALSE)
    seurat_obj <- RunPCA(seurat_obj, verbose = TRUE, npcs = 30)
    seurat_obj <- RunUMAP(seurat_obj, dims = 1:30, verbose = FALSE)
    seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30, verbose = FALSE)
    seurat_obj <- FindClusters(seurat_obj, verbose = TRUE, resolution = 0.3)
    
    # Return both Seurat object and the summary of ambient RNA correction
    return(list(
      seurat_obj = seurat_obj,
      high_background_genes = high_background_genes
    ))
  }, error = function(e) {
    message("Error in processing ", sample_name, ": ", e$message)
    return(NULL)  # Return NULL if an error occurs
  })
}

# Parallelize loading, correcting, and preprocessing the datasets
plan(multicore, workers = 8)  # Adjust the number of workers based on your system
results_list <- future_lapply(sample_dirs, load_correct_and_preprocess_sample, future.seed = TRUE)

# Filter out NULL results due to errors
results_list <- Filter(Negate(is.null), results_list)

saveRDS(results_list, file="/home/hkates/blue_garrett/Campbell-Thompson/scRNA/NS3809/results/seurat/objects/3a_soupx_seurats.v2.Rds")

# Extract Seurat objects and additional outputs
seurat_list <- lapply(results_list, `[[`, "seurat_obj")
high_background_genes_list <- lapply(results_list, `[[`, "high_background_genes")
rm(results_list)

# Merge seurats
seurat_merged <- merge(seurat_list[[1]], y = c(seurat_list[2:13]), project = "HuBMAP_panc")
seurat_merged$percent.mt <- PercentageFeatureSet(seurat_merged, pattern = "^MT-")
seurat_merged <- SCTransform(seurat_merged, vars.to.regress = "percent.mt", verbose = FALSE,conserve.memory=TRUE,variable.features.n = 3000)
seurat_merged <- RunPCA(seurat_merged, verbose = FALSE)
seurat_merged <- RunUMAP(seurat_merged, dims = 1:30, verbose = FALSE)
seurat_merged <- FindNeighbors(seurat_merged, dims = 1:30, verbose = FALSE)
seurat_merged <- FindClusters(seurat_merged, verbose = FALSE, resolution = 0.3,future.seed=TRUE)
saveRDS(seurat_merged, file="/home/hkates/blue_garrett/Campbell-Thompson/scRNA/NS3809/results/seurat/objects/3a_soupx_seurat_merged.v2.Rds")

########################################################
####### Plot marker genes from Elgamal et al. 2024 #####
########################################################

#cell type marker genes
markers <- read.table("NS3809/scripts/References/marker_genes.csv",header=FALSE)
colnames(markers) <- c("CellType","Feature")
top_variable_genes <- seurat_merged@assays[["SCT"]]@var.features
# To get the top 500 variable genes
top_500_variable_genes <- head(top_variable_genes, 500)

# Filter the markers to keep only those in the top 500 variable genes
variable_markers <- markers %>% filter(Feature %in% top_500_variable_genes)

# Group by CellType and select the top 9 markers while keeping the original order
top_markers_by_cell_type <- variable_markers %>%
  group_by(CellType) %>%
  slice_head(n = 9)
top_markers_by_cell_type$CellType <- gsub("$"," Cells",top_markers_by_cell_type$CellType)

#############################################
####### Plot established marker genes #######
#############################################
#cell type marker genes
est_markers <- read.table("NS3809/scripts/References/Established_markers.txt",header=TRUE,sep="\t")
est_markers <- as.data.frame(cbind(est_markers[,2],est_markers[,1]))
colnames(est_markers) <- c("CellType","Feature")

# Filter the markers to keep only those in the top 500 variable genes
variable_est_markers <- est_markers[est_markers$Feature %in% seurat_merged@assays[["SCT"]]@var.features,]

# Define the output directory
output_dir <- "/home/hkates/blue_garrett/Campbell-Thompson/scRNA/NS3809/results/seurat/plots"

library(patchwork)

# Define the function
generate_feature_plots <- function(seurat_obj, features_df, features_name, output_dir) {
  # De-frame the features_df to use its name in titles and filenames
  features_df_name <- deparse(substitute(features_df))
  
  # Initialize an empty list to store the plots
  est_plots <- list()
  
  # Loop through each CellType, generate FeaturePlots, save them, and store them in the list
  for(cell_type in unique(features_df$CellType)) {
    features <- features_df %>% filter(CellType == cell_type) %>% pull(Feature)
    
    # Dynamically set ncol based on the number of features
    n_features <- length(features)
    if (n_features < 3) {
      ncol <- n_features
    } else if (n_features == 4) {
      ncol <- 2
    } else {
      ncol <- 3
    }
    
    # Generate the plot list to allow for a shared title
    plot_list <- FeaturePlot(
      seurat_obj,
      features = features,
      combine = FALSE,
      max.cutoff = 'q98'
    )
    
    # Combine the plots using patchwork and add the title with features_df name
    combined_plot <- wrap_plots(plot_list, ncol = ncol) + 
      plot_annotation(title = paste(cell_type, features_name,sep=" "))
    
    # Define the output file path with features_df name included
    output_file <- file.path(output_dir, paste0(cell_type, "_", gsub(" ","_",features_name), "_FeaturePlot.v2.png"))
    
    # Save the plot as a high-resolution PNG file
    ggsave(output_file, plot = combined_plot, width = 10, height = 10, units = "in", dpi = 300)
    
    # Store the combined plot in the list with the features_df name
    est_plots[[paste0(cell_type, "_", features_name)]] <- combined_plot
  }
  
  # Return the list of plots
  return(est_plots)
}

# Generate plots for established markers
est_plots <- generate_feature_plots(seurat_merged, variable_est_markers,"variable established markers","/home/hkates/blue_garrett/Campbell-Thompson/scRNA/NS3809/results/seurat/plots")
# Generate plots for additional markers
addtl_feature_plots <- generate_feature_plots(seurat_merged,top_markers_by_cell_type ,"variable Elgamal et al. top cluster markers","/home/hkates/blue_garrett/Campbell-Thompson/scRNA/NS3809/results/seurat/plots")

clusterplot <- DimPlot(seurat_merged)+ggtitle ("Seurat Clusters (98,946 cells from 13 merged samples)")
ggsave(file.path(output_dir, paste0("ClusterPlot.v2.png")), plot = clusterplot, dpi = 300, width = 10, height = 8)

samplesplot <- DimPlot(seurat_merged,group.by="sample")+ggtitle ("Samples (98,946 cells from 13 merged donor samples)")
ggsave(file.path(output_dir, paste0("SamplesPlot.v2.png")), plot = samplesplot, dpi = 300, width = 10, height = 8)

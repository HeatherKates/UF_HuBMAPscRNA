# Load Required Libraries
library(Seurat)
library(future)
library(future.apply)

# Define the base path where your data is located
base_path <- "/home/hkates/blue_garrett/Campbell-Thompson/scRNA/NS3809/results/cellranger/"

# List all directories matching the pattern
sample_dirs <- list.dirs(base_path, full.names = TRUE, recursive = TRUE)
sample_dirs <- sample_dirs[grep("NS-3809_.*/outs/filtered_feature_bc_matrix$", sample_dirs)]

# Function to load, filter, and normalize each dataset
load_and_preprocess_sample <- function(dir) {
  data <- Read10X(data.dir = dir)
  sample_name <- basename(dirname(dirname(dir)))  # Extract sample name from directory path
  seurat_obj <- CreateSeuratObject(counts = data, project = sample_name, min.cells = 3, min.features = 200)
  
  # Calculate the percentage of mitochondrial genes
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  
  # Filter cells based on quality metrics
  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
  
  # Normalize the data
  seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
  
  # Find variable features
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
  
  return(seurat_obj)
}

# Parallelize loading and preprocessing the datasets
plan(multicore, workers = 8)  # Adjust the number of workers based on your system
seurat_list <- future_lapply(sample_dirs, load_and_preprocess_sample, future.seed = TRUE)

# Merge all Seurat objects into one
combined_seurat <- Reduce(function(x, y) merge(x, y), seurat_list)

# Save the combined Seurat object
saveRDS(combined_seurat, file="/home/hkates/blue_garrett/Campbell-Thompson/scRNA/NS3809/results/seurat/combined_seurat.Rds")


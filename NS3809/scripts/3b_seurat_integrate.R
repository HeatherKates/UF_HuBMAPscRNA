library(Seurat)
library(future.apply)

# Increase the future.globals.maxSize limit (set to 8 GiB here, adjust if needed)
options(future.globals.maxSize = 8 * 1024^3)  # 8 GiB

# Set up parallel processing
plan(multicore, workers = 8)

# Load the combined Seurat object
combined_seurat <- readRDS("/home/hkates/blue_garrett/Campbell-Thompson/scRNA/NS3809/results/seurat/combined_seurat.Rds")

# Split the combined Seurat object by the original identity (sample) labels
seurat_list_split <- SplitObject(combined_seurat, split.by = "orig.ident")

# Select the features to use for integration
features <- SelectIntegrationFeatures(object.list = seurat_list_split)

# Parallelize finding integration anchors
anchors <- FindIntegrationAnchors(object.list = seurat_list_split, anchor.features = features, 
                                  k.anchor = 10)  # Increase k.anchor if necessary

# Parallelize integrating the data
integrated_seurat <- IntegrateData(anchorset = anchors, k.weight = 20)  # Lower k.weight

# Save the integrated Seurat object
saveRDS(integrated_seurat, file="/home/hkates/blue_garrett/Campbell-Thompson/scRNA/NS3809/results/seurat/3b_integrated_seurat.Rds")

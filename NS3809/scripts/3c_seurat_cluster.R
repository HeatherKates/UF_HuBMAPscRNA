library(Seurat)
library(future.apply)
#read data
integrated_seurat=readRDS("/home/hkates/blue_garrett/Campbell-Thompson/scRNA/NS3809/results/seurat/3b_integrated_seurat.Rds")

# Scaling, PCA, and Clustering
integrated_seurat <- ScaleData(integrated_seurat)
integrated_seurat <- RunPCA(integrated_seurat, npcs = 30)
integrated_seurat <- RunUMAP(integrated_seurat, dims = 1:30)
integrated_seurat <- FindNeighbors(integrated_seurat, dims = 1:30)
integrated_seurat <- FindClusters(integrated_seurat, resolution = 0.5)

# Visualize the integrated data
p=DimPlot(integrated_seurat, reduction = "umap", label = TRUE)

saveRDS(file="/home/hkates/blue_garrett/Campbell-Thompson/scRNA/NS3809/results/seurat/3c_clustered.integrated_seurat.Rds",integrated_seurat)
saveRDS(file="/home/hkates/blue_garrett/Campbell-Thompson/scRNA/NS3809/results/seurat/3c_clustered.dimplot.Rds",p)


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

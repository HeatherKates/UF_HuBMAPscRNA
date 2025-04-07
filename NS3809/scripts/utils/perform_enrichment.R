
  # Function to perform enrichment and save additional outputs for each cluster
  perform_enrichment_with_outputs <- function(seurat_obj) {
    all_idents <- unique(Idents(seurat_obj))  # Get all cluster IDs
    
    enrichment_results <- list()  # Store all results
    
    # Loop through each ident (cluster)
    for (ident in all_idents) {
      
      # Identify cluster markers
      cluster_markers <- FindMarkers(seurat_obj, ident.1 = ident, pos.only = TRUE) %>%
        filter(p_val_adj < 0.05) %>%  # Filter significant markers
        filter(avg_log2FC > 0)  # Filter for upregulated genes
      
      if (nrow(cluster_markers) == 0) {
        next  # Skip if no significant markers
      }
      
      gene_list <- rownames(cluster_markers)
      
      # Create a list to store the three dataframes for this cluster/ident
      cluster_results <- list()
      
      # Store the filtered cluster markers
      cluster_results$cluster_markers <- cluster_markers
      
      # Perform tissue enrichment using TissueEnrich
      if (length(gene_list) > 0) {
        my.gs <- GeneSet(geneIds = gene_list, organism = "Homo Sapiens", geneIdType = SymbolIdentifier()) 
        my.output <- teEnrichment(inputGenes = my.gs)
        
        # Convert tissue enrichment result to dataframe before filtering
        teEnrichmentOutput <- as.data.frame(my.output[[1]]@assays@data@listData[[1]])
        
        # Store the raw tissue enrichment output before filtering
        cluster_results$teEnrichmentOutput <- teEnrichmentOutput
      }
      
      # Perform cell type enrichment using rcellmarker
      if (length(gene_list) > 0) {
        celltype_res <- cells(gene_list, species = "human", keytype = "SYMBOL")
        celltype_res_df <- celltype_res@result
        
        # Store the raw cell type enrichment output before filtering
        cluster_results$celltype_res_df <- celltype_res_df
      }
      
      # Save the results for this ident (cluster)
      enrichment_results[[as.character(ident)]] <- cluster_results
    }
    
    return(enrichment_results)
  }



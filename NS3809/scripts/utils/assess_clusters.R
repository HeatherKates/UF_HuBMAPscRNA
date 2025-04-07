#Example usage
#lowrescluster_plots <- plot_highlight_idents(seurat_merged)
## Save the generated ident plots to a multi-page PDF
#save_idents_plots_pdf(lowrescluster_plots, file_name = paste0(output_dir,"/lowres_cluster_highlight_plots.pdf"),2)
##Plot a heatmap of dataset x cluster enrichment
#enrichment_heatmap(seurat_merged, "seurat_clusters", "dataset")

# Function to dynamically generate side-by-side plots for all levels of the current Idents() in the Seurat object
plot_highlight_idents <- function(seurat_obj) {
  # Get current identities and their levels
  cluster_ids <- as.factor(Idents(seurat_obj))  # Explicitly convert to factor to handle levels properly
  ident_levels <- levels(cluster_ids)  # Explicitly get the levels of the identities
  
  # Dynamically retrieve the metadata column name corresponding to the current Idents()
  ident_name <- colnames(seurat_obj@meta.data)[grep(pattern = levels(cluster_ids)[[1]], x = seurat_obj@meta.data)]
  
  # Fallback in case ident_name is not found
  ident_name <- ifelse(length(ident_name) == 0, "Ident", ident_name)
  
  # Define colors used by Seurat (default color palette)
  ident_colors <- scales::hue_pal()(length(ident_levels))
  
  # Create a list to store the combined plots
  combined_plots <- list()
  
  # Generate the original DimPlot (without highlighting)
  p1 <- DimPlot(seurat_obj)
  
  # Loop through each level of the current ident (e.g., each sample or cluster)
  for (ident_level in ident_levels) {
    # Get the color for the current ident level
    ident_color <- ident_colors[which(ident_levels == ident_level)]
    
    # Dynamically create a named list for cells.highlight with the ident level name
    highlight_cells <- setNames(list(WhichCells(seurat_obj, idents = ident_level)), paste0(ident_name, ": ", ident_level))
    
    # Generate the second plot, highlighting the current ident level with dynamic labels
    p2 <- DimPlot(seurat_obj, 
                  cells.highlight = highlight_cells, 
                  cols.highlight = ident_color)
    
    # Combine the two plots side by side
    combined_plot <- p1 | p2
    
    # Store the combined plot in the list, named by the ident level
    combined_plots[[paste0(ident_name, "_", ident_level)]] <- combined_plot
  }
  
  # Return the list of combined plots
  return(combined_plots)
}

# Function to plot multiple ident levels per page and save to a multi-page PDF
save_idents_plots_pdf <- function(ident_plots, file_name = "ident_highlight_plots.pdf", plots_per_page = 2) {
  # Split the list of plots into chunks of 'plots_per_page'
  plot_chunks <- split(ident_plots, ceiling(seq_along(ident_plots) / plots_per_page))
  
  # Create a PDF file
  pdf(file = file_name, width = 10, height = 12)  # Adjust the size as needed
  
  # Loop through each chunk (group of plots)
  for (chunk in plot_chunks) {
    # Combine the plots for the current chunk using patchwork (2x2 grid)
    combined_plot <- wrap_plots(chunk, ncol = 1, nrow = 2)
    
    # Print the combined plot to the PDF
    print(combined_plot)
  }
  
  # Close the PDF file
  dev.off()
}
enrichment_heatmap <- function(seurat_obj, var1, var2) {
  # Step 1: Create a contingency table summarizing the distribution of var1 (e.g., clusters) across var2 (e.g., datasets)
  distribution_table <- table(seurat_obj@meta.data[[var1]], seurat_obj@meta.data[[var2]])
  
  # Step 2: Calculate total proportion of cells for each var2 (e.g., dataset) across the entire Seurat object
  total_cells <- nrow(seurat_obj@meta.data)
  total_cells_per_var2 <- table(seurat_obj@meta.data[[var2]])
  proportion_total_var2 <- total_cells_per_var2 / total_cells
  
  # Step 3: Calculate the proportion of cells for each var2 in each var1 (e.g., dataset in each cluster)
  total_cells_per_var1 <- rowSums(distribution_table)
  proportion_per_var1 <- sweep(distribution_table, 1, total_cells_per_var1, "/")
  
  # Step 4: Compute enrichment factor (proportion in var1 / overall proportion)
  enrichment_factor <- sweep(proportion_per_var1, 2, proportion_total_var2, "/")
  
  # Step 5: Convert the enrichment_factor table to a data frame for plotting
  enrichment_df <- as.data.frame.matrix(enrichment_factor)
  enrichment_df[[var1]] <- rownames(enrichment_df)
  
  # Melt the data frame for ggplot (long format)
  enrichment_melt <- melt(enrichment_df, id.vars = var1, variable.name = var2, value.name = "enrichment")
  
  # Step 6: Visualize the enrichment factor as a heatmap
  ggplot(enrichment_melt, aes_string(x = var1, y = var2, fill = "enrichment")) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 1, 
                         limit = c(0, max(enrichment_melt$enrichment, na.rm = TRUE)), 
                         space = "Lab", name = "Enrichment") +
    theme_minimal() +
    labs(title = paste("Enrichment of", var2, "in", var1),
         x = var1,
         y = var2,
         fill = "Enrichment Factor") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

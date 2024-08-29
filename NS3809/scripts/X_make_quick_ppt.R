library(officer)
library(magrittr)

# Define the directory containing the PNG files
png_dir <- "/home/hkates/blue_garrett/Campbell-Thompson/scRNA/NS3809/results/seurat/plots"

# List all PNG files in the directory
png_files <-list.files(png_dir, pattern = "\\.png$", full.names = TRUE)

## Create a new PowerPoint object
ppt <- read_pptx()

# Add a title slide with the title text
ppt <- add_slide(ppt, layout = "Title and Content", master = "Office Theme") %>%
  ph_with(value = fpar(ftext("Initial analysis of snRNA from 16 samples across six donors",fp_text(bold = TRUE, font.size = 24))), location = ph_location_type(type = "title"))

# Add the bullet points
#bullet_points <- c(
#  "1. cellranger count performs alignment, filtering, barcode counting, and UMI counting\n2. SoupX to remove ambient RNA contamination\n3. Seurat v5 processing on corrected counts\n     1. Per sample: QC > Normalize > PCA > UMAP > Compute nearest neighbors > Identify Clusters\n     2. Combine all samples > Normalize > PCA > UMAP > Compute nearest neighbors > Identify Clusters\n4. Assess batch effect (no obvious batch effect, no integration)\n5. Plot expression of cell type markers in UMAP space\n     1. Established markers (as listed in Elgamal et al. 2023)\n     2. Top cluster markers from Elgamal et al. 2023\n")
#ppt <- ph_with(ppt,value=fpar(ftext(bullet_points,fp_text(bold = FALSE, font.size = 20))),location=ph_location_type(type="body"))
  
# Add hyperlinks

# Loop through each PNG file and add it as a slide
for (png_file in png_files) {
  ppt <- add_slide(ppt, layout = "Blank", master = "Office Theme") %>%
    ph_with(value = external_img(png_file), location = ph_location_fullsize())
}

# Define the output PowerPoint file path
output_pptx <- paste0(png_dir,"/seurat_presentation.pptx")

# Save the PowerPoint presentation
print(ppt, target = output_pptx)

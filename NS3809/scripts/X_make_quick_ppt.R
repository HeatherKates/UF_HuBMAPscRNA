library(officer)
library(magrittr)

# Define the directory containing the PNG files
png_dir <- paste0("results/plots/",dependent_variable,"/")

# List all PNG files in the directory
png_files <-list.files(png_dir, pattern = "\\.png$", full.names = TRUE)

## Create a new PowerPoint object
ppt <- read_pptx()

# Add a title slide with the title text
ppt <- add_slide(ppt, layout = "Title and Content", master = "Office Theme") %>%
  ph_with(value = fpar(ftext(paste("Effect of predictor variables on",dependent_variable),fp_text(bold = TRUE, font.size = 24))), location = ph_location_type(type = "title"))

# Loop through each PNG file and add it as a slide
for (png_file in png_files) {
  ppt <- add_slide(ppt, layout = "Blank", master = "Office Theme") %>%
    ph_with(value = external_img(png_file), location = ph_location_fullsize())
}

# Define the output PowerPoint file path
output_pptx <- paste0(png_dir,dependent_variable_plots,".pptx")

# Save the PowerPoint presentation
print(ppt, target = output_pptx)

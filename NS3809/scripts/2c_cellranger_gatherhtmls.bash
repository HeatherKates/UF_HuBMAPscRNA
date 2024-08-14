#!/bin/bash

# Define the base path where your data is located
base_path="/home/hkates/blue_garrett/Campbell-Thompson/scRNA/NS3809/results/cellranger/"

# Define the target directory where renamed files will be stored
target_dir="/home/hkates/blue_garrett/Campbell-Thompson/scRNA/NS3809/results/web_summary_htmls"

# Create the target directory if it doesn't exist
mkdir -p "$target_dir"

# Find all web_summary.html files in the specified directories
for sample_dir in "$base_path"NS-3809_*; do
    # Define the full path to the web_summary.html file
    web_summary_file="$sample_dir/outs/web_summary.html"
    
    # Check if the web_summary.html file exists
    if [ -f "$web_summary_file" ]; then
        # Get the name of the sample directory (e.g., NS-3809_10_P4_1)
        sample_name=$(basename "$sample_dir")
        
        # Define the new filename with the sample name included
        new_filename="${sample_name}_web_summary.html"
        
        # Copy and rename the file to the target directory
        cp "$web_summary_file" "$target_dir/$new_filename"
    fi
done

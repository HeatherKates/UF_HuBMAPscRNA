#!/bin/bash

# Directory containing the FastQC zip files
FASTQC_DIR="/blue/zhangw/hkates/Tanzia_RNAseq/results/fastqc"

# Output CSV file
OUTPUT_CSV="${FASTQC_DIR}/fastqc_summary.csv"

# Header for the CSV file
echo "Sample,Total Sequences,Sequences flagged as poor quality,Sequence length,GC content,Basic Statistics,Per base sequence quality,Per tile sequence quality,Per sequence quality scores,Per base sequence content,Per sequence GC content,Per base N content,Sequence Length Distribution,Sequence Duplication Levels,Overrepresented sequences,Adapter Content" > $OUTPUT_CSV

# Loop through each FastQC zip file
for zip_file in ${FASTQC_DIR}/*.zip
do
  # Extract the sample name from the zip file name
  sample_name=$(basename $zip_file | sed 's/_fastqc.zip//')

  # Unzip the fastqc_data.txt file
  unzip -p $zip_file */fastqc_data.txt > temp_fastqc_data.txt
  unzip -p $zip_file */summary.txt > temp_summary.txt

  # Extract relevant metrics from fastqc_data.txt
  total_sequences=$(grep "Total Sequences" temp_fastqc_data.txt | awk '{print $3}')
  poor_quality=$(grep "Sequences flagged as poor quality" temp_fastqc_data.txt | awk '{print $6}')
  sequence_length=$(grep "Sequence length" temp_fastqc_data.txt | awk '{print $3}')
  gc_content=$(grep "%GC" temp_fastqc_data.txt | awk '{print $2}')

  # Extract summary.txt metrics
  basic_statistics=$(grep "Basic Statistics" temp_summary.txt | awk '{print $1}')
  per_base_quality=$(grep "Per base sequence quality" temp_summary.txt | awk '{print $1}')
  per_tile_quality=$(grep "Per tile sequence quality" temp_summary.txt | awk '{print $1}')
  per_seq_quality=$(grep "Per sequence quality scores" temp_summary.txt | awk '{print $1}')
  per_base_content=$(grep "Per base sequence content" temp_summary.txt | awk '{print $1}')
  per_seq_gc=$(grep "Per sequence GC content" temp_summary.txt | awk '{print $1}')
  per_base_n=$(grep "Per base N content" temp_summary.txt | awk '{print $1}')
  seq_length_dist=$(grep "Sequence Length Distribution" temp_summary.txt | awk '{print $1}')
  seq_dup_levels=$(grep "Sequence Duplication Levels" temp_summary.txt | awk '{print $1}')
  overrep_seq=$(grep "Overrepresented sequences" temp_summary.txt | awk '{print $1}')
  adapter_content=$(grep "Adapter Content" temp_summary.txt | awk '{print $1}')

  # Append the metrics to the CSV file
  echo "$sample_name,$total_sequences,$poor_quality,$sequence_length,$gc_content,$basic_statistics,$per_base_quality,$per_tile_quality,$per_seq_quality,$per_base_content,$per_seq_gc,$per_base_n,$seq_length_dist,$seq_dup_levels,$overrep_seq,$adapter_content" >> $OUTPUT_CSV

  # Clean up temporary files
  rm temp_fastqc_data.txt
  rm temp_summary.txt
done

echo "Summary CSV file created at $OUTPUT_CSV"


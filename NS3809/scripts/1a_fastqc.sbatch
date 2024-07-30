#!/bin/bash
#SBATCH --job-name=fastqc.%j
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8GB
#SBATCH --qos=timgarrett-b
#SBATCH --account=timgarrett
#SBATCH -t 02:00:00
#SBATCH --output=/blue/timgarrett/hkates/Campbell-Thompson/scRNA/NS3809/logs/fastqc.%j.out
#SBATCH --error=/blue/timgarrett/hkates/Campbell-Thompson/scRNA/NS3809/logs/fastqc.%j.err
#SBATCH --array=1-16

echo "This is task $SLURM_ARRAY_TASK_ID"

echo -e "\nInfo: Starting a job on $(date) on $(hostname) in $(pwd).\n"

###################################################################
# This file is to perform FastQC quality control
###################################################################

module load fastqc/0.12.1

# Define the base directories for lanes
BASE_DIR1="/orange/timgarrett/Campbell_Thompson/NS3809/NS-3809_NS3810-MBrusko/NS-3809-MBrusko_10B-22K5WnLT3-Lane1-3"
BASE_DIR2="/orange/timgarrett/Campbell_Thompson/NS3809/NS-3809_NS3810-MBrusko/NS-3809-MBrusko_10B-22K5WnLT3-Lane4-6"

# Define the concatenated files directory
CONCAT_DIR="/orange/timgarrett/Campbell_Thompson/NS3809/Concatenated"

# Create the concatenated directory if it doesn't exist
mkdir -p ${CONCAT_DIR}

# Create an array of unique sample identifiers dynamically
SAMPLES=($(ls ${BASE_DIR1} | cut -d "_" -f 1,2,3 | sort | uniq))

# Get the sample for the current task
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID-1]}

# Define the output directory
OUT_DIR="/blue/timgarrett/hkates/Campbell-Thompson/scRNA/NS3809/results/fastqc"

# Concatenate files from all lanes
cat ${BASE_DIR1}/${SAMPLE}*/P*_R1_001.fastq.gz ${BASE_DIR2}/${SAMPLE}*/P*_R1_001.fastq.gz > ${CONCAT_DIR}/${SAMPLE}_combined_R1.fastq.gz
cat ${BASE_DIR1}/${SAMPLE}*/P*_R2_001.fastq.gz ${BASE_DIR2}/${SAMPLE}*/P*_R2_001.fastq.gz > ${CONCAT_DIR}/${SAMPLE}_combined_R2.fastq.gz

# Run FastQC with specified threads
fastqc --threads 2 ${CONCAT_DIR}/${SAMPLE}_combined_R1.fastq.gz ${CONCAT_DIR}/${SAMPLE}_combined_R2.fastq.gz -o ${OUT_DIR}

echo -e "\nInfo: FastQC job completed on $(date) on $(hostname).\n"
##End


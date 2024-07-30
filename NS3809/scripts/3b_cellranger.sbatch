#!/bin/bash
#SBATCH --job-name=cellranger.%j
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64GB
#SBATCH --qos=timgarrett-b
#SBATCH --account=timgarrett
#SBATCH -t 24:00:00
#SBATCH --output=/blue/timgarrett/hkates/Campbell-Thompson/scRNA/NS3809/logs/cellranger.%j.out
#SBATCH --error=/blue/timgarrett/hkates/Campbell-Thompson/scRNA/NS3809/logs/cellranger.%j.err
#SBATCH --array=1-16

echo "This is task $SLURM_ARRAY_TASK_ID"

echo -e "\nInfo: Starting a job on $(date) on $(hostname) in $(pwd).\n"

###################################################################
# This file is to perform quantification using Cell Ranger
###################################################################

module load cellranger/6.1.1

# Define the base directory
BASE_DIR="/blue/timgarrett/hkates/Campbell-Thompson/scRNA/NS3809/results/trimgalore"

# Create an array of unique sample identifiers dynamically
SAMPLES=($(ls | cut -d "_" -f 1,2,3 | sort | uniq))

# Get the sample for the current task
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID-1]}

# Define the output directory
OUT_DIR="/blue/timgarrett/hkates/Campbell-Thompson/scRNA/NS3809/results/cellranger"

# Define the path to the reference transcriptome
REF_DIR="/path/to/refdata-gex-GRCh38-2024-A"  # Update this path to the decompressed directory

# Run Cell Ranger
cellranger count --id=${SAMPLE} \
                 --fastqs=${BASE_DIR} \
                 --sample=${SAMPLE} \
                 --transcriptome=${REF_DIR} \
                 --localcores=8 \
                 --localmem=64 \
                 --output-dir=${OUT_DIR}

echo -e "\nInfo: Cell Ranger job completed on $(date) on $(hostname).\n"
##End


#!/bin/bash
#SBATCH --job-name=cellranger.%j
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64GB
#SBATCH --qos=timgarrett-b
#SBATCH --account=timgarrett
#SBATCH -t 18:00:00
#SBATCH --output=/blue/timgarrett/hkates/Campbell-Thompson/scRNA/NS3809/logs/cellranger.%j.out
#SBATCH --error=/blue/timgarrett/hkates/Campbell-Thompson/scRNA/NS3809/logs/cellranger.%j.err
#SBATCH --array=1-16

echo "This is task $SLURM_ARRAY_TASK_ID"

echo -e "\nInfo: Starting a job on $(date) on $(hostname) in $(pwd).\n"

###################################################################
# This file is to perform quantification using Cell Ranger
###################################################################

module load cellranger/8.0.0

# Define the base directories for raw FASTQ files
BASE_DIR1="/orange/timgarrett/Campbell_Thompson/NS3809/NS-3809_NS3810-MBrusko/NS-3809-MBrusko_10B-22K5WnLT3-Lane1-3"
BASE_DIR2="/orange/timgarrett/Campbell_Thompson/NS3809/NS-3809_NS3810-MBrusko/NS-3809-MBrusko_10B-22K5WnLT3-Lane4-6"

# Create an array of unique sample identifiers dynamically
SAMPLES=($(ls ${BASE_DIR1} | cut -d "_" -f 1,2,3,4 | sort | uniq))

# Ensure SLURM_ARRAY_TASK_ID starts at 1 (correct usage)
SAMPLE_INDEX=$(($SLURM_ARRAY_TASK_ID - 1))
SAMPLE=${SAMPLES[$SAMPLE_INDEX]}

# Define a unique output directory for each sample
OUT_DIR="/blue/timgarrett/hkates/Campbell-Thompson/scRNA/NS3809/results/cellranger/${SAMPLE}_${SLURM_ARRAY_TASK_ID}"

# Define the path to the reference transcriptome
REF_DIR="/blue/timgarrett/hkates/Campbell-Thompson/scRNA/NS3809/scripts/refdata-gex-GRCh38-2024-A"

# Create a directory for symbolic links for the current sample
LINK_DIR="/blue/timgarrett/hkates/Campbell-Thompson/scRNA/NS3809/symlinks/${SAMPLE}_${SLURM_ARRAY_TASK_ID}"
mkdir -p ${LINK_DIR}

# Remove any existing symbolic links in the directory
rm -f ${LINK_DIR}/*

# Create symbolic links for the original FASTQ files in the link directory
for lane in 1 2 3 4 5 6; do
    ln -s ${BASE_DIR1}/${SAMPLE}*L00${lane}/*_R1_001.fastq.gz ${LINK_DIR}
    ln -s ${BASE_DIR1}/${SAMPLE}*L00${lane}/*_R2_001.fastq.gz ${LINK_DIR}
    ln -s ${BASE_DIR2}/${SAMPLE}*L00${lane}/*_R1_001.fastq.gz ${LINK_DIR}
    ln -s ${BASE_DIR2}/${SAMPLE}*L00${lane}/*_R2_001.fastq.gz ${LINK_DIR}
done

# Extract the actual sample name from one of the FASTQ files
ACTUAL_SAMPLE=$(ls ${LINK_DIR}/*_R1_001.fastq.gz | head -n 1 | xargs -n 1 basename | cut -d '_' -f 1)

# Remove existing output directory for the sample to avoid conflicts
if [ -d ${OUT_DIR} ]; then
    rm -rf ${OUT_DIR}
fi

# Run Cell Ranger
cellranger count --id=${SAMPLE}_${SLURM_ARRAY_TASK_ID} \
                 --fastqs=${LINK_DIR} \
                 --sample=${ACTUAL_SAMPLE} \
                 --transcriptome=${REF_DIR} \
                 --localcores=8 \
                 --localmem=64 \
                 --create-bam=true \
                 --output-dir=${OUT_DIR}

# Clean up the symbolic link directory if desired (optional)
# rm -rf ${LINK_DIR}

echo -e "\nInfo: Cell Ranger job completed on $(date) on $(hostname).\n"
##End


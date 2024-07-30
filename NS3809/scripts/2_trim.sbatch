#!/bin/bash
#SBATCH --job-name=trim_galore.%j
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16GB
#SBATCH --qos=timgarrett-b
#SBATCH --account=timgarrett
#SBATCH -t 16:00:00
#SBATCH --output=/blue/timgarrett/hkates/Campbell-Thompson/scRNA/NS3809/logs/trim_galore.%j.out
#SBATCH --error=/blue/timgarrett/hkates/Campbell-Thompson/scRNA/NS3809/logs/trim_galore.%j.err
#SBATCH --array=1-16

echo "This is task $SLURM_ARRAY_TASK_ID"

echo -e "\nInfo: Starting a job on $(date) on $(hostname) in $(pwd).\n"

###################################################################
# This file is to perform trimming and filtering using Trim Galore!
###################################################################

module load trim_galore/0.6.10

# Define the base directory for concatenated files
BASE_DIR="/orange/timgarrett/Campbell_Thompson/NS3809/Concatenated"

# Create an array of unique sample identifiers dynamically
SAMPLES=($(ls ${BASE_DIR} | cut -d "_" -f 1,2,3 | sort | uniq))

# Get the sample for the current task
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID-1]}

# Define the output directory
OUT_DIR="/blue/timgarrett/hkates/Campbell-Thompson/scRNA/NS3809/results/trimgalore"

# Define the path to the FastQC summary file
FASTQC_SUMMARY="/blue/timgarrett/hkates/Campbell-Thompson/scRNA/NS3809/results/fastqc/fastqc_summary.csv"

# Function to determine trim parameters based on FastQC results
get_trim_params() {
  local sample=$1
  local summary_file=$2

  # Extract the necessary fields for the sample
  local result=$(grep "^${sample}," $summary_file)

  # Initialize trim parameters
  local params=""

  # Check for issues in the FastQC results and set parameters accordingly
  if [[ $result == *",FAIL,"* ]] || [[ $result == *",WARN,"* ]]; then
    params="--length 30"
    if [[ $result == *",Adapter Content,FAIL,"* ]] || [[ $result == *",Adapter Content,WARN,"* ]]; then
      params="$params --adapter AGATCGGAAGAGC"
    fi
    if [[ $result == *",Per base sequence quality,FAIL,"* ]] || [[ $result == *",Per base sequence quality,WARN,"* ]]; then
      params="$params --quality 20"
    fi
  fi

  # Return the determined parameters
  echo $params
}

# Get the parameters for the current sample
TRIM_PARAMS=$(get_trim_params $SAMPLE $FASTQC_SUMMARY)

# Run Trim Galore with the dynamically set parameters
trim_galore --cores 4 $TRIM_PARAMS --paired ${BASE_DIR}/${SAMPLE}_combined_R1.fastq.gz ${BASE_DIR}/${SAMPLE}_combined_R2.fastq.gz -o ${OUT_DIR}

echo -e "\nInfo: Trim Galore! job completed on $(date) on $(hostname).\n"
##End


#!/bin/bash
#SBATCH --job-name=soupx_and_seurat
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=125GB
#SBATCH --qos=timgarrett
#SBATCH --account=timgarrett
#SBATCH -t 05:00:00
#SBATCH --output=/blue/timgarrett/hkates/Campbell-Thompson/scRNA/NS3809/logs/soupx_3a.%j.out
#SBATCH --error=/blue/timgarrett/hkates/Campbell-Thompson/scRNA/NS3809/logs/soupx_3a.%j.err

module load R/4.4
Rscript 3a_soupx_create_seurat.R

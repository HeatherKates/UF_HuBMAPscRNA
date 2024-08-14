#!/bin/bash
#SBATCH --job-name=seurat_read_and_combine
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64GB
#SBATCH --qos=timgarrett-b
#SBATCH --account=timgarrett
#SBATCH -t 02:00:00
#SBATCH --output=/blue/timgarrett/hkates/Campbell-Thompson/scRNA/NS3809/logs/seurat_3a.out
#SBATCH --error=/blue/timgarrett/hkates/Campbell-Thompson/scRNA/NS3809/logs/seurat_3a.%j.err

module load R/4.2
Rscript 3a*R

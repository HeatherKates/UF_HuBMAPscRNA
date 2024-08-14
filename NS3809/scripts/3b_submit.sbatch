#!/bin/bash
#SBATCH --job-name=seurat_integrate
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=100GB
#SBATCH --qos=timgarrett-b
#SBATCH --account=timgarrett
#SBATCH -t 02:00:00
#SBATCH --output=/blue/timgarrett/hkates/Campbell-Thompson/scRNA/NS3809/logs/seurat_3b.out
#SBATCH --error=/blue/timgarrett/hkates/Campbell-Thompson/scRNA/NS3809/logs/seurat_3b.err

module load R/4.2
Rscript 3b*R

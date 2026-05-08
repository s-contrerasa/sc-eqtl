#!/bin/bash
#SBATCH --job-name=scrna_qc
#SBATCH --account=gao824
#SBATCH --partition=cpu
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --output=01_scrna_quality_control/logs/qc_%j.out
#SBATCH --error=01_scrna_quality_control/logs/qc_%j.err

module load r/4.4.1
Rscript 01_scrna_quality_control/scripts/scrna_qc.R

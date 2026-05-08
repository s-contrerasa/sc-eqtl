#!/bin/bash
#SBATCH --job-name=filter_celltypes
#SBATCH --account=gao824
#SBATCH --partition=cpu
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --output=02_expression_normalization/logs/filter_%j.out
#SBATCH --error=02_expression_normalization/logs/filter_%j.err

module load r/4.4.1
Rscript 02_expression_normalization/scripts/filter_celltypes.R

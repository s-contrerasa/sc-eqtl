#!/bin/bash
#SBATCH --job-name=sctransform
#SBATCH --account=gao824
#SBATCH --partition=cpu
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --output=02_expression_normalization/logs/normalize_%j.out
#SBATCH --error=02_expression_normalization/logs/normalize_%j.err

module load r/4.4.1
Rscript 02_expression_normalization/scripts/sctransform_normalize.R

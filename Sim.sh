#!/bin/bash
#SBATCH --job-name=KernelDA
#SBATCH -o Results-XX/Rep_%a.Rout
#SBATCH --array=1-1200
#SBATCH --mail-user=xxx@xyz
#SBATCH --mail-type=END
#SBATCH --account=abc
#SBATCH --qos=abc
#SBATCH --mem=2gb
#SBATCH -t 96:00:00

ml R/4.2
R CMD BATCH --vanilla Main-XX.R  Results-XX/Rep_${SLURM_ARRAY_TASK_ID}.Rout

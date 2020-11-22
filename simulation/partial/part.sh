#!/bin/bash
#SBATCH --account=def-tmparker
#SBATCH --mem-per-cpu=250M
#SBATCH --time=00-03:00 # time (DD-HH:MM)
#SBATCH --array=31,47

module load gcc r
Rscript partialID_normal.R $SLURM_ARRAY_TASK_ID


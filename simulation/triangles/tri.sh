#!/bin/bash
#SBATCH --account=def-tmparker
#SBATCH --mem-per-cpu=250M
#SBATCH --time=0-15:00 # time (DD-HH:MM)

module load gcc
module load r

Rscript pointID_tri.R


#!/bin/bash
#SBATCH --job-name=kzmi_grid
#SBATCH --partition=sun.q
#SBATCH --account=sun 
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --error=%A.err
#SBATCH --time=0-23:59:59 ## time format is DD-HH:MM:SS
#SBATCH --output=%A.out


srun home/kazuumi/lus/B0/a.out

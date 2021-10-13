#!/bin/sh

#SBATCH --mail-type=FAIL
#SBATCH --mail-user=beckandy@umich.edu
#SBATCH --mem=1500MB
#SBATCH --ntasks=1
#SBATCH --time 12:00:00
#SBATCH --job-name=2pos_cts
#SBATCH --partition=nomosix
#SBATCH --array=1-9
#SBATCH --requeue
#SBATCH -e slurm/step7-%J.err
#SBATCH -o slurm/step7-%J.out


srun Rscript /net/snowwhite/home/beckandy/research/BRIDGES_redo/src/step8_2pos_counts.R ${SLURM_ARRAY_TASK_ID}

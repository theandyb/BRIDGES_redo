#!/bin/sh

#SBATCH --mail-type=FAIL
#SBATCH --mail-user=beckandy@umich.edu
#SBATCH --mem=1500MB
#SBATCH --ntasks=1
#SBATCH --time 12:00:00
#SBATCH --job-name=2pos_cts
#SBATCH --array=1-9
#SBATCH --requeue
#SBATCH -e slurm/step8-%J.err
#SBATCH -o slurm/step8-%J.out


Rscript /net/snowwhite/home/beckandy/research/BRIDGES_redo/src/2_pos_model_data.R ${SLURM_ARRAY_TASK_ID}

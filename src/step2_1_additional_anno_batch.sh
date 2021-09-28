#!/bin/sh
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=beckandy@umich.edu
#SBATCH --ntasks=1
#SBATCH --mem=10GB
#SBATCH --time=01:30:00
#SBATCH --job-name=ALL
#SBATCH --partition=nomosix
#SBATCH --array=1-22
#SBATCH --requeue
#SBATCH -e slurm/step2_1_addl-%J.err
#SBATCH -o slurm/step2_1_addl-%J.out

echo ${SLURM_ARRAY_TASK_ID}

srun Rscript /net/snowwhite/home/beckandy/research/BRIDGES_redo/src/step2_additional_anno.R /net/snowwhite/home/beckandy/research/BRIDGES_redo/output/singletons/chr${SLURM_ARRAY_TASK_ID}_filtered_annotated

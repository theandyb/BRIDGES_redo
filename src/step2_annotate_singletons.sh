#!/bin/sh

#SBATCH --mail-type=FAIL
#SBATCH --mail-user=beckandy@umich.edu
#SBATCH --mem=4GB
#SBATCH --ntasks=1
#SBATCH --time 02:00:00
#SBATCH --job-name=annotate
#SBATCH --partition=nomosix
#SBATCH --array=1-22
#SBATCH --requeue
#SBATCH -e slurm/step2-%J.err
#SBATCH -o slurm/step2-%J.out

python /net/snowwhite/home/beckandy/research/BRIDGES_redo/src/step2_append_motif.py -c ${SLURM_ARRAY_TASK_ID} -s /net/snowwhite/home/beckandy/research/BRIDGES_redo/output/singletons/chr${SLURM_ARRAY_TASK_ID}_filtered.singletons -o /net/snowwhite/home/beckandy/research/BRIDGES_redo/output/singletons/chr${SLURM_ARRAY_TASK_ID}_filtered_annotated.csv

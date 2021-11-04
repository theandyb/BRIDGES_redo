#!/bin/sh

#SBATCH --mail-type=FAIL
#SBATCH --mail-user=beckandy@umich.edu
#SBATCH --mem=1500MB
#SBATCH --ntasks=1
#SBATCH --time 01:00:00
#SBATCH --job-name=BRIDGES
#SBATCH --partition=nomosix
#SBATCH --array=1-4180
#SBATCH --requeue
#SBATCH -e slurm/step6-%J.err
#SBATCH -o slurm/step6-%J.out


srun python /net/snowwhite/home/beckandy/research/BRIDGES_redo/src/gw_2_count_6cats.py \
-o /net/snowwhite/home/beckandy/research/BRIDGES_redo/output/gw_2_count/6cats \
-j ${SLURM_ARRAY_TASK_ID}

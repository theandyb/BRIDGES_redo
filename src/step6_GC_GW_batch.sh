#!/bin/sh

#SBATCH --mail-type=FAIL
#SBATCH --mail-user=beckandy@umich.edu
#SBATCH --mem=1500MB
#SBATCH --ntasks=1
#SBATCH --time 12:00:00
#SBATCH --job-name=GC_GW
#SBATCH --partition=nomosix
#SBATCH --array=1-22
#SBATCH --requeue
#SBATCH -e slurm/step6_GC-%J.err
#SBATCH -o slurm/step6_GC-%J.out


srun python /net/snowwhite/home/beckandy/research/BRIDGES_redo/src/step6_GC_GW_2_count.py -o /net/snowwhite/home/beckandy/research/BRIDGES_redo/output/gw_2_count/cpg -c ${SLURM_ARRAY_TASK_ID}

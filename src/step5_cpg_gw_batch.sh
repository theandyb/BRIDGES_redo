#!/bin/sh

#SBATCH --mail-type=FAIL
#SBATCH --mail-user=beckandy@umich.edu
#SBATCH --mem=1500MB
#SBATCH --ntasks=1
#SBATCH --time 01:10:00
#SBATCH --job-name=CG_C_1
#SBATCH --partition=nomosix
#SBATCH --array=1-22
#SBATCH --requeue
#SBATCH -e slurm/step5_cpg-%J.err
#SBATCH -o slurm/step5_cpg-%J.out

python /net/snowwhite/home/beckandy/research/BRIDGES_redo/src/step5_cpg_gw_1_count.py -o /net/snowwhite/home/beckandy/research/BRIDGES_redo/output/gw_1_count/cpg -c ${SLURM_ARRAY_TASK_ID}

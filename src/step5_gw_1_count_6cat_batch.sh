#!/bin/sh
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=beckandy@umich.edu
#SBATCH --ntasks=1
#SBATCH --mem=2GB
#SBATCH --time=06:00:00
#SBATCH --job-name=gw_1_6
#SBATCH --partition=nomosix
#SBATCH --array=1-22
#SBATCH --requeue
#SBATCH -e slurm/step5_1-%J.err
#SBATCH -o slurm/step5_1-%J.out

python /net/snowwhite/home/beckandy/research/BRIDGES_redo/src/gw_1_count_6cats.py ${SLURM_ARRAY_TASK_ID}

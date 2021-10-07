#!/bin/sh
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=beckandy@umich.edu
#SBATCH --ntasks=1
#SBATCH --mem=12GB
#SBATCH --time=01:30:00
#SBATCH --job-name=at_new
#SBATCH --partition=nomosix
#SBATCH --array=1-22
#SBATCH --requeue
#SBATCH -e slurm/step3_1_at-%J.err
#SBATCH -o slurm/step3_1_at-%J.out

echo ${SLURM_ARRAY_TASK_ID}

srun python /net/snowwhite/home/beckandy/research/BRIDGES_redo/src/step3_1_reverse_comp_control.py -i /net/snowwhite/home/beckandy/research/BRIDGES_redo/output/control2/at/chr${SLURM_ARRAY_TASK_ID} -n A

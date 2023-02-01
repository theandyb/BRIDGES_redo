#!/bin/sh
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=beckandy@umich.edu
#SBATCH --ntasks=1
#SBATCH --mem=1GB
#SBATCH --time=03:00:00
#SBATCH --job-name=all_gc
#SBATCH --array=1-22
#SBATCH --requeue
#SBATCH -e slurm/all_gc-%A_%a.err
#SBATCH -o slurm/all_gc-%A_%a.out

echo "GC ${SLURM_ARRAY_TASK_ID}"
python /net/snowwhite/home/beckandy/research/BRIDGES_redo/src/control_sample_all_gc.py \
-s /net/snowwhite/home/beckandy/research/BRIDGES_redo/output/singletons/chr${SLURM_ARRAY_TASK_ID}_annotated.csv \
-f /net/snowwhite/home/beckandy/research/BRIDGES_redo/reference_data/masked_ref.fa \
-o /net/snowwhite/home/beckandy/research/BRIDGES_redo/output/controls/chr${SLURM_ARRAY_TASK_ID}_gc_all.csv \
-n 5 ${SLURM_ARRAY_TASK_ID}

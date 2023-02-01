#!/bin/sh
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=beckandy@umich.edu
#SBATCH --ntasks=1
#SBATCH --mem=1GB
#SBATCH --time=03:00:00
#SBATCH --job-name=sample_GC
#SBATCH --array=1-22
#SBATCH --requeue
#SBATCH -e slurm/step3_GC-%A_%a.err
#SBATCH -o slurm/step3_GC-%A_%a.out

echo "GC ${SLURM_ARRAY_TASK_ID} $POP"
python /net/snowwhite/home/beckandy/research/BRIDGES_redo/src/control_sample_gc.py \
-s /net/snowwhite/home/beckandy/research/BRIDGES_redo/output/singletons/chr${SLURM_ARRAY_TASK_ID}_annotated.csv \
-f /net/snowwhite/home/beckandy/research/BRIDGES_redo/reference_data/masked_ref.fa \
-o /net/snowwhite/home/beckandy/research/BRIDGES_redo/output/controls/chr${SLURM_ARRAY_TASK_ID}_gc.csv \
-n 5 ${SLURM_ARRAY_TASK_ID}

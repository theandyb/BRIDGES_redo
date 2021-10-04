#!/bin/sh
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=beckandy@umich.edu
#SBATCH --ntasks=1
#SBATCH --mem=1GB
#SBATCH --time=01:30:00
#SBATCH --job-name=sample_GC
#SBATCH --partition=nomosix
#SBATCH --array=1-22
#SBATCH --requeue
#SBATCH -e slurm/step3_iCpG-%J.err
#SBATCH -o slurm/step3_iCpG-%J.out

echo "GC ${SLURM_ARRAY_TASK_ID}"
python /net/snowwhite/home/beckandy/research/BRIDGES_redo/src/step3_sample_gc_ignore_cpg.py -s /net/snowwhite/home/beckandy/research/BRIDGES_redo/output/singletons/chr${SLURM_ARRAY_TASK_ID}_filtered_annotated_addl.csv -f /net/snowwhite/home/beckandy/research/BRIDGES_redo/reference_data/human_g1k_v37.fasta -o /net/snowwhite/home/beckandy/research/BRIDGES_redo/output/control2/gc/ignore_cpg/chr${SLURM_ARRAY_TASK_ID}.csv -n 5 ${SLURM_ARRAY_TASK_ID}


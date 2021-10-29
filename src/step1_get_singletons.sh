#!/bin/sh

#SBATCH --mail-type=FAIL
#SBATCH --mail-user=beckandy@umich.edu
#SBATCH --mem=700MB
#SBATCH --ntasks=1
#SBATCH --time 04:00:00
#SBATCH --job-name=step1
#SBATCH --partition=nomosix
#SBATCH --array=1-22
#SBATCH --requeue
#SBATCH -e slurm/step1_all-%J.err
#SBATCH -o slurm/step1_all-%J.out

vcftools --gzvcf /net/bipolar/lockeae/final_freeze/snps/vcfs/chr${SLURM_ARRAY_TASK_ID}/chr${SLURM_ARRAY_TASK_ID}.filtered.modified.vcf.gz \
--keep /net/snowwhite/home/beckandy/research/BRIDGES_redo/reference_data/keep_ids.csv \
--remove-indels --remove-filtered-all --singletons \
--out /net/snowwhite/home/beckandy/research/BRIDGES_redo/output/singletons/chr${SLURM_ARRAY_TASK_ID}.txt
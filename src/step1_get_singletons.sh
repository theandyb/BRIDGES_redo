#!/bin/sh

#SBATCH --mail-type=FAIL
#SBATCH --mail-user=beckandy@umich.edu
#SBATCH --mem=700MB
#SBATCH --ntasks=1
#SBATCH --time 12:00:00
#SBATCH --job-name=step1
#SBATCH --partition=nomosix
#SBATCH --array=1-22
#SBATCH --requeue
#SBATCH -e slurm/step1_all-%J.err
#SBATCH -o slurm/step1_all-%J.out

bcftools reheader -h /net/snowwhite/home/beckandy/research/BRIDGES_redo/reference_data/header_${SLURM_ARRAY_TASK_ID}.txt /net/bipolar/lockeae/final_freeze/snps/vcfs/chr${SLURM_ARRAY_TASK_ID}/chr${SLURM_ARRAY_TASK_ID}.filtered.modified.vcf.gz |\
bcftools plugin fill-AN-AC -O v |\
bcftools view --types 'snps' -f 'FILTER="PASS"' -C 1 -c 1 -O v |\
bcftools view -S /net/snowwhite/home/beckandy/research/BRIDGES_redo/reference_data/keep_ids.csv |\
vcftools --vcf - --singletons -c > /net/snowwhite/home/beckandy/research/BRIDGES_redo/output/singletons/chr${SLURM_ARRAY_TASK_ID}.singletons

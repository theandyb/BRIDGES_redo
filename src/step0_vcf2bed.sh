#!/bin/sh

#SBATCH --mail-type=FAIL
#SBATCH --mail-user=beckandy@umich.edu
#SBATCH --mem=1GB
#SBATCH --ntasks=1
#SBATCH --time 05:00:00
#SBATCH --job-name=step0
#SBATCH --array=1-22
#SBATCH --requeue
#SBATCH -e slurm/step0-%A_%a.err
#SBATCH -o slurm/step0-%A_%a.out

bcftools reheader -h /net/snowwhite/home/beckandy/research/BRIDGES_redo/reference_data/header_${SLURM_ARRAY_TASK_ID}.txt \
  /net/bipolar/lockeae/final_freeze/snps/vcfs/chr${SLURM_ARRAY_TASK_ID}/chr${SLURM_ARRAY_TASK_ID}.filtered.modified.vcf.gz |\
  bcftools view -G --types 'snps' -f 'FILTER="PASS"' -Ov |\
  bcftools query -f '%CHROM\t%POS\t%POS\n' | awk '{print($1"\t"$2 - 1"\t"$3)}' >\
  /net/snowwhite/home/beckandy/research/BRIDGES_redo/reference_data/vcfbed/chr${SLURM_ARRAY_TASK_ID}.bed

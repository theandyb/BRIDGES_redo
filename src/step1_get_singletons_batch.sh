#!/bin/sh

#SBATCH --mail-type=FAIL
#SBATCH --mail-user=beckandy@umich.edu
#SBATCH --mem=700MB
#SBATCH --ntasks=1
#SBATCH --time 06:00:00
#SBATCH --job-name=singletons
#SBATCH --partition=nomosix
#SBATCH --array=1-22
#SBATCH --requeue
#SBATCH -e slurm/step1_%J.out
#SBATCH -o slurm/step1_%J.err
#SBATCH --open-mode=append

# PURPOSE: Run the '--singletons' command on the BRIDGES vcf data in order to generate files with 
#          chromosome, position, ref, alt, and ID


vcftools --gzvcf /net/bipolar/lockeae/final_freeze/snps/vcfs/chr${SLURM_ARRAY_TASK_ID}/chr${SLURM_ARRAY_TASK_ID}.filtered.modified.vcf.gz --remove-indels --remove-filtered-all --singletons --out /net/snowwhite/home/beckandy/research/BRIDGES_redo/output/singletons/chr${SLURM_ARRAY_TASK_ID}

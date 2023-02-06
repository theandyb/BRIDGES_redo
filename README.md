# Introduction

The point of this repository is to re-do my analyses of the germline mutation rate using singletons as a proxy. I want to ensure that my results using the BRIDGES dataset are both correct and reproducible. To that end I am making a few changes:

1. I'd like to stream-line the pipeline as much as possible. While a tool like snakemake would ensure this, for now I'm going to maintain a connected set of "steps" encoded in scripts to generate my results. (I really should learn snakemake though).
2. If I get run over by a bus, someone else should be able to read this document, and based on that (and access to the data) they should be able to reproduce my results.

# Step 0: Directory Structure

```{bash}

mkdir output
mkdir output/control
mkdir output/control/at
mkdir output/control/gc
mkdir output/control2
mkdir output/control2/at
mkdir output/control2/gc
mkdir output/gw_1_count
mkdir output/gw_1_count/cpg
mkdir output/gw_2_count
mkdir output/gw_2_count/at
mkdir output/gw_2_count/cpg
mkdir output/single_pos_df
mkdir output/singletons
mkdir output/all_count_2_pos

```

# Step 0.5: Setting up an environment

```{bash}
# assuming you've install mambaforge 
# https://mamba.readthedocs.io/en/latest/installation.html

mamba create -n bridges pyfaidx r-base r-tidyverse pandas r-doparallel r-writexl regex biopython -c conda-forge -c bioconda

```

# Step 0: Mask reference genome

For our analyses we will want to mask sites in the reference genome which are variant sites in the BRIDGES data set. In order to do this, we will:

1. Convert the VCF to a list of variant sites in BED format using `src/step0_vcf2bed.sh`
2. Mask the reference genome using the `maskfasta` option of `bedtools`

After running `step0_vcf2bed.sh`, we can then combine the bed files into one using the following commands:

```
for i in `seq 1 22`; do
cat chr${i}.bed >> bridges_sites.bed
done
```

and then mask the reference genome:

```
bedtools maskfasta -fi /net/snowwhite/home/beckandy/research/BRIDGES_redo/reference_data/human_g1k_v37.fasta -bed /net/snowwhite/home/beckandy/research/BRIDGES_redo/reference_data/vcfbed/bridges_sites.bed -fo /net/snowwhite/home/beckandy/research/BRIDGES_redo/reference_data/masked_ref.fa
```

which we can index with samtools:

```
samtools faidx /net/snowwhite/home/beckandy/research/BRIDGES_redo/reference_data/masked_ref.fa
```

# Step 1: Generate singleton files

For this step we simply rely on the vcftools `--singletons` option. This is done via a batch script (`src/step1_get_singletons_batch.sh`).

# Step 2: Annotate singletons

The script which appends each singleton with its motif is `src/append_motif.py`, and the batch script `step2_annotate_batch.sh` will submit slurm jobs for each chromosome. The output files will be headerless csvs with the following columns:

1. Chromosome
2. Position
3. Original Motif
4. Simple subtype (REF>ALT)
5. ALT
6. Sample ID
7. REF
8. Full Motif (motif and its reverse complement)
9. Condensed sub-type

# Step 3: Sample control distribution
The scripts to sample 5 controls per singleton are `src/control_sample_at.py` and `src/control_sample_gc.py`. Scripts to submit batch jobs to slurm are available in files with the prefix `step3_sample` (also in the src directory).

NOTE: if using a reference genome other than hg38, you might need to change the variable ref_prefix in the main function definition of the two scripts to match the chromosome names in the fasta file (i.e. in hg37, chromosomes are named only by their number, where as in hg38 each chromosome is prefixed with "chr", e.g. chr1, chr2, ...)

# Step 4: Per-chromosome files -> per subtype files

In the above steps, we've generated files per chromosome. In this step we'll compile singletons and controls across all chromosomes.

Within `output/singletons`, run the following commands (note: this will take a few minutes to run):

```{bash}
# Generate per-subtype singleton files

awk -F, '{if($9 == "AT_CG")print(substr($8,1,21))}' chr*_annotated.csv > AT_CG.txt
awk -F, '{if($9 == "AT_GC")print(substr($8,1,21))}' chr*_annotated.csv > AT_GC.txt
awk -F, '{if($9 == "AT_TA")print(substr($8,1,21))}' chr*_annotated.csv > AT_TA.txt
awk -F, '{if($9 == "GC_AT")print(substr($8,1,21))}' chr*_annotated.csv > GC_AT.txt
awk -F, '{if($9 == "GC_TA")print(substr($8,1,21))}' chr*_annotated.csv > GC_TA.txt
awk -F, '{if($9 == "GC_CG")print(substr($8,1,21))}' chr*_annotated.csv > GC_CG.txt
awk -F, '{if($9 == "cpg_GC_AT")print(substr($8,1,21))}' chr*_annotated.csv > cpg_GC_AT.txt
awk -F, '{if($9 == "cpg_GC_TA")print(substr($8,1,21))}' chr*_annotated.csv > cpg_GC_TA.txt
awk -F, '{if($9 == "cpg_GC_CG")print(substr($8,1,21))}' chr*_annotated.csv> cpg_GC_CG.txt
```

Within `output/controls`, run the following commands to yield per-subtype files:

```{bash}
awk -F, '{if($4 == "AT_CG")print(substr($8,1,21))}' chr*_at.csv > AT_CG.txt
awk -F, '{if($4 == "AT_GC")print(substr($8,1,21))}' chr*_at.csv > AT_GC.txt
awk -F, '{if($4 == "AT_TA")print(substr($8,1,21))}' chr*_at.csv > AT_TA.txt
awk -F, '{if($4 == "GC_AT")print(substr($8,1,21))}' chr*_gc.csv > GC_AT.txt
awk -F, '{if($4 == "GC_TA")print(substr($8,1,21))}' chr*_gc.csv > GC_TA.txt
awk -F, '{if($4 == "GC_CG")print(substr($8,1,21))}' chr*_gc.csv > GC_CG.txt
awk -F, '{if($4 == "cpg_GC_AT")print(substr($8,1,21))}' chr*_gc.csv > cpg_GC_AT.txt
awk -F, '{if($4 == "cpg_GC_TA")print(substr($8,1,21))}' chr*_gc.csv > cpg_GC_TA.txt
awk -F, '{if($4 == "cpg_GC_CG")print(substr($8,1,21))}' chr*_gc.csv> cpg_GC_CG.txt
awk -F, '{if($4 == "GC_AT" || $4 == "cpg_GC_AT")print(substr($8,1,21))}' chr*_gc_all.csv > all_GC_AT.txt
awk -F, '{if($4 == "GC_TA" || $4 == "cpg_GC_AT")print(substr($8,1,21))}' chr*_gc_all.csv > all_GC_TA.txt
awk -F, '{if($4 == "GC_CG" || $4 == "cpg_GC_AT")print(substr($8,1,21))}' chr*_gc_all.csv > all_GC_CG.txt
```

If wanting to perform analyses using nearest and furthest controls:

```
awk -F, '{if($4 == "AT_CG")print(substr($8,1,21))}' chr*_at.csv.min > AT_CG_min.txt
awk -F, '{if($4 == "AT_GC")print(substr($8,1,21))}' chr*_at.csv.min > AT_GC_min.txt
awk -F, '{if($4 == "AT_TA")print(substr($8,1,21))}' chr*_at.csv.min > AT_TA_min.txt
awk -F, '{if($4 == "GC_AT")print(substr($8,1,21))}' chr*_gc.csv.min > GC_AT_min.txt
awk -F, '{if($4 == "GC_TA")print(substr($8,1,21))}' chr*_gc.csv.min > GC_TA_min.txt
awk -F, '{if($4 == "GC_CG")print(substr($8,1,21))}' chr*_gc.csv.min > GC_CG_min.txt
awk -F, '{if($4 == "cpg_GC_AT")print(substr($8,1,21))}' chr*_gc.csv.min > cpg_GC_AT_min.txt
awk -F, '{if($4 == "cpg_GC_TA")print(substr($8,1,21))}' chr*_gc.csv.min > cpg_GC_TA_min.txt
awk -F, '{if($4 == "cpg_GC_CG")print(substr($8,1,21))}' chr*_gc.csv.min > cpg_GC_CG_min.txt

awk -F, '{if($4 == "AT_CG")print(substr($8,1,21))}' chr*_at.csv.max > AT_CG_max.txt
awk -F, '{if($4 == "AT_GC")print(substr($8,1,21))}' chr*_at.csv.max > AT_GC_max.txt
awk -F, '{if($4 == "AT_TA")print(substr($8,1,21))}' chr*_at.csv.max > AT_TA_max.txt
awk -F, '{if($4 == "GC_AT")print(substr($8,1,21))}' chr*_gc.csv.max > GC_AT_max.txt
awk -F, '{if($4 == "GC_TA")print(substr($8,1,21))}' chr*_gc.csv.max > GC_TA_max.txt
awk -F, '{if($4 == "GC_CG")print(substr($8,1,21))}' chr*_gc.csv.max > GC_CG_max.txt
awk -F, '{if($4 == "cpg_GC_AT")print(substr($8,1,21))}' chr*_gc.csv.max > cpg_GC_AT_max.txt
awk -F, '{if($4 == "cpg_GC_TA")print(substr($8,1,21))}' chr*_gc.csv.max > cpg_GC_TA_max.txt
awk -F, '{if($4 == "cpg_GC_CG")print(substr($8,1,21))}' chr*_gc.csv.max > cpg_GC_CG_max.txt
```

# Step 5: Genome-wide Background Rates - Single Position Models

The scripts to generate counts based on the reference genome are `gw_1_count3cats.py` and `gw_1_count_6cats.py`. Batch scripts to submit jobs to slurm are `step5_gw_1_count_3cat_batch.sh` and `step5_gw_1_count_6cat_batch.sh`.

After this, run the script `step5_combine_chromosomes.R` to combine the per-chromosome per-reference files into per-reference files.


# Step 6: Genome-wide Background Rates - Two Postion Models

The script to generate counts based on the reference genome is `gw_2_count_6cats.py`, with batch script to submit jobs to slurm `step6_gw_2_6cats.sh`.
# Step 7: Generate Single Position Model Tables

For each sub-type, we aggregate all the data needed to fit the models at each flanking position. Each position will have a table, and each nucleotide will have a value for singleton counts, control distribution counts, and genome-wide counts. We will also compute expectations for the singletons based on rates from the control and genome-wide counts, and use these expectations to compute chi-square residuals (which we can then sum across the nucleotides to obtain the chi-square goodness of fit statistic for that position).

This is done by the script `step7_single_position_models.R`

# Step 8: Generate Two Position Count Tables

For each sub-type and pair of positions, we obtain the 16 counts for the singletons, controls, and genome-wide background rate. This is done by the script `step8_2pos_batch.sh`

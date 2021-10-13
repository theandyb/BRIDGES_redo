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

```

# Step 0.5: Setting up an environment

```{bash}
# assuming you've install mambaforge 
# https://mamba.readthedocs.io/en/latest/installation.html

mamba create -n bridges pyfaidx r-base r-tidyverse pandas r-doparallel r-writexl regex biopython -c conda-forge -c bioconda

```

# Step 1: Generate singleton files

For this step we simply rely on the vcftools `--singletons` option. This is done via a batch script (`src/step1_get_singletons_batch.sh`). We also need to filter out the subjects that were excluded in Carlson, et al (2017); this I'll do via an R script (`src/step1_filter_subjects.R`) 

# Step 2: Annotate singletons

For this step, we will loop over all the singletons and pull the 21-mer motif centered at each position. We do this in two stages; first, we run `step2_append_motif.py` to annotate each position with the 21-mer motif and the simple subtype. We then run `step2_additional_anno.R` to take reverse-complement (when necessary) and produce the "full" sub-type (i.e. A>C is now AT>CG, C>T is now GC>AT, etc).

The batch scripts for this step are `step2_annotate_singletons.sh` and `step2_1_additional_anno_batch.sh`

# Step 3: Sample control distribution

I do this two different ways; the first I match only by the REF allele; these results are in the `output/count/` directory. The second approach, which I **no longer** think is more correct, is to match As and Ts to either As or Ts; CpGs are matched to CG, but we allow for either the "C position" or the "G position" to be the center of the motif (in the latter, the reverse-complement has a CpG at the center); non-CpGs are matched to either \[ACT\]C\[ACT\] or \[AGT\]G\[AGT\] (that is, non-CpGs). The latter samples are stored in `output/control2`

The batch scripts for approach 1 are `step3_sample_gc_batch.sh` and `step3_sample_at_batch.sh`

The batch scripts for approach 2 are `step3_sample_gc_batch_2.sh` and `step3_sample_at_batch_2.sh`

I also do a version of the sampling where I match all GC_ mutations to any C or G in the reference genome (ignoring CpG status). This should allow me to demonstrate the known CpG effect at the +1.

We then need to do similar processing with these motifs, since we'll want to take reverse-complements when the center is either a G or T. The batch script for this procedure is `step3_1_reverse_comp_control.sh`.

# Step 4: Per-chromosome files -> per subtype files

This is done with a set of command-line operations you can find in the file `step4_subtype_files.sh`

# Step 5: Get genome-wide rates for single position models

For this analysis, we employ a sliding-window approach to generate the 20 count tables used for each sub-type. Iterating over each position in the genome:

1. Is the nucleotide in the set \{A,C,G,T\}? If not, move to next position
2. Is the nucleotide an A or C? If so, set switch = 1, if not set switch = -1
3. For all relative positions in -10:-1, 1:10
    * Get x, the nucleotide at position rp from current nucleotide
    * If x is not in the set \{A,C,G,T\} move to next rp
    * If switch = -1, set x = rc(x), where rc is the reverse-complement map
    * Table index ix = rp * direction (i.e. if we're taking reverse-complement, -10 is now +10, etc)
    * Lookup count for nucleotide x in table ix and add 1
    
The batch script to generate the counts for each chromosome and rp is `step5_gw_1_count_batch.sh`. An R script to combine results across chromosomes (`step5_add_chromosome.R`) is also in the src directory.

We also would like the C (and G) counts for CpG and non-CpGs separately. This is done by the batch script `step5_cpg_gw_batch.sh` (note: the script it calls does not perform the reverse-complement opperation, and yields 4 tables: CpG-C, CpG-G, C, and G). This is done by the batch script `step5_cpg_gw_batch.sh`. This generates per-chromosome files, which we then combine using the script `step5_cpg_combine_chrom.R`.

# Step 6: Get genome-wide rates for 2 position models

This is similar to step 5, but here we stratify counts based on the nucleotides at all pairs of flanking positions in the +/- 10 bp window. Here we don't worry about taking reverse complements while generating the counts (i.e. for each pair of positions we'll have 4 tables instead of 2).

The batch scripts for generating the counts are `step6_GC_GW_batch.sh` and `step6_AT_GW_batch.sh`. The script `step6_generate_model_tables.R` then creates the tables we need to fit the models at each pair of positions for the 9 subtypes.

# Step 7: Generate Single Position Model Tables

For each sub-type, we aggregate all the data needed to fit the models at each flanking position. Each position will have a table, and each nucleotide will have a value for singleton counts, control distribution counts, and genome-wide counts. We will also compute expectations for the singletons based on rates from the control and genome-wide counts, and use these expectations to compute chi-square residuals (which we can then sum across the nucleotides to obtain the chi-square goodness of fit statistic for that position).

This is done by the script `step7_single_position_models.R`

# Step 8: Generate Two Position Count Tables

For each sub-type and pair of positions, we obtain the 16 counts for the singletons, controls, and genome-wide background rate. This is done by the script `step8_2pos_batch.sh`
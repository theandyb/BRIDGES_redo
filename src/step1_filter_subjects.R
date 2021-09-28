library(tidyverse)
library(foreach)
library(doParallel)

keep_ids <- read_csv("/net/snowwhite/home/beckandy/research/BRIDGES_redo/reference_data/keep_ids.csv", col_types = cols())
keep_list <- keep_ids$ID

singleton_dir <- "/net/snowwhite/home/beckandy/research/BRIDGES_redo/output/singletons/"

process_chrom <- function(i){
	f_name <- paste0(singleton_dir, "chr", i, ".singletons")

	df <- read_tsv(f_name, col_types = cols()) %>%
		filter(INDV %in% keep_list)

	out_name <- paste0(singleton_dir, "chr", i, "_filtered.singletons")
	write_tsv(df, file = out_name)
	return(0)
}

registerDoParallel(22)
foreach(i=1:22) %dopar% {
  process_chrom(i)
}

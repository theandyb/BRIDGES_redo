# Combine genome-wide two-flanking-position counts across chromosomes
# We'll also perform the reverse-complement operation here to merge G -> C and T -> A
library(tidyverse)

#' Load counts for two flanking positions for a focal nucleotide in a chromosome
#' 
#' @param nucleotide The central nucleotide (A, C, G, T)
#' @param p1 First flanking position in +/- 10 bp window (p1 != 0)
#' @param p2 Second flanking position in +/- 10 bp window (p2 > p1)
#' @param chromosome which chromosome to load counts for
#' @return Data frame with counts for each dinucleotide pair at positions p1 and p2
load_positions_chrom <- function(nucleotide, p1, p2, chromosome){
  
  nucs <- c("A", "C", "G", "T")
  
  data_dir <- "/net/snowwhite/home/beckandy/research/BRIDGES_redo/output/gw_2_count/"
  if(str_starts(nucleotide, "[at]")){
    sub_dir <- "at/"
  } else {
    sub_dir <- "cpg/"
  }
  
  f_name <- paste0(data_dir, sub_dir, nucleotide, "_chrom", chromosome, "_p", p1, "_q", p2, ".csv")
  df <- read_csv(f_name, col_types = cols()) %>%
    mutate(p1 = str_sub(Nucs, 1, 1), p2 = str_sub(Nucs, 2, 2)) %>%
    select(p1, p2, n) %>%
    filter(p1 %in% nucs, p2 %in% nucs) %>%
    arrange(p1, p2)
  
  return(df)
  
}

#' Load counts for two flanking positions across all autosomes
#' 
#' @param nucleotide The central nucleotide (A, C, G, T)
#' @param p1 First flanking position in +/- 10 bp window (p1 != 0)
#' @param p2 Second flanking position in +/- 10 bp window (p2 > p1)
#' @return Data frame with counts for each dinucleotide pair at positions p1 and p2
load_positions <- function(nucleotide, p1, p2){
  
  df <- load_positions_chrom(nucleotide, p1, p2, 1)
  
  for(chromosome in c(2:22)){
    df2 <- load_positions_chrom(nucleotide, p1, p2, chromosome)
    
    df <- full_join(df, df2, by = c("p1", "p2")) %>%
      replace_na(list("n.x" = 0, "n.y" = 0)) %>%
      mutate(n = n.x + n.y) %>%
      select(p1, p2, n)
  }
  return(df)
}

#' Single nucleotide reverse-complement map
#' 
#' @param nucleotides Vector of nucleotides which we want to map each individual nucleotide to its reverse complement
#' @return Reverse complement of the nucleotide
simple_rc <- function(nucleotides){
  map_vec <- list("A" = "T",
                  "C" = "G",
                  "G" = "C",
                  "T" = "A")
  rc_nucs <- sapply(nucleotides, function(x){map_vec[[x]]})
  return(rc_nucs)
} 

#' Load all counts at two flanking positions for a specific subtype
#' 
#' @param subtype One of the 9 basic subtypes
#' @param p1 First flanking position in +/- 10 bp window (p1 != 0)
#' @param p2 Second flanking position in +/- 10 bp window (p2 > p1)
#' @return Data frame with counts for each dinucleotide pair at positions p1 and p2
load_subtype_counts <- function(subtype, p1, p2){
  if(str_starts(subtype, "AT")){
    nuc_1 <- "a"
    nuc_2 <- "t"
  } else if(str_starts(subtype, "cpg")){
    nuc_1 <- "cpg_c"
    nuc_2 <- "cpg_g"
  } else {
    nuc_1 <- "non_c"
    nuc_2 <- "non_g"
  }
  
  df1 <- load_positions(nuc_1, p1, p2)
  df2 <- load_positions(nuc_2, -p2, -p1) %>%
    mutate(p1rc = simple_rc(p2),
           p2rc = simple_rc(p1)) %>%
    select(p1rc, p2rc, n) %>%
    rename(p1 = p1rc, p2 = p2rc)
  df <- full_join(df1, df2, by = c("p1", "p2")) %>%
    replace_na(list("n.x" = 0, "n.y" = 0)) %>%
    mutate(n = n.x + n.y) %>%
    select(p1, p2, n)
  
  return(df)
} 

# Let's generate the combined tables for the two position models
# Note that all three AT_NN subtypes use the same background genome-wide rates, etc
for(subtype in c("AT", "GC", "cpg_GC")){
  for(p1 in c(-10:-1,1:9)){
    for(p2 in c((p1+1):10)){
      if(p2 == 0) next
      df <- load_subtype_counts(subtype, p1, p2)
      f_name <- paste0("/net/snowwhite/home/beckandy/research/BRIDGES_redo/output/gw_2_count/", 
                       subtype,
                       "_p", p1,
                       "_q", p2,
                       ".csv")
      write_csv(df, f_name)
    }
  }
}

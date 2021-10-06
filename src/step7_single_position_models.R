# Gather data for the single position models

library(tidyverse)

#' Get the singleton counts for a sub-type at a flanking position
#' 
#' @param subtype One of the nine simple sub types
#' @param rp Relative position which to get singleton counts from
#' @return data frame with counts for each nucleotide
get_singleton_pos <- function(subtype, rp){
  nucs <- c("A", "C", "G", "T")
  
  pos <- rp + 11
  
  singleton_dir <- "/net/snowwhite/home/beckandy/research/BRIDGES_redo/output/singletons/"
  fName <- paste0(singleton_dir, subtype, ".txt")
  
  awk_cmd <- paste0("awk '{count[substr($1,", pos,", 1)]++}END{for(key in count) print(key, \"\\t\", count[key])}' ",
                    fName)
  
  df <- vroom::vroom(pipe(awk_cmd), col_names = c("Nuc", "singletons"), delim = "\t", show_col_types = FALSE) %>%
    filter(Nuc %in% nucs )
  return(df)
}

#' Get nucleotide counts at position relative to control observation
#' 
#' @param subtype One of the nine basic subtypes
#' @param rp Position relative to focal site to get counts for
#' @return Data.frame with control counts stratified by nucleotide at flanking position
get_control_pos <- function(subtype, rp){
  nucs <- c("A", "C", "G", "T")
  
  pos <- rp + 11
  
  if(str_starts(subtype, "AT")){
    control_dir <- "/net/snowwhite/home/beckandy/research/BRIDGES_redo/output/control2/at/"
  } else {
    control_dir <- "/net/snowwhite/home/beckandy/research/BRIDGES_redo/output/control2/cpg/gc/"
  }
  fName <- paste0(control_dir, subtype, ".txt")
  
  awk_cmd <- paste0("awk '{count[substr($1,", pos,", 1)]++}END{for(key in count) print(key, \"\\t\", count[key])}' ",
                    fName)
  
  df <- vroom::vroom(pipe(awk_cmd), col_names = c("Nuc", "controls"), delim = "\t", show_col_types = FALSE) %>%
    filter(Nuc %in% nucs)
  return(df)
}

#' Get genome-wide counts of nucleotides flanking A or C
#' 
#' @param nucleotide Either A, cpg, or non
#' @param rp Flanking position at which we stratify the genome-wide count of A or C
#' @return data.frame with genome-wide counts of nucleotides flanking focal site
get_gw_position <- function(nucleotide, rp){
  
  
  data_dir <- "/net/snowwhite/home/beckandy/research/BRIDGES_redo/output/gw_1_count/"
  
  if(nucleotide == "A"){
    f_name <- paste0(data_dir, nucleotide, "_rp", rp, ".csv")
  } else {
    f_name <- paste0(data_dir, "cpg/" , nucleotide, "_rp", rp, ".csv")
  }
  df <- read_csv(f_name, col_types = cols())
  colnames(df) <- c("Nuc", "n_gw")
  return(df)
  
}

#' Get singleton, control, and genome-wide counts for a sub-type at a flanking position
#' 
#' @param subtype One of the nine basic sub-types
get_all_position
library(tidyverse)
args = commandArgs(trailingOnly=TRUE) # Expected: number 1 - 9
subtype_ix <- as.numeric(args[1])

subtypes <- c("AT_CG", "AT_GC", "AT_TA",
              "GC_AT", "GC_TA", "GC_CG", 
              "cpg_GC_AT", "cpg_GC_TA", "cpg_GC_CG")
subtype <- subtypes[subtype_ix]

print(paste0("Subtype = ", subtype))

#############################################################
# FUNCTIONS
#############################################################

#' Get the singleton counts for a sub-type at a flanking position
#' 
#' @param subtype One of the nine simple sub types
#' @param rp1 First relative position which to get singleton counts from
#' @param rp1 Second relative position which to get singleton counts from
#' @return data frame with counts for each nucleotide
get_singleton_counts <- function(subtype, rp1, rp2){
  nucs <- c("A", "C", "G", "T")
  
  pos1 <- rp1 + 11
  pos2 <- rp2 + 11
  
  singleton_dir <- "/net/snowwhite/home/beckandy/research/BRIDGES_redo/output/singletons/"
  fName <- paste0(singleton_dir, subtype, ".txt")
  
  awk_cmd <- paste0("awk '{count[substr($1,", pos1,", 1)substr($1,", pos2,", 1)]++}END{for(key in count) print(key, \"\\t\", count[key])}' ",
                    fName)
  
  df <- vroom::vroom(pipe(awk_cmd), col_names = c("Nuc", "singletons"), delim = "\t", show_col_types = FALSE) %>%
    mutate(p1 = str_sub(Nuc, 1, 1), p2 = str_sub(Nuc, 2, 2)) %>%
    filter(p1 %in% nucs, p2 %in% nucs) %>%
    select(p1, p2, singletons) %>%
    arrange(p1, p2)
  
  return(df)
}

#' Get nucleotide counts at position pairs relative to control observation
#' 
#' @param subtype One of the nine basic subtypes
#' @param rp First Position relative to focal site to get counts for
#' @param rp Second Position relative to focal site to get counts for
#' @return Data.frame with control counts stratified by nucleotide at flanking position
get_control_counts <- function(subtype, rp1, rp2){
  nucs <- c("A", "C", "G", "T")
  
  pos1 <- rp1 + 11
  pos2 <- rp2 + 11
  
  control_dir <- "/net/snowwhite/home/beckandy/research/BRIDGES_redo/output/controls/"

  fName <- paste0(control_dir, subtype, ".txt")
  
  awk_cmd <- paste0("awk '{count[substr($1,", pos1,", 1)substr($1,", pos2,", 1)]++}END{for(key in count) print(key, \"\\t\", count[key])}' ",
                    fName)
  
  df <- vroom::vroom(pipe(awk_cmd), col_names = c("Nuc", "controls"), delim = "\t", show_col_types = FALSE) %>%
    mutate(p1 = str_sub(Nuc, 1, 1), p2 = str_sub(Nuc, 2, 2)) %>%
    filter(p1 %in% nucs, p2 %in% nucs) %>%
    select(p1, p2, controls) %>%
    arrange(p1, p2)
  return(df)
}

#' Get genomewide counts for a given subtype
#' 
#' @param subtype Subtype
#' @param p1 First relative position
#' @param p2 Second relative position (p2 > p1)
#' @return Data.frame with nucleotides at two positions and their counts
# get_gw_counts <- function(subtype, p1, p2){
#   
#   input_dir <- "/net/snowwhite/home/beckandy/research/BRIDGES_redo/output/gw_2_count/"
#   if(str_starts(subtype, "cpg")){
#     prefix <- "cpg_GC"
#   } else if(str_starts(subtype, "GC")){
#     prefix <- "GC"
#   } else {
#     prefix <- "AT"
#   }
#   
#   file_name <- paste0(input_dir, prefix, "_p", p1, "_q", p2, ".csv")
#   df <- read_csv(file_name, col_types = cols()) %>%
#     rename(n_gw = n)
#   
#   return(df)
#   
# }

#' Get all counts for a given sub-type
#' 
#' @param subtype Sub-type
#' @param p1 First relative position
#' @param p2 Second relative position (p2 > p1)
#' @return Data.frame with nucleotides at two positions and their counts
get_all_counts <- function(subtype, p1, p2){
  df_s <- get_singleton_counts(subtype, p1, p2)
  df_c <- get_control_counts(subtype, p1, p2)
  #df_gw <- get_gw_counts(subtype, p1, p2)
  
  df <- full_join(df_s, df_c, by = c("p1", "p2")) %>%
    replace_na(list("singletons" = 0, "controls" = 0))
    #full_join(df_gw, by = c("p1", "p2")) %>%
    #replace_na(list("singletons" = 0, "controls" = 0))
  
  return(df)
}

######################################################################
# GET COUNTS ACROSS ALL PAIRS OF POSITIONS
######################################################################

out_dir <- "/net/snowwhite/home/beckandy/research/BRIDGES_redo/output/all_count_2_pos/"

for(p1 in c(-10:-1,1:9)){
  for(p2 in (p1 + 1):10){
    if(p2 == 0) next
    df <- get_all_counts(subtype, p1, p2)
    write_csv(df, paste0(out_dir, subtype, "_p", p1, "_q", p2, ".csv"))
  }
}

print("Done!")
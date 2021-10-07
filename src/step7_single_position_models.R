# Gather data for the single position models

library(tidyverse)
library(foreach)
library(doParallel)
library(writexl)

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
get_control_new_pos <- function(subtype, rp){
  nucs <- c("A", "C", "G", "T")
  
  pos <- rp + 11
  
  if(str_starts(subtype, "AT")){
    control_dir <- "/net/snowwhite/home/beckandy/research/BRIDGES_redo/output/control2/at/"
  } else {
    control_dir <- "/net/snowwhite/home/beckandy/research/BRIDGES_redo/output/control2/gc/"
  }
  fName <- paste0(control_dir, subtype, ".txt")
  
  awk_cmd <- paste0("awk '{count[substr($1,", pos,", 1)]++}END{for(key in count) print(key, \"\\t\", count[key])}' ",
                    fName)
  
  df <- vroom::vroom(pipe(awk_cmd), col_names = c("Nuc", "controls"), delim = "\t", show_col_types = FALSE) %>%
    filter(Nuc %in% nucs) %>%
    rename(controls_2 = controls)
  return(df)
}

#' Get nucleotide counts at position relative to control observation
#' 
#' @param subtype One of the nine basic subtypes
#' @param rp Position relative to focal site to get counts for
#' @return Data.frame with control counts stratified by nucleotide at flanking position
get_control_old_pos <- function(subtype, rp){
  nucs <- c("A", "C", "G", "T")
  
  pos <- rp + 11
  
  if(str_starts(subtype, "AT")){
    control_dir <- "/net/snowwhite/home/beckandy/research/BRIDGES_redo/output/control/at/"
  } else {
    control_dir <- "/net/snowwhite/home/beckandy/research/BRIDGES_redo/output/control/gc/"
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
#' @param rp Flanking position at which we stratify the genome-wide count of A or C
#' @return data.frame with singleton, control, and genome-wide counts of nucleotides flanking focal site
get_all_position <- function(subtype, rp){
  if(str_starts(subtype, "AT")){
    nuc <- "A"
  } else if(str_starts(subtype, "GC")){
    nuc <- "non_c"
  } else {
    nuc <- "cpg_c"
  }
  
  df_s <- get_singleton_pos(subtype, rp)
  df_c <- get_control_old_pos(subtype, rp)
  df_c_2 <- get_control_new_pos(subtype, rp)
  df_gw <- get_gw_position(nuc, rp)
  
  df <- full_join(df_s, df_c, by = "Nuc") %>%
    full_join(df_c_2, by = "Nuc") %>%
    full_join(df_gw, by = "Nuc") %>%
    replace_na(list("singletons" = 0, "controls" = 0, "controls_2" = 0, "n_gw" = 0)) %>%
    mutate(pct_gw = n_gw / sum(n_gw),
           pct_c = controls / sum(controls),
           pct_c_2 = controls_2 / sum(controls_2),
           pct_s = singletons / sum(singletons)) %>%
    mutate(exp_gw = sum(singletons) * pct_gw,
           exp_ct = sum(singletons) * pct_c,
           exp_ct_2 = sum(singletons) * pct_c_2) %>%
    mutate(chi_sq_gw = ((singletons - exp_gw)^2)/exp_gw,
           chi_sq_ct = ((singletons - exp_ct)^2)/exp_ct,
           chi_sq_ct_2 = ((singletons - exp_ct_2)^2)/exp_ct_2)
  
  df$rp <- rp
  return(df)
  
}

subtypes <- c("AT_CG", "AT_GC", "AT_TA",
              "GC_AT", "GC_TA", "GC_CG",
              "cpg_GC_AT", "cpg_GC_TA", "cpg_GC_CG")
results <- vector("list", length(subtypes))
names(results) <- subtypes

registerDoParallel(20)
for(st in subtypes){
  message(st)
  results[[st]] <- foreach(x = c(-10:-1,1:10), .combine = 'rbind') %dopar% {
    get_all_position(st, x)
  }
}
stopImplicitCluster()

out_dir <- "/net/snowwhite/home/beckandy/research/BRIDGES_redo/output/single_pos_df/"
for(st in subtypes){
  message(st)
  df <- results[[st]]
  for(i in c(-10:-1,1:10)){
    df2 <- df %>%
      filter(rp == i) %>%
      select(-rp)
    out_file <- paste0(out_dir, st, "_rp", i, ".csv")
    write_csv(df2, out_file)
  }
}

excel_results <- vector("list", length(subtypes))
names(excel_results) <- subtypes
for(st in subtypes){
  df <- results[[st]]
  df_f <- df %>%
    filter(rp == i) %>%
    select(-rp) %>%
    select(Nuc, singletons, controls, n_gw, chi_sq_gw, chi_sq_ct) %>%
    gather(category, n, -Nuc)
  colnames(df_f) <- c("Nuc", "category", "n-10")
  for(i in c(-9:-1,1:10)){
    df2 <- df %>%
      filter(rp == i) %>%
      select(-rp) %>%
      select(Nuc, singletons, controls, n_gw, chi_sq_gw, chi_sq_ct) %>%
      gather(category, n, -Nuc)
    colnames(df2) <- c("Nuc", "category", paste0("n",i))
    df_f <- full_join(df_f, df2, by = c("Nuc", "category")) %>%
      mutate_all(~replace(., is.na(.), 0))
  }
  excel_results[[st]] <- df_f
}

write_xlsx(excel_results, paste0(out_dir, "all_sp_bridges.xlsx"))


library(tidyverse)

data_dir <- "/net/snowwhite/home/beckandy/research/BRIDGES_redo/output/gw_1_count/cpg/"

load_cpg_nuc_chrom_rp <- function(nuc, rp, cpg_stat, chrom,  data_dir = "/net/snowwhite/home/beckandy/research/BRIDGES_redo/output/gw_1_count/cpg/"){
  
  f_name <- paste0(data_dir, cpg_stat ,"_", nuc, "_chrom", chrom, "_rp", rp, ".csv")
  df <- read_csv(f_name, col_types = cols())
  
  return(df)
}

load_cpg_nuc_rp_all <- function(nuc, rp, cpg_stat) {
  nucs <- c("A", "C", "G", "T")
  df <- load_cpg_nuc_chrom_rp(nuc, rp, cpg_stat, 1)
  
  for(i in c(2:22)){
    df2 <- load_cpg_nuc_chrom_rp(nuc, rp, cpg_stat, i)
    
    df <- full_join(df, df2, by = "Nuc") %>%
      replace_na(list("n.x" = 0, "n.y" = 0)) %>%
      mutate(n = n.x + n.y) %>%
      select(Nuc, n)
  }
  
  df <- df %>%
    filter(Nuc %in% nucs) %>%
    arrange(Nuc)
  
  return(df)
}

for(status in c("cpg", "non")){
  for(nuc in c("c", "g")){
    for(rp in c(-10:-1,1:10)){
      df <- load_cpg_nuc_rp_all(nuc, rp, status)
      write_csv(df, paste0(data_dir, status, "_rp", rp, ".csv"))
    }
  }
}
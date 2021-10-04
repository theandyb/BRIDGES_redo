library(tidyverse)

input_dir <- "/net/snowwhite/home/beckandy/research/BRIDGES_redo/output/gw_1_count"


load_df <- function(nuc, chr, rp){
  df <- read_csv(paste0(input_dir,"/chr", chr, "_", nuc, "_rp", rp, ".csv"), col_types = cols())
  return(df)
}

get_rp <- function(nuc, rp){
  df <- load_df(nuc, 1, rp)
  
  for(chrom in c(2:22)){
    df2 <- load_df(nuc, chrom, rp)
    df <- df %>%
      full_join(df2, by="Nuc") %>%
      replace_na(list("N.x" = 0, "N.y" = 0)) %>%
      mutate(N = N.x + N.y) %>%
      select(Nuc, N)
  }
  
  return(df)
}

for(rp in c(-10:-1, 1:10)){
  print(rp)
  df <- get_rp("A", rp)
  write_csv(df, paste0(input_dir, "/A_rp", rp, ".csv"))
}

for(rp in c(-10:-1, 1:10)){
  print(rp)
  df <- get_rp("C", rp)
  write_csv(df, paste0(input_dir, "/C_rp", rp, ".csv"))
}

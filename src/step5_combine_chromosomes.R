library(tidyverse)

data_dir <- "/net/snowwhite/home/beckandy/research/BRIDGES_redo/output/gw_1_count/3cat"

for(nuc in c("A", "C", "cpg_C")){
  for(rp in c(-10:-1,1:10)){
    f_name <- paste0(data_dir, "/chr", 1, "_", nuc, "_rp", rp, ".csv")
    out_name <- paste0(data_dir, "/", nuc, "_rp", rp, ".csv")
    df <- read_csv(f_name, col_types = cols())
    for(i in c(2:22)){
      f_name2 <- paste0(data_dir, "/chr", i, "_", nuc, "_rp", rp, ".csv")
      df2 <- read_csv(f_name2, col_types = cols())
      df %>% full_join(df, df2, by = "Nuc") %>%
        replace_na(list("N.x" = 0, "N.y" = 0)) %>%
        mutate(N = N.x + N.y) %>%
        select(Nuc, N)
    }
    write_csv(df, out_name)
  }
}

data_dir <- "/net/snowwhite/home/beckandy/research/BRIDGES_redo/output/gw_1_count/6cat"

for(nuc in c("A", "C", "G", "T", "cpg_C", "cpg_G")){
  for(rp in c(-10:-1,1:10)){
    f_name <- paste0(data_dir, "/chr", 1, "_", nuc, "_rp", rp, ".csv")
    out_name <- paste0(data_dir, "/", nuc, "_rp", rp, ".csv")
    df <- read_csv(f_name, col_types = cols())
    for(i in c(2:22)){
      f_name2 <- paste0(data_dir, "/chr", i, "_", nuc, "_rp", rp, ".csv")
      df2 <- read_csv(f_name2, col_types = cols())
      df %>% full_join(df, df2, by = "Nuc") %>%
        replace_na(list("N.x" = 0, "N.y" = 0)) %>%
        mutate(N = N.x + N.y) %>%
        select(Nuc, N)
    }
    write_csv(df, out_name)
  }
}
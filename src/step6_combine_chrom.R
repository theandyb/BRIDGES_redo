library(tidyverse)

data_dir <- "/net/snowwhite/home/beckandy/research/BRIDGES_redo/output/gw_2_count/6cats/"

for(nuc in c("a", "t", "non_c", "non_g", "cpg_c", "cpg_g")){
  print(nuc)
  for(p in c(-10:-1, 1:9)){
    for(q in (p+1):10){
      if(q == 0) next
      print(paste0("p: ", p, "; q: ", q))
      fname <- paste0(data_dir, nuc, "_chrom", 1, "_p", p, "_q", q, ".csv")
      df <- read_csv(fname, col_types = cols())

      for(chrom in c(2:22)){
        fname2 <- paste0(data_dir, nuc, "_chrom", chrom, "_p", p, "_q", q, ".csv")
        df2 <- read_csv(fname, col_types = cols())

        df <- full_join(df, df2, by = "Nucs") %>%
          replace_na(list("n.x" = 0, "n.y" = 0)) %>%
          mutate(n = n.x + n.y) %>%
          select(Nucs, n)
      }

      oname <- paste0(data_dir, nuc, "_p", p, "_q", q, ".csv")
      write_csv(df, oname)

    }
  }
}

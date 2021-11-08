library(tidyverse)

data_dir <- "/net/snowwhite/home/beckandy/research/BRIDGES_redo/output/gw_2_count/6cats/"
out_dir <- "/net/snowwhite/home/beckandy/research/BRIDGES_redo/output/gw_2_count/3cats/"

rc <- function(z){
  rc1 <- function(zz){
    s <- strsplit(zz, split = "")[[1]]
    s <- rev(s)
    dchars <- strsplit("ACGTMRWSYKVHDBNI", split = "")[[1]]
    comps <- strsplit("TGCAKYWSRMBDHVNI", split = "")[[1]]
    s <- s[s %in% dchars] # remove spaces etc
    s <- dchars[match(s, comps)]
    s <- paste0(s, collapse = "")
    return(s)
  }
  z <- toupper(z)
  tmpnames <- names(z)
  res <- unname(sapply(z, rc1))
  if(!is.null(attr(z, "quality"))){
    strev <- function(x) sapply(lapply(lapply(unname(x), charToRaw), rev), rawToChar)
    attr(res, "quality") <- unname(sapply(attr(z, "quality"), strev))
  }
  names(res) <- tmpnames
  return(res)
}

# A and T
for(i in c(-10:-1, 1:9)){
  for(j in (i+1):10){
    if(j == 0) next
    a_name <- paste0(data_dir, "a_p", i, "_q", j, ".csv")
    t_name <- paste0(data_dir, "t_p", -j, "_q", -i, ".csv")
    a_df <- read_csv(a_name, col_types = cols())
    t_df <- read_csv(t_name, col_types = cols()) %>%
      mutate(Nucs = rc(Nucs))

    at_df <- full_join(a_df, t_df, by="Nucs") %>%
      replace_na(list("n.x" = 0, "n.y" = 0)) %>%
      mutate(n = n.x + n.y) %>%
      select(Nucs, n)

    out_name <- paste0(out_dir, "at_p", i, "_q", j, ".csv")
    write_csv(at_df, out_name)
  }
}

# C and G
for(i in c(-10:-1, 1:9)){
  for(j in (i+1):10){
    if(j == 0) next
    c_name <- paste0(data_dir, "c_p", i, "_q", j, ".csv")
    g_name <- paste0(data_dir, "g_p", -j, "_q", -i, ".csv")
    c_df <- read_csv(c_name, col_types = cols())
    g_df <- read_csv(g_name, col_types = cols()) %>%
      mutate(Nucs = rc(Nucs))

    cg_df <- full_join(c_df, g_df, by="Nucs") %>%
      replace_na(list("n.x" = 0, "n.y" = 0)) %>%
      mutate(n = n.x + n.y) %>%
      select(Nucs, n)

    out_name <- paste0(out_dir, "gc_p", i, "_q", j, ".csv")
    write_csv(cg_df, out_name)
  }
}

# CPG C and G
for(i in c(-10:-1, 1:9)){
  for(j in (i+1):10){
    if(j == 0) next
    c_name <- paste0(data_dir, "cpg_c_p", i, "_q", j, ".csv")
    g_name <- paste0(data_dir, "cpg_g_p", -j, "_q", -i, ".csv")
    c_df <- read_csv(c_name, col_types = cols())
    g_df <- read_csv(g_name, col_types = cols()) %>%
      mutate(Nucs = rc(Nucs))

    cg_df <- full_join(c_df, g_df, by="Nucs") %>%
      replace_na(list("n.x" = 0, "n.y" = 0)) %>%
      mutate(n = n.x + n.y) %>%
      select(Nucs, n)

    out_name <- paste0(out_dir, "cpg_gc_p", i, "_q", j, ".csv")
    write_csv(cg_df, out_name)
  }
}

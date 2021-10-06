library(tidyverse)

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

args = commandArgs(trailingOnly=TRUE)
fname <- paste0(args[1], ".csv")
oname <- paste0(args[1], "_rc.csv")

print(paste0("Input: ", fname))
print(paste0("Output: ", oname))

df <- read_csv(fname,
              col_names = c("chrom", "pos", "motif", "subtype", "alt","window","distance"),
              col_types = cols())
df <- df %>% 
    mutate(center = str_sub(motif,11,11)) %>%
    mutate(motif2 = ifelse(center %in% c("A","C"), motif, rc(motif))) %>%
    select(chrom, pos, motif2, subtype, window, distance) %>%
    rename(motif = motif2)

write_csv(df, oname)
print("Done!")
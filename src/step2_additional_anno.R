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
oname <- paste0(args[1], "_addl.csv")

print(paste0("Input: ", fname))
print(paste0("Output: ", oname))

df <- read_csv(fname,
              col_names = c("chrom", "pos", "motif", "subtype", "alt","id"),
              col_types = cols(
                  chrom = col_character(),
                  pos = col_double(),
                  motif = col_character(),
                  subtype = col_character(),
                  alt = col_character(),
                  id = col_character()
                ))
df <- df %>% mutate(REF = str_sub(subtype, 1, 1)) %>% 
    mutate(ALT = str_sub(subtype, 3, 3)) %>%
    mutate(subtype2 = ifelse(REF %in% c("A","G"), 
                             paste0(REF, rc(REF), ">", ALT,rc(ALT)), 
                             paste0(rc(REF), REF, ">", rc(ALT),ALT))) %>%
    mutate(motif2 = ifelse(REF %in% c("A","C"), 
                           paste0(motif,"(", rc(motif), ")"), 
                           paste0(rc(motif), "(", motif, ")"))) %>% 
    mutate(subtype2 = ifelse(str_sub(motif2,11,12) == "CG", paste0("cpg_", subtype2), subtype2))

write_csv(df, oname)
print("Done!")
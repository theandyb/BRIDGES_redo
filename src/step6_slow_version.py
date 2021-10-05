from pyfaidx import Fasta
import pandas as pd
import sys
ref_file = "/net/snowwhite/home/beckandy/research/BRIDGES_redo/reference_data/human_g1k_v37.fasta"

def count_string(nuc_string, window, result_tables):
  rc_dict = {"A":"T", "C":"G", "G":"C", "T":"A"}
  nucs = ["A", "C", "G", "T"]
  for i in range(len(nuc_string)):
      nuc = nuc_string[i]
      if i % 1000000 == 0:
        print(i)
      if not(nuc in nucs):
          continue
      if nuc in ["A", "C"]:
          direction = 1
      else:
          direction = -1
          nuc = rc_dict[nuc]
      for j in range(-window,(window)):
          if j == 0: 
              continue
          if j + i < 0: 
              continue
          if j + i >= len(nuc_string): 
              continue
          new_nuc_j = nuc_string[(i+j)]
          if not(new_nuc_j in nucs):
              continue
          if direction == -1:
              new_nuc_j = rc_dict[new_nuc_j]
          ix_j = j * direction
          for k in range(j+1,(window+1)):
            if k == 0:
              continue
            if i + k >= len(nuc_string):
              continue
            new_nuc_k = nuc_string[(i+k)]
            if not(new_nuc_k in nucs):
              continue
            if direction == -1:
              new_nuc_k = rc_dict[new_nuc_k]
            ix_k = k * direction
            result_tables[nuc][ix_j + window][ix_k + window]["{}{}".format(new_nuc_j,new_nuc_k)] += 1
  return 0

def table_df(tables, nuc, ix_j, ix_k):
  df = pd.DataFrame.from_dict(tables[nuc][ix_j][ix_k], orient='index')
  df.index.name = 'Nuc'
  df.reset_index(inplace=True)
  df = df.rename({0:'N'}, axis = 1)
  return df

def main(chrom):
  ref_file = "/net/snowwhite/home/beckandy/research/BRIDGES_redo/reference_data/human_g1k_v37.fasta"
  out_dir = "/net/snowwhite/home/beckandy/research/BRIDGES_redo/output/gw_2_count/"
  fasta_obj = Fasta(ref_file)
  window = 10 # how many bases up/down are counting flanking positions?
  
  empty_table = {"AA": 0, "AC": 0, "AG": 0, "AT": 0, \
  "CA": 0, "CC": 0, "CG": 0, "CT": 0, \
  "GA": 0, "GC": 0, "GG": 0, "GT": 0, \
  "TA": 0, "TC": 0, "TG": 0, "TT": 0,}
  
  result_tables = {"A": [[empty_table for i in range(((window*2)+1))] for i in range(((window * 2)) + 1)],
  "C": [[empty_table for i in range(((window*2)+1))] for i in range(((window * 2)) + 1)]}
  
  print("Chromosome {}...".format(chrom))
  seq = fasta_obj["{}".format(chrom)]
  seqstr = seq[0:len(seq)].seq
  count_string(seqstr, window, result_tables)
  print("Done!")
    
  print("Writing files...")
  for i in range(-window,(window)):
    if i == 0:
        continue
    for j in range((i+1), (window+1)):
      if j == 0:
        continue
      ix_i = i + window
      ix_j = j + window
      df_a = table_df(result_tables, "A", ix_i, ix_j)
      df_c = table_df(result_tables, "C", ix_i, ix_j)
      a_file = "{}chr{}_A_rp{}_{}.csv".format(out_dir, chrom, i, j)
      c_file = "{}chr{}_C_rp{}_{}.csv".format(out_dir, chrom, i, j)
      df_a.to_csv(a_file, index=False)
      df_c.to_csv(c_file, index=False)
  
  return 0

if __name__ == "__main__":
  chrom = sys.argv[1]
  main(chrom)

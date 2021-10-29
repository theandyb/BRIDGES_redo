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
      if nuc_string[i:(i+2)] == "CG":
        cat_prefix = "cpg_"
      elif i > 0:
        if nuc_string[(i-1):(i+1)] == "CG":
          cat_prefix = "cpg_"
        else:
          cat_prefix = ""
      else:
        cat_prefix = ""
      direction = 1
      #if nuc in ["A", "C"]:
      #    direction = 1
      #else:
      #    direction = -1
      #    nuc = rc_dict[nuc]
      for rp in range(-window,(window+1)):
          if rp == 0: 
              continue
          if rp + i < 0:
            continue
          if rp + i >= len(nuc_string):
            continue
          new_nuc = nuc_string[(i+rp)]
          if not(new_nuc in nucs):
              continue
          if direction == -1:
              new_nuc = rc_dict[new_nuc]
          ix = rp * direction
          try:
            result_tables["{}{}".format(cat_prefix, nuc)][ix + window][new_nuc] += 1
          except KeyError:
            print("Bad table name!")
            print("Offending nucleotide: {}".format(new_nuc))
  return 0

def table_df(tables, nuc, ix):
  df = pd.DataFrame.from_dict(tables[nuc][ix], orient='index')
  df.index.name = 'Nuc'
  df.reset_index(inplace=True)
  df = df.rename({0:'N'}, axis = 1)
  return df

def main(chrom, chrom_prefix = ""):
  ref_file = "/net/snowwhite/home/beckandy/research/BRIDGES_redo/reference_data/human_g1k_v37.fasta"
  out_dir = "/net/snowwhite/home/beckandy/research/BRIDGES_redo/output/gw_1_count/6cat/"
  fasta_obj = Fasta(ref_file)
  window = 10 # how many bases up/down are counting flanking positions?
  
  result_tables = {"A": [{"A":0, "C": 0, "G": 0, "T": 0} for i in range(((window*2)+1))],
             "C": [{"A":0, "C": 0, "G": 0, "T": 0} for i in range(((window*2)+1))],
             "G": [{"A":0, "C": 0, "G": 0, "T": 0} for i in range(((window*2)+1))],
             "T": [{"A":0, "C": 0, "G": 0, "T": 0} for i in range(((window*2)+1))],
             "cpg_C": [{"A":0, "C": 0, "G": 0, "T": 0} for i in range(((window*2)+1))],
             "cpg_G":  [{"A":0, "C": 0, "G": 0, "T": 0} for i in range(((window*2)+1))]}
  
  print("Chromosome {}...".format(chrom))
  seq = fasta_obj["{}{}".format(chrom_prefix, chrom)]
  seqstr = seq[0:len(seq)].seq
  count_string(seqstr, window, result_tables)
  print("Done!")
    
  print("Writing files...")
  for rp in range(-window,(window+1)):
    if rp == 0:
        continue
    ix = rp + window
    df_a = table_df(result_tables, "A", ix)
    df_c = table_df(result_tables, "C", ix)
    df_g = table_df(result_tables, "G", ix)
    df_t = table_df(result_tables, "T", ix)
    df_cpg_c = table_df(result_tables, "cpg_C", ix)
    df_cpg_g = table_df(result_tables, "cpg_G", ix)
    
    a_file = "{}chr{}_A_rp{}.csv".format(out_dir, chrom, rp)
    c_file = "{}chr{}_C_rp{}.csv".format(out_dir, chrom, rp)
    g_file = "{}chr{}_G_rp{}.csv".format(out_dir, chrom, rp)
    t_file = "{}chr{}_T_rp{}.csv".format(out_dir, chrom, rp)
    cpg_c_file = "{}chr{}_cpg_C_rp{}.csv".format(out_dir, chrom, rp)
    cpg_g_file = "{}chr{}_cpg_G_rp{}.csv".format(out_dir, chrom, rp)
    
    df_a.to_csv(a_file, index=False)
    df_c.to_csv(c_file, index=False)
    df_g.to_csv(g_file, index=False)
    df_t.to_csv(t_file, index=False)
    df_cpg_c.to_csv(cpg_c_file, index=False)
    df_cpg_g.to_csv(cpg_g_file, index=False)
    
  
  return 0

if __name__ == "__main__":
  chrom = sys.argv[1]
  main(chrom)

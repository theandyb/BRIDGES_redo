from pyfaidx import Fasta
import pandas as pd
from Bio.Seq import Seq
import argparse
import numpy as np

chroms = [x for x in range(1,23)]
p1 = [x for x in range(-10,0)] + [x for x in range(1,10)]
pairs = {}

for c in chroms:
  for i in p1:
    for j in range((i+1),11):
      if j == 0: continue
      pairs['{}_{}_{}'.format(c, i, j)] = 0

keys = list(pairs)

def gen_empty_dict(nucs =  ["A","C","G","T"]):
  """Create dict of dicts in which count results will be stored"""
  combs = ['AA','AC', 'AG', 'AT','CA','CC', 'CG', 'CT','GA','GC', 'GG', 'GT','TA','TC', 'TG', 'TT']
  result = {}
  for y in nucs:
    result[y] = {}
    for c in combs:
        result[y][c] = 0
  return result

def main(chrom_prefix = ""):
  """Generate genome-wide counts of nucleotides at postions relative to CpGs and C/Gs in non-CpG contexts"""
  
  parser = argparse.ArgumentParser(description="Genome-wide counts for CpG/non-CpG single position models.")
  parser.add_argument("-j", "--jobid", help="which of the 4180 jobs", required=True)
  parser.add_argument("-o", "--output", help="Path to output", required=True)
  args = parser.parse_args()
  
  nucs = ['A', 'C', 'G', 'T']
  
  ref = Fasta("/net/snowwhite/home/beckandy/research/BRIDGES_redo/reference_data/human_g1k_v37.fasta", as_raw=True)
  key_id = int(args.jobid) - 1 # 0 index :/
  out_dir = args.output
  
  job_list = keys[key_id].split("_")
  chrom = job_list[0]
  p1 = int(job_list[1])
  p2 = int(job_list[2])
  
  # Load in the chromosome's sequence as a string
  chr_seq = np.asarray(ref['{}{}'.format(chrom_prefix, chrom)])
  non_cpg = gen_empty_dict()
  cpg_res = gen_empty_dict(['C', 'G'])
  
  end = len(chr_seq)
  for i in range(0, end):
    if i % 1000000 == 0:
      print(i)
    if p1 + i < 0 or p2 + i < 0: continue
    if p1 + i >= end or p2 + i >= end: break
    index_nuc = str(chr_seq[i], encoding='ascii')
    p1_nuc = str(chr_seq[i+p1], encoding='ascii')
    p2_nuc = str(chr_seq[i+p2], encoding='ascii')
    if index_nuc not in nucs: continue
    if p1_nuc not in nucs: continue
    if p2_nuc not in nucs: continue
    if str(chr_seq[i:i+2], encoding='ascii') == "CG" or str(chr_seq[i-1:i+1], encoding='ascii') == "CG":
      cpg_bool = True
      cpg_res[index_nuc]["{}{}".format(p1_nuc, p2_nuc)] += 1
    else:
      cpg_bool = False
      non_cpg[index_nuc]["{}{}".format(p1_nuc, p2_nuc)] += 1
    
  # Write outputs
  cpg_c_df = pd.DataFrame.from_dict(cpg_res["C"], orient='index', columns=['n'])
  cpg_g_df = pd.DataFrame.from_dict(cpg_res["G"], orient='index', columns=['n'])
  a_df = pd.DataFrame.from_dict(non_cpg["A"], orient='index', columns=['n'])
  c_df = pd.DataFrame.from_dict(non_cpg["C"], orient='index', columns=['n'])
  g_df = pd.DataFrame.from_dict(non_cpg["G"], orient='index', columns=['n'])
  t_df = pd.DataFrame.from_dict(non_cpg["T"], orient='index', columns=['n'])
  
  out_file_cpg_c = out_dir + "/cpg_c_chrom{}_p{}_q{}.csv".format(chrom, p1, p2)
  out_file_non_c = out_dir + "/non_c_chrom{}_p{}_q{}.csv".format(chrom, p1, p2)
  out_file_cpg_g = out_dir + "/cpg_g_chrom{}_p{}_q{}.csv".format(chrom, p1, p2)
  out_file_non_g = out_dir + "/non_g_chrom{}_p{}_q{}.csv".format(chrom, p1, p2)
  out_file_a = out_dir + "/a_chrom{}_p{}_q{}.csv".format(chrom, p1, p2)
  out_file_t = out_dir + "/t_chrom{}_p{}_q{}.csv".format(chrom, p1, p2)
  
  cpg_c_df.to_csv(out_file_cpg_c, index_label="Nucs")
  c_df.to_csv(out_file_non_c, index_label="Nucs")
  cpg_g_df.to_csv(out_file_cpg_g, index_label="Nucs")
  g_df.to_csv(out_file_non_g, index_label="Nucs")
  a_df.to_csv(out_file_a, index_label="Nucs")
  t_df.to_csv(out_file_t, index_label="Nucs")
    
  return 0

if __name__ == "__main__":
    main()

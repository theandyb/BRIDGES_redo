from pyfaidx import Fasta
import pandas as pd
from Bio.Seq import Seq
import argparse

def gen_empty_dict():
    """Create dict of dicts in which count results will be stored"""
    result = {}
    for y in ["C","G"]:
        result[y] = {}
        for x in range(-10,11):
            if x==0: continue
            result[y][x] = {}
    return result

def main():
    """Generate genome-wide counts of nucleotides at postions relative to CpGs and C/Gs in non-CpG contexts"""
    
    parser = argparse.ArgumentParser(description="Genome-wide counts for CpG/non-CpG single position models.")
    parser.add_argument("-c", "--chromosome", help="Chromsome to get count from", required=True)
    parser.add_argument("-o", "--output", help="Path to output", required=True)
    args = parser.parse_args()
    
    ref = Fasta("/net/snowwhite/home/beckandy/research/BRIDGES_redo/reference_data/human_g1k_v37.fasta")
    chrom = args.chromosome
    out_dir = args.output
    
    # Load in the chromosome's sequence as a string
    chr_seq = ref['{}'.format(chrom)]
    chr_string = chr_seq[0:len(chr_seq)].seq
    
    #chr_string = "ACGGCTGGCNNNCGNATACCCGCGCGCGCGCGCGCGACCT"
    # Object to save our results
    cpg_res = gen_empty_dict()
    non_res = gen_empty_dict()
    # Traverse the chromosome!
    start = 0
    end = len(chr_string)
    for i in range(start, end):
        nuc = chr_string[i]
        if nuc != "C" and nuc != "G":
            continue
        # Grab 21-mer centered at current position
        # Handle chromosome 17 edge case
        if i < 10:
            pad_left = -i
        else:
            pad_left = -10
        if i + 11 >= end:
            pad_right = end - i - 1 # how many positions do we have available upstream of the nucleotide of interest?
        else:
            pad_right = 10
        motif_i = -pad_left # index of the current nucleotide in motif we're grabbing, which will usually be 21-mer (catching edge cases here)
        motif = chr_string[(i+pad_left):(i+pad_right+1)]
        cpg_bool = motif[motif_i:(motif_i+2)] == "CG" or motif[(motif_i-1):(motif_i+1)] == "CG"
        # Now add counts to corresponding nucleotide counts for each relative position!
        for pos in range(0, len(motif)):
            rp = pos - motif_i
            if rp == 0: continue
            x = motif[pos]
            if cpg_bool:
                if x in cpg_res[nuc][rp]:
                    cpg_res[nuc][rp][x] += 1
                else:
                    cpg_res[nuc][rp][x] = 1
            else:
                if x in non_res[nuc][rp]:
                    non_res[nuc][rp][x] += 1
                else:
                    non_res[nuc][rp][x] = 1
    
    # Save off results: table per relative position
    for key in cpg_res["C"]:
        cpg_df = pd.DataFrame.from_dict(cpg_res["C"][key], orient='index', columns=['n'])
        non_df = pd.DataFrame.from_dict(non_res["C"][key], orient='index', columns=['n'])
        out_file_cpg = out_dir + "/cpg_c_chrom{}_rp{}.csv".format(chrom, key)
        out_file_non = out_dir + "/non_c_chrom{}_rp{}.csv".format(chrom, key)
        cpg_df.to_csv(out_file_cpg, index_label="Nuc")
        non_df.to_csv(out_file_non, index_label="Nuc")
    for key in cpg_res["G"]:
        cpg_df = pd.DataFrame.from_dict(cpg_res["G"][key], orient='index', columns=['n'])
        non_df = pd.DataFrame.from_dict(non_res["G"][key], orient='index', columns=['n'])
        out_file_cpg = out_dir + "/cpg_g_chrom{}_rp{}.csv".format(chrom, key)
        out_file_non = out_dir + "/non_g_chrom{}_rp{}.csv".format(chrom, key)
        cpg_df.to_csv(out_file_cpg, index_label="Nuc")
        non_df.to_csv(out_file_non, index_label="Nuc")
    
if __name__ == "__main__":
    main()

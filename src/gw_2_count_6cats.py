from pyfaidx import Fasta
import pandas as pd
from Bio.Seq import Seq
import argparse

def gen_empty_dict(nucs =  ["A","C","G","T"]):
    """Create dict of dicts in which count results will be stored"""
    combs = ['AA','AC', 'AG', 'AT','CA','CC', 'CG', 'CT','GA','GC', 'GG', 'GT','TA','TC', 'TG', 'TT']
    result = {}
    for y in nucs:
        result[y] = {}
        for x in range(-10,10):
            if x==0: continue
            result[y][x] = {}
            for z in range((x+1),11):
                if y == 0: continue
                result[y][x][z] = {}
                for c in combs:
                    result[y][x][z][c] = 0
    return result

def main(chrom_prefix = ""):
    """Generate genome-wide counts of nucleotides at postions relative to CpGs and C/Gs in non-CpG contexts"""
    
    parser = argparse.ArgumentParser(description="Genome-wide counts for CpG/non-CpG single position models.")
    parser.add_argument("-c", "--chromosome", help="Chromsome to get count from", required=True)
    parser.add_argument("-o", "--output", help="Path to output", required=True)
    args = parser.parse_args()
    
    nucs = ['A', 'C', 'G', 'T']
    
    ref = Fasta("/net/snowwhite/home/beckandy/research/BRIDGES_redo/reference_data/human_g1k_v37.fasta")
    chrom = args.chromosome
    out_dir = args.output

    # Load in the chromosome's sequence as a string
    chr_seq = ref['{}{}'.format(chrom_prefix, chrom)]
    chr_string = chr_seq[0:len(chr_seq)].seq

    # Objects to save our results
    cpg_res = gen_empty_dict(["C","G"])
    non_res = gen_empty_dict()

    # Traverse the chromosome!
    start = 0
    end = len(chr_string)

    for i in range(start, end):
        if i % 1000000 == 0:
            print(i)
        nuc = chr_string[i]
        if nuc != "C" and nuc != "G" and nuc != "A" and nuc != "T":
            continue
        cpg_bool = chr_string[i:(i+2)] == "CG" or  chr_string[(i - 1):(i+1)] == "CG"# check if CpG
        # Add counts
        # j and k are positions relative to our nucleotide of interest
        for rp1 in range(-10,10):
            if rp1 == 0: continue
            if i + rp1 < 0: continue
            nuc1 = chr_string[i + rp1]
            if nuc1 not in nucs: continue
            for rp2 in range((rp1+1), 10):
                if i + rp2 >= end: continue
                if rp2 == 0: continue
                nuc2 = chr_string[i + rp2]
                if nuc2 not in nucs: continue
                sub_mot = nuc1 + nuc2
                if cpg_bool:
                    cpg_res[nuc][rp1][rp2][sub_mot] += 1
                else:
                    non_res[nuc][rp1][rp2][sub_mot] += 1

    # Save off results: table per relative position
    for j in range(-10,10):
        if j==0: continue
        for k in range((j+1),11):
            if k == 0: continue
            cpg_c_df = pd.DataFrame.from_dict(cpg_res["C"][j][k], orient='index', columns=['n'])
            non_c_df = pd.DataFrame.from_dict(non_res["C"][j][k], orient='index', columns=['n'])
            cpg_g_df = pd.DataFrame.from_dict(cpg_res["G"][j][k], orient='index', columns=['n'])
            non_g_df = pd.DataFrame.from_dict(non_res["G"][j][k], orient='index', columns=['n'])
            a_df = pd.DataFrame.from_dict(non_res["A"][j][k], orient='index', columns=['n'])
            t_df = pd.DataFrame.from_dict(non_res["T"][j][k], orient='index', columns=['n'])
            #Output file names
            out_file_cpg_c = out_dir + "/cpg_c_chrom{}_p{}_q{}.csv".format(chrom, j, k)
            out_file_non_c = out_dir + "/non_c_chrom{}_p{}_q{}.csv".format(chrom, j, k)
            out_file_cpg_g = out_dir + "/cpg_g_chrom{}_p{}_q{}.csv".format(chrom, j, k)
            out_file_non_g = out_dir + "/non_g_chrom{}_p{}_q{}.csv".format(chrom, j, k)
            out_file_a = out_dir + "/a_chrom{}_p{}_q{}.csv".format(chrom, j, k)
            out_file_t = out_dir + "/t_chrom{}_p{}_q{}.csv".format(chrom, j, k)
            #write to output
            cpg_c_df.to_csv(out_file_cpg_c, index_label="Nucs")
            non_c_df.to_csv(out_file_non_c, index_label="Nucs")
            cpg_g_df.to_csv(out_file_cpg_g, index_label="Nucs")
            non_g_df.to_csv(out_file_non_g, index_label="Nucs")
            a_df.to_csv(out_file_a, index_label="Nucs")
            t_df.to_csv(out_file_t, index_label="Nucs")
    return 0
    
if __name__ == "__main__":
    main()

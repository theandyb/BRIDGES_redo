from pyfaidx import Fasta
import pandas as pd
from Bio.Seq import Seq
import argparse

def gen_empty_dict():
    """Create dict of dicts in which count results will be stored"""
    result = {}
    for y in ["A","T"]:
        result[y] = {}
        for x in range(-10,10):
            if x==0: continue
            result[y][x] = {}
            for z in range((x+1),11):
                if z == 0: continue
                result[y][x][z] = {}
    return result

def main():
    """Generate genome-wide counts of nucleotides at postions relative to CpGs and C/Gs in non-CpG contexts"""
    
    parser = argparse.ArgumentParser(description="Genome-wide counts for CpG/non-CpG single position models.")
    parser.add_argument("-c", "--chromosome", help="Chromsome to get count from", required=True)
    parser.add_argument("-o", "--output", help="Path to output", required=True)
    args = parser.parse_args()
    
    ref = Fasta("/net/snowwhite/home/beckandy/research/germline_full_project/reference_data/human_g1k_v37/human_g1k_v37.fasta")
    chrom = args.chromosome
    out_dir = args.output

    # Load in the chromosome's sequence as a string
    chr_seq = ref['{}'.format(chrom)]
    chr_string = chr_seq[0:len(chr_seq)].seq

    # Objects to save our results
    results = gen_empty_dict()

    # Traverse the chromosome!
    #chr_string = "ACGGCGCGGACCGCGCGGGGGCCCGCGCGCCCAGGG" # Test
    start = 0
    end = len(chr_string)

    for i in range(start, end):
        nuc = chr_string[i]
        if nuc != "A" and nuc != "T":
            continue
        # Grab 21-mer centered at current position
        # Handle chromosome 17 edge case
        # case 1: not enough space upstream
        if i < 10:
            pad_left = -i
            motif_i = i # where in the motif we extract is our "center"
        else:
            pad_left = -10
            motif_i = 10 # in python this is the 11th character
        # Case 2: not enough space downstream
        # This is where I've been frazled: len(str) returns the number of characters, but str[len(str)] is out of bounds 
        # but if slicing str[1:len(str)] works since right is non-inclusive 
        # So we need to check how many characters (up to 10) we can snab up from i
        if i + 11 >= end:
            pad_right = end - 1 - i
        else:
            pad_right = 10
        motif = chr_string[(i+pad_left):(i+pad_right+1)]
        # Now for another pain point: where is our "center" if we have to reverse-complement an unbalanced motif? 
        # Add counts
        # j and k are positions relative to our nucleotide of interest
        for j in range(0,len(motif)-1):
            rp1 = j - motif_i
            if rp1 == 0: continue
            for k in range((j+1), len(motif)):
                rp2 = k - motif_i
                if rp2 == 0: continue
                sub_mot = motif[j] + motif[k]
                if sub_mot in results[nuc][rp1][rp2]:
                    results[nuc][rp1][rp2][sub_mot] += 1
                else:
                    results[nuc][rp1][rp2][sub_mot] = 1
        
    # Save off results: table per relative position
    for j in range(-10,10):
        if j==0: continue
        for k in range((j+1),11):
            if k == 0: continue
            df = pd.DataFrame.from_dict(results["A"][j][k], orient='index', columns=['n'])
            out_file = out_dir + "/a_chrom{}_p{}_q{}.csv".format(chrom, j, k)
            df.to_csv(out_file, index_label="Nucs")
    for j in range(-10,10):
        if j==0: continue
        for k in range((j+1),11):
            if k == 0: continue
            df = pd.DataFrame.from_dict(results["T"][j][k], orient='index', columns=['n'])
            out_file = out_dir + "/t_chrom{}_p{}_q{}.csv".format(chrom, j, k)
            df.to_csv(out_file, index_label="Nucs")
    
if __name__ == "__main__":
    main()
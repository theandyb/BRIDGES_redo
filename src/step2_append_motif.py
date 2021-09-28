from pyfaidx import Fasta
import argparse
import pandas as pd

def get_motif(seqstr, pos, bp = 10):
    return seqstr[(pos - bp):(pos + bp + 1)]

parser = argparse.ArgumentParser(description="Annotate genomic locations with 21-mer motif")
parser.add_argument("-c", "--chrom", help="Which chromosome?", required=True)
parser.add_argument("-s", "--singleton", help="Location of singleton file", required=True)
parser.add_argument("-o", "--output", help="Location of output file", required=True)
args = parser.parse_args()

chrom = args.chrom
singleton_file = args.singleton
output_file = args.output

ref_file = "/net/snowwhite/home/beckandy/research/BRIDGES_redo/reference_data/human_g1k_v37.fasta"
print("Reading reference file...")
fasta_obj = Fasta(ref_file)
seq = fasta_obj["{}".format(chrom)]
seqstr = seq[0:len(seq)].seq
print("Reference read!")

output_list = []

with open(singleton_file) as fp:
    line = fp.readline()
    line = fp.readline()
    cnt = 1
    while line:
        data = line.strip().split("\t") # CHROM, POS, S/D, ALLELE, ID
        chrom = data[0]
        pos = int(data[1])
        alt = data[3]
        subject = data[4]
        motif = get_motif(seqstr, pos-1)
        cat = "{}>{}".format(seqstr[(pos-1)], alt)
        entry = {
            'chrom' : chrom,
            'pos' : pos,
            'motif' : motif,
            'cat': cat,
            'alt': alt,
            'subject': subject
        }
        output_list.append(entry)
        if cnt % 10000 == 0:
            pd.DataFrame(output_list).to_csv(output_file, index = None, header=False, mode='a')
            output_list = []
        line = fp.readline()

if output_list:
    pd.DataFrame(output_list).to_csv(output_file, index = None, header=False, mode='a')

from Bio.Seq import Seq
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description="Reverse complement control files")
parser.add_argument("-i", "--input", help="Location of input file (without extension)", required=True)
parser.add_argument("-n", "--nuc", help="Which REF do we not RC?", required=True)
args = parser.parse_args()

input_file = "{}.csv".format(args.input)
output_file = "{}_rc.csv".format(args.input)
check_nuc = args.nuc

df = pd.read_csv(input_file, names = ['Chrom', 'Pos', 'motif', 'subtype', 'REF', 'window', 'distance'])
df['motif2'] = df.apply(lambda x: x['motif'] if x['motif'][10] == check_nuc else str(Seq(x['motif']).reverse_complement()), axis = 1)

df.to_csv(output_file, index=False)

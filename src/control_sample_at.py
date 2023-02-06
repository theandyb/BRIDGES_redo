from pyfaidx import Fasta
import regex as re
import random
import pandas as pd
import argparse
from Bio.Seq import Seq

def full_motif(motif, ref):
    motif_seq = Seq(motif)
    motif_rc = motif_seq.reverse_complement().__str__()
    if ref in ['A','C']:
        final = "{}({})".format(motif, motif_rc)
    else:
        final = "{}({})".format(motif_rc, motif)
    return final

def main(ref_prefix = ""):
    parser = argparse.ArgumentParser(description="Sample control distribution.")
    parser.add_argument("-s", "--singleton", help="Location of singleton file", required=True)
    parser.add_argument("-f", "--fasta", help="FASTA file to grab sequence from", required=True)
    parser.add_argument("-o", "--output", help="Path to output", required=True)
    parser.add_argument("-n", "--nSample", help="Number of controls per singleton", type=int, default = 1)
    parser.add_argument("chrom", help="Chromosome we are sampling from")
    args = parser.parse_args()

    ref_file = args.fasta 
    singleton_file = args.singleton
    chrom = args.chrom
    nSample = args.nSample
    output_list = []
    # Addition: 10-Nov-22
    out_file = args.output
    min_file = out_file + ".min"
    max_file = out_file + ".max"
    min_list = []
    max_list = []
    # Create fasta object
    print("Reading fasta...")
    fasta_obj = Fasta(ref_file)
    seq = fasta_obj["{}{}".format(ref_prefix, chrom)] # hg37 does not have chr prefix, but hg38 does
    seqstr = seq[0:len(seq)].seq
    print("FASTA read!")
    # Iterate over singletons file
    print("Sampling control observations for singletons...")
    counter = 1
    with open(singleton_file) as fp:
        line = fp.readline()
        while line:
            # columns are: 0:chrom, 1:pos, 2:motif, 3:subtype, 4:alt, 5:id, 6:REF, 7:motif2, 8:subtype2
            content = line.strip().split(",")
            pos = int(content[1])
            ref = content[6]
            motif = content[7][:21]
            cat = content[8]
            if(not cat.startswith("A")):
                line = fp.readline()
                continue
            new_entry = sample_control(chrom, pos, ref, cat, seqstr, nSample)
            # New: get control nearest and further from singleton
            
            output_list.extend(new_entry)
            # Find control with minimum and maximum distance
            min_dist = min([t.get("distance") for t in new_entry])
            max_dist = max([t.get("distance") for t in new_entry])
            min_entry = [t for t in new_entry if t.get('distance') == min_dist]
            max_entry = [t for t in new_entry if t.get('distance') == max_dist ]
            min_list.extend(min_entry)
            max_list.extend(max_entry)
            line = fp.readline()
            counter += 1
            if counter % 10000 == 0:
                print(counter)
                df = pd.DataFrame(output_list)
                df.to_csv(out_file, index = None, header=False, mode='a')
                output_list = []
                # New: add min and max files
                pd.DataFrame(min_list).to_csv(min_file, index = None, header=False, mode='a')
                min_list = []
                pd.DataFrame(max_list).to_csv(max_file, index = None, header=False, mode='a')
                max_list = []
    print("Done sampling...")
    if output_list:
        pd.DataFrame(output_list).to_csv(out_file, index = None, header=False, mode='a')
        pd.DataFrame(min_list).to_csv(min_file, index = None, header=False, mode='a')
        pd.DataFrame(max_list).to_csv(max_file, index = None, header=False, mode='a')
    print("Done!")

def sample_control(chrom, pos, ref, cat, seq, nSample, window=150, bp=10):
  sites = []
  newlist = []
  while(len(sites) < nSample + 5):
    #lseg_lb = max((pos-1-window-bp), 0)
    #useg_ub = upBound = min(len(seq), pos + window + bp)
    #subseq = seq[lseg_lb:useg_ub]
    subseq = seq[(pos - 1 - window):(pos+window)]
    #subseq = re.sub(r'^N+', '', subseq)
    #search_pat = "(?=([ATCG]{%d}%s[ATCG]{%d}))" % (bp, ref, bp)
    #sites = [m.start() for m in re.finditer(search_pat, subseq)]
    sites = [m.start() for m in re.finditer(ref, subseq)]
    sites = [s for s in sites if (s > bp+window+1 or s < window-bp-1)]
    window += 50 #expand window in edge case where mut_site is only ref_allele in window
  window -= 50
  while len(newlist) < nSample:
    ix = random.choice(sites)
    #chrom_ix = ix - window + pos + bp
    chrom_ix = ix - window + pos
    try:
      newSeq = seq[(chrom_ix - bp - 1):(chrom_ix+bp)].upper()
      distance = abs(ix - window)
      motif2 = full_motif(newSeq, ref)
      entry = {
        'chrom' : chrom,
        'pos' : pos,
        'motif' : newSeq,
        'cat': cat,
        'ref': ref,
        'window': window,
        'distance': distance,
        'motif2':motif2,
        'spos': chrom_ix
      }
      newlist.append(entry)
    except ValueError:
      print("Edge Case: position " + str(pos))
    sites.remove(ix)
  return newlist



if __name__ == "__main__":
    #random.seed( 8675 ) # threeeeee ohhhhh niiiiii-eee-iiiiine
    random.seed( 1776 ) # round two
    main()

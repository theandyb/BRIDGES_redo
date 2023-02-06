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
    min_list = []
    max_list = []
    min_file = out_file + ".min"
    max_file = out_file + ".max"
    # Create fasta object
    fasta_obj = Fasta(ref_file)
    seq = fasta_obj["{}{}".format(ref_prefix, chrom)] # hg37 has chromosomes indexed without the text "chr", whereas hg38 does 
    seqstr = seq[0:len(seq)].seq
    print("FASTA read!")
    # Iterate over singletons file
    print("Sampling control observations for singletons...")
    counter = 1
    bad_sites = 0
    with open(singleton_file) as fp:
        line = fp.readline()
        while line:
            content = line.strip().split(",")
            # columns are: 0:chrom, 1:pos, 2:motif, 3:subtype, 4:alt, 5:id, 6:REF, 7:motir2, 8:subtype2
            pos = int(content[1])
            ref = content[6]
            motif = content[7][:21]
            cat = content[8]
            
            if(cat.startswith("A")):
                line = fp.readline()
                continue
            elif(cat.startswith("cpg")):
                cpg_bool = True
            else:
                cpg_bool = False
            
            new_entry = sample_control(chrom, pos, ref, cat, nSample, cpg_bool, seqstr)
            if new_entry == 0:
              bad_sites += 1
            else:
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
            if counter % 10000 == 0 and output_list:
                print(counter)
                pd.DataFrame(output_list).to_csv(out_file, index = None, header=False, mode='a')
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
    print("Number singletons without matched controls: " + str(bad_sites))
    print("Total singletons: " + str(counter))

def sample_control(chrom, pos, ref, cat, nSample, cpg_bool, seqstr, window=150, bp=10):
  sites = []
  newlist = []
  # NEW CODE FOR SAMPLING CpG SITES (OR NON)
  if(cpg_bool):
    search_str = "CG"
  else:
    if ref == "C":
      search_str = "C[ACT]"
    else:
      search_str = "[AGT]G"
  while len(sites) < nSample + 1:
    subseq = seqstr[(pos - 1 - window):(pos+window)]
    #subseq = re.sub(r'^N+', '', subseq) # Trim Ns at beginning or end of sequence
    #search_pat = "(?=([ATCG]{%d}%s[ATCG]{%d}))" % (bp, search_str, (bp - 1))
    #sites = [m.start() for m in re.finditer(search_pat, subseq)] 
    sites = [m.start() for m in re.finditer(search_str, subseq, overlapped=True)] # These two lines were changed for CpG/non-CpG sampling
    sites = [s for s in sites if (s > bp+window+1 or s < window-bp-1)]# Identify all possible motifs to sample from
    window += 100
  window -= 100
  while len(newlist) < nSample:
    if(len(sites)==0):
      print("Bad pos: {}".format(pos))
    ix = random.choice(sites)
    if search_str == "[AGT]G":
      ix = ix + 1
    #chrom_ix = ix - window + pos + bp
    chrom_ix = ix - window + pos
    try:
      newSeq = seqstr[(chrom_ix - bp - 1):(chrom_ix+bp)].upper()
      motif2 = full_motif(newSeq, newSeq[bp])
      distance = abs(ix - window)
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
    except:
      print("Sumthin' bad happened maaaaaan")
    finally:
      if search_str == "[AGT]G":
        ix = ix - 1
      sites.remove(ix)
  return newlist



if __name__ == "__main__":
    #random.seed( 8675 ) # threeeeee ohhhhh niiiiii-eee-iiiiine
    random.seed( 1776 ) # round two
    main()

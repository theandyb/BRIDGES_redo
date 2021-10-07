from pyfaidx import Fasta
import regex as re
import random
import pandas as pd
import argparse


def main():
    parser = argparse.ArgumentParser(description="Sample control distribution.")
    parser.add_argument("-s", "--singleton", help="Location of singleton file", required=True)
    parser.add_argument("-f", "--fasta", help="FASTA file to grab sequence from", required=True)
    parser.add_argument("-o", "--output", help="Path to output", required=True)
    parser.add_argument("-n", "--nSample", help="Number of controls per singleton", type=int, default = 1)
    parser.add_argument("chrom", help="Chromosome we are sampling from")
    args = parser.parse_args()

    ref_file = args.fasta #"reference_data/human_g1k_v37/human_g1k_v37.fasta"
    singleton_file = args.singleton #"/net/snowwhite/home/beckandy/research/smaug-redux/summaries/singletons.full.summary"
    chrom = args.chrom
    nSample = args.nSample
    output_list = []
    # Create fasta object
    fasta_obj = Fasta(ref_file)
    seq = fasta_obj["{}".format(chrom)] # hg37 has chromosomes indexed without the text "chr"
    seqstr = seq[0:len(seq)].seq
    print("FASTA read!")
    # Iterate over singletons file
    print("Sampling control observations for singletons...")
    counter = 1
    with open(singleton_file) as fp:
        next(fp)
        line = fp.readline() # first line is a header
        print(line)
        while line:
            content = line.strip().split(",")
            # columns are: 0:chrom, 1:pos, 2:motif, 3:subtype, 4:alt, 5:id, 6:REF, 7:ALT, 8:subtype2, 9:motif2
            pos = int(content[1])
            ref = content[6]
            motif = content[9][:21]
            cat = content[8]
            
            if(cat.startswith("A")):
                line = fp.readline()
                continue
            elif(cat.startswith("cpg")):
                cpg_bool = True
            else:
                cpg_bool = False
            
            output_list.extend(sample_control(chrom, pos, ref, cat, seq, nSample, cpg_bool, seqstr))
            line = fp.readline()
            counter += 1
            if counter % 10000 == 0:
                print(counter)
                df = pd.DataFrame(output_list)
                df.to_csv(args.output, index = None, header=False, mode='a')
                output_list = []
    print("Done sampling...")
    if output_list:
        pd.DataFrame(output_list).to_csv(args.output, index = None, header=False, mode='a')
    print("Done!")

def sample_control(chrom, pos, ref, cat, seq, nSample, cpg_bool, seqstr, window=150, bp=10):
    sites1 = []
    sites2 = []
    newlist = []
    # NEW CODE FOR SAMPLING CpG SITES (OR NON)
    if(cpg_bool):
        search_str = "CG"
    else:
        if ref == "C":
            search_str = "C[ACT]"
        else:
            search_str = "[AGT]G"
        
    while(len(sites1) + len(sites2) < nSample):
        # Break 300bp into two halves, separated by 21-mer centered at POS
        lseg_lb = max((pos-1-window-bp), 0)
        lseg_ub = pos - bp - 1
        useg_lb = pos + bp
        useg_ub = upBound = min(len(seq), pos + window + bp)
        subseq1 = seqstr[lseg_lb:lseg_ub]
        subseq2 = seqstr[useg_lb:useg_ub]
        subseq1 = re.sub(r'^N+', '', subseq1) # Trim Ns at beginning or end of sequence
        subseq2 = re.sub(r'N+$', '', subseq2)
        sites1 = [m.start() for m in re.finditer(search_str, subseq1, overlapped=True)] # These two lines were changed for CpG/non-CpG sampling
        sites2 = [m.start() for m in re.finditer(search_str, subseq2, overlapped=True)]
        sites1 = [s for s in sites1 if (s >= bp+1 and s < (len(subseq1)-bp-1))] # Identify all possible motifs to sample from
        sites2 = [s for s in sites2 if (s >= bp+1 and s < (len(subseq2)-bp-1))]
        window += 50 #expand window in edge case where mut_site is only ref_allele in window
    window -= 50
    while len(newlist) < nSample:
        flip = random.randint(0, 1)
        if ((flip == 0 and len(sites1)>0) or (len(sites2)==0)):
            subseq = subseq1
            sites = sites1
            c_direction = -1
        else:
            subseq = subseq2
            sites = sites2
            c_direction = 1
        if(len(sites)==0):
            print("Bad pos: {}".format(pos))
        ix = random.choice(sites)
        if ref == "C":
            newSeq = subseq[(ix - bp):(ix+bp+1)]
        else:
            newSeq = subseq[(ix + 1 - bp):(ix+bp+2)]
        while not re.search("[NATCG]{21}", newSeq): # THIS CHANGED!
            #print(pos)
            sites.remove(ix)
            ix = random.choice(sites)
            if ref == "C":
                newSeq = subseq[(ix - bp):(ix+bp+1)]
            else: 
                newSeq = subseq[(ix + 1 - bp):(ix+bp+2)]
        if c_direction == -1:
            distance = window + bp - ix
        else:
            distance = window + ix
        entry = {
            'chrom' : chrom,
            'pos' : pos,
            'motif' : newSeq,
            'cat': cat,
            'ref': ref,
            'window': window,
            'distance': distance
        }
        newlist.append(entry)
        if ((flip == 0 and len(sites1)>0) or (len(sites2)==0)):
            sites1.remove(ix)
        else:
            sites2.remove(ix)
    return newlist



if __name__ == "__main__":
    #random.seed( 8675 ) # threeeeee ohhhhh niiiiii-eee-iiiiine
    random.seed( 1776 ) # round two
    main()

#!/usr/bin/env python
import collections
import itertools
import argparse
import os

#usage: python ~/github/Bacterial_rRNA_db/src/Bacterial_rRNA_db/GenusSeqDiffStats.py --indir NCBI_accession_16SrRNA --outdir NCBI_accession_16SrRNA_compare


__author__ = "Zhikun Wu"
__email__ = "598466208@qq.com"
__date__ = "2018.06.11"

def read_muscle_clw_file(align_file):
    """
    align_file example:
    MUSCLE (3.8) multiple sequence alignment


    Ureaplasma_diversum         ATTTTTAAGAGTTTGATCCTGGCTCAGGATTAACGCTGGCGGCATGCCTAATACATGCAA
    Ureaplasma_urealyticum      ATTTTAAAGAGTTTGATCCTGGCTCAGGATTAACGCTGGCGGCATGCCTAATACATGCAA
    Ureaplasma_parvum           ATTTTAAAGAGTTTGATCCTGGCTCAGGATTAACGCTGGCGGCATGCCTAATACATGCAA
                                ***** ******************************************************


    """
    SpeciesSeq = collections.defaultdict()
    in_h = open(align_file, "r")
    head = in_h.readline().strip()
    if head.startswith(">"):
        SpeciesSeq = None
    else:
        for line in in_h:
            line = line.strip()
            if line:
                if not line.startswith("*"):
                    lines = line.split()
                    species, seq = lines
                    if species not in SpeciesSeq:
                        SpeciesSeq[species] = seq
                    else:
                        SpeciesSeq[species] += seq
    in_h.close()
    return SpeciesSeq

def compare_two_seq(seq1, seq2):
    """
    seq1: AGGATTAGATACCCTAGTAGTCCACACCGTAAACGATCATCATTAAATGTCGGCTCGCTT
    seq2: AGGATTAGATACCCTAGTAGTCCACACCGTAAACGATCATCATTAAATGTCGGCTCGA--
          *********************************************************
          
    """
    count = 0
    if len(seq1) == len(seq2):
        for i in range(len(seq1)):
            if seq1[i] != seq2[i]:
                count += 1
    else:
        print("Please check whether the lengths of seq1 and seq2 are equal:\n%s\n%s\n" % (seq1, seq2))
    return count

def compare_seq_of_dict(align_file, out_file):
    SpeciesSeq = read_muscle_clw_file(align_file)
    ### no output file if SpeciesSeq is None
    if SpeciesSeq != None:
        out_h = open(out_file, "w")
        Species = sorted(list(SpeciesSeq.keys()))
        combinations = itertools.combinations(Species, 2)
        for c in combinations:
            s1, s2 = c
            seq1 = SpeciesSeq[s1]
            seq2 = SpeciesSeq[s2]
            count = compare_two_seq(seq1, seq2)
            out_h.write("%s\t%s\t%s\n" % (s1, s2, count))
        out_h.close()

def add_dash(indir):
    if not indir.endswith("/"):
        indir += "/"
    return indir

def out_to_dir(indir, outdir):
    files = os.listdir(indir)
    indir = add_dash(indir)
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    outdir = add_dash(outdir)
    for f in files:
        in_file = indir + f
        if in_file.endswith(".clw"):
            out_file = outdir + f
            out_file = out_file.rstrip(".clw")
            compare_seq_of_dict(in_file, out_file)


def main():
    parser = argparse.ArgumentParser(description="Get the different number of two sequences based on the muscle aligned file.")
    parser.add_argument("-i", "--indir", help="The input file.")
    parser.add_argument("-o", "--outdir", help="The output file.")
    args = parser.parse_args()
    out_to_dir(args.indir, args.outdir)
    # compare_seq_of_dict(args.indir, args.outdir)

if __name__ == "__main__":
    main()

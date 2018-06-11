#!/usr/bin/env python
from tinyfasta import FastaParser
import collections
import argparse
import os

#usage: python SameGenusSeq.py --input NCBI_accession_16SrRNA.fasta --outdir NCBI_accession_16SrRNA

__author__ = "Zhikun Wu"
__email__ = "598466208@qq.com"
__date__ = "2018.06.11"

def same_genus_sequence(seq_file):
    """
    seq_file example:
    >GCF_000010525.1,NC_009937.1,680506,681988      Bacteria;Proteobacteria;Alphaproteobacteria;Rhizobiales;Xanthobacteraceae;Azorhizobium;Azorhizobium caulinodans
    CAACTTGAGAGTTTGATCCTGGCTCAGAGCGAACGCTGGCGGCAGGCCTAACACATGCAAGTCGAGCGCCCAGCAATGGGAGCGGCAGACGGGTGAGTAACGCGTGGGGATGTGCCCAATGGTGCGGAATAACCCAGGGAAA
    """
    GenusSeq = collections.defaultdict(list)
    for record in FastaParser(seq_file):
        desc = str(record.description)
        seq = str(record.sequence)
        descs = desc.split("\t")
        # accession = descs[0]
        genus = '_'.join(descs[-1].split(";")[:-1])
        accession = '>' + '_'.join(descs[-1].split(";")[-1].split())
        GenusSeq[genus].append([accession, seq])
    return GenusSeq


def out_seq(records, name, outdir):
    if not os.path.exists(outdir) or not os.path.isdir(outdir):
        os.makedirs(outdir)
    if not outdir.endswith("/"):
        outdir += "/"
    file = outdir + name
    out_h = open(file, "w")
    for r in records:
        out_h.write("%s\n" % "\n".join(r))
    out_h.close()
    ### mutiple alignment
    if not file.endswith("None"):
        cmd = "muscle -in %s -out %s -clw" % (file, file + ".clw")
        os.system(cmd)

def out_genus_seq(seq_file, out_dir):
    GenusSeq = same_genus_sequence(seq_file)
    for g in GenusSeq:
        seqs = GenusSeq[g]
        out_seq(seqs, g, out_dir)

def main():
    parser = argparse.ArgumentParser(description="Out the sequence into the file for the same genus.")
    parser.add_argument("-i", "--input", help="The input file.")
    parser.add_argument("-o", "--outdir", help="The output directory.")
    args = parser.parse_args()
    out_genus_seq(args.input, args.outdir)

if __name__ == "__main__":
    main()

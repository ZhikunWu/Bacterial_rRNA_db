#!/usr/bin/env python
from tinyfasta import FastaParser
import collections
import argparse

#usage: python ~/github/Bacterial_rRNA_db/src/Bacterial_rRNA_db/SpeciesNum.py --fasta ../database/NCBI_accession_16SrRNA.fasta --mapping fastqjoin.join_record.txt --out fastqjoin.join_taxonomy.txt

__author__ = "Zhikun Wu"
__email__ = "598466208@qq.com"
__date__ = "2018.06.12"

def get_accession_taxonomy(ref_seq):
    """
    ref_seq example:
    >GCF_000010525.1,NC_009937.1,680506,681988      Bacteria;Proteobacteria;Alphaproteobacteria;Rhizobiales;Xanthobacteraceae;Azorhizobium;Azorhizobium caulinodans
    CAACTTGAGAGTTTGATCCTGGCTCAGAGCGAACGCTGGCGGCAGGCCTAACACATGCAAGTCGAGCGCCCAGCAATGGGAGC
    """
    AccTaxo = {}
    for record in FastaParser(ref_seq):
        desc = str(record.description)
        desc = desc.lstrip(">")
        descs = desc.split("\t")
        accession, taxonomy = descs
        AccTaxo[accession] = taxonomy
    return AccTaxo

def parse_mapping_taxonomy(ref_seq, mapping_record, out_file):
    """
    mapping_record example:
    E00602:91:HK2HGCCXY:8:1101:25175:2663   GCF_000154805.1,NZ_DS562845.1,220941,222348
    E00602:91:HK2HGCCXY:8:1101:17868:2698   GCF_000157915.1,NZ_EQ973630.1,5061,3531
    """
    TaxoCount = collections.Counter()
    AccTaxo = get_accession_taxonomy(ref_seq)
    in_h = open(mapping_record, "r")
    for line in in_h:
        lines = line.strip().split("\t")
        query, subject = lines
        if subject in AccTaxo:
            taxo = AccTaxo[subject]
            TaxoCount[taxo] += 1
        else:
            print("Please check whether the accession %s have the corresponding taxonomy in the file %s." % (subject, ref_seq))
    in_h.close()
    ### outpu the taxonomy and the counts
    out_h = open(out_file, "w")
    out_h.write("Taxonomy\tCount\n")
    taxos = sorted(list(TaxoCount.keys()))
    for t in taxos:
        count = TaxoCount[t]
        out_h.write("%s\t%d\n" % (t, count))
    out_h.close()


def main():
    parser = argparse.ArgumentParser(description="Get the taxonomy number based on the mapping file and fasta file containing taxonomy.")
    parser.add_argument("-f", "--fasta", help="The input fasta file containing the accession and taxonomy in the description.")
    parser.add_argument("-m", "--mapping", help="The input mapping record file.")
    parser.add_argument("-o", "--out", help="The output file.")
    args = parser.parse_args()
    parse_mapping_taxonomy(args.fasta, args.mapping, args.out)

if __name__ == "__main__":
    main()





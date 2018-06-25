#!/usr/bin/env python
import argparse
from tinyfasta import FastaParser

#usage: python ~/github/Bacterial_rRNA_db/src/Bacterial_rRNA_db/GreenGeneSpecies.py --taxonomy taxonomy/99_otu_taxonomy.txt --sequence rep_set/99_otus.fasta --species 99_otu_species.txt --genus 99_otu_genus.txt 

__author__ = "Zhikun Wu"
__email__ = "598466208@qq.com"
__date__ = "2018.06.25"

"""
Get the taxonomy with genus and species leves and corresponding sequenece based on the taxonomy id and sequence in Greengene database.
"""


def trans_taxos(taxo):
    """
    argv:
        k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Staphylococcaceae; g__Staphylococcus; s__
    return:
        Bacteria;Firmicutes;Bacilli;Bacillales;Staphylococcaceae;Staphylococcus
    """
    taxos = taxo.split("; ")
    new_taxos = [t.split("__")[-1].strip("[").strip("]") for t in taxos]
    new_taxo = ";".join(new_taxos)
    new_taxo = new_taxo.rstrip(";")
    return new_taxo

def get_taxonomy_seq_id(otu_taxo):
    seqIDSpecies = {}
    seqIDGenus = {}
    taxo_h = open(otu_taxo, "r")
    for line in taxo_h:
        lines = line.strip().split("\t")
        seqID, Taxonomy = lines
        Taxos = Taxonomy.split("; ")
        genus = Taxos[5]
        species = Taxos[6]
        if not genus.endswith("g__"):
            if not species.endswith("s__"):
                seqIDSpecies[seqID] = trans_taxos(Taxonomy)
            else:
                seqIDGenus[seqID] = trans_taxos(Taxonomy)
    taxo_h.close()
    return seqIDGenus, seqIDSpecies

def get_seq_taxo(otu_taxo, otu_seq, spe_out, gen_out):
    """
    argv:
        otu_taxo:
        228054  k__Bacteria; p__Cyanobacteria; c__Synechococcophycideae; o__Synechococcales; f__Synechococcaceae; g__Synechococcus; s__
        228057  k__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rickettsiales; f__Pelagibacteraceae; g__; s__
        
        otu_seq:
        >1111886
        AACGAACGCTGGCGGCATGCCTAACACATGCAAGTCGAACGAGACCTTCGGGTCTAGTGGCGCACGGGTGCGTAACGCGTGGGAATCTGCCCTTGGGTACGGAATAACAGTTAGAAATGACTGCTAATACCGTATA
        >1111883
        GCTGGCGGCGTGCCTAACACATGTAAGTCGAACGGGACTGGGGGCAACTCCAGTTCAGTGGCAGACGGGTGCGTAACACGTGAGCAACTTGTCCGACGGCGGGGGATAGCCGGCCCAACGGCCGGGTAATACCGCG

    return: 
        >Bacteria;Proteobacteria;Gammaproteobacteria;Xanthomonadales;Xanthomonadaceae;Pseudoxanthomonas;Pseudoxanthomonas mexicana
        AGAGTTTGATCATGGCTCAGAGTGAACGCTGGCGGTAGGCCTAACACATGCAAGTCGAACGGCAGCACAGGAAAGCTTGCTCTCTGGGTGGCGAGTGGCGGACGGGTGAGGAATACATCGGAAT
    """
    seqIDGenus, seqIDSpecies = get_taxonomy_seq_id(otu_taxo)
    spe_h = open(spe_out, "w")
    gen_h = open(gen_out, "w")
    for record in FastaParser(otu_seq):
        desc = str(record.description).lstrip(">")
        seq = str(record.sequence)
        if desc in seqIDSpecies:
            taxo = seqIDSpecies[desc]
            taxos = taxo.split(";")
            new_taxo = ";".join(taxos[:6]) + ";" + taxos[5] + " " + taxos[6]
            spe_h.write(">%s\n%s\n" % (new_taxo, seq))
        elif desc in seqIDGenus:
            taxo = seqIDGenus[desc]
            gen_h.write(">%s\n%s\n" % (taxo, seq))

def main():
    parser = argparse.ArgumentParser(description = "Get the taxonomy with genus and species leves and corresponding sequenece based on the taxonomy id and sequence in Greengene database.")
    parser.add_argument("-t", "--taxonomy", help="The input taxonomy file of Greengene database.")
    parser.add_argument("-s", "--sequence", help="The sequence file of Greengene database.")
    parser.add_argument("-sp", "--species", help="The output sequence with species level.")
    parser.add_argument("-ge", "--genus", help="The output sequence with genus level.")
    args = parser.parse_args()
    get_seq_taxo(args.taxonomy, args.sequence, args.species, args.genus)

if __name__ == "__main__":
    main()


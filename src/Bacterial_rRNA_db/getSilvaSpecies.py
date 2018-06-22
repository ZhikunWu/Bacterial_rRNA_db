#!/usr/bin/env python
from tinyfasta import FastaParser
import argparse
import re
import sqlite3

#usage: python ~/github/Bacterial_rRNA_db/src/Bacterial_rRNA_db/getSilvaSpecies.py --input SILVA_132_SSUParc_tax_silva_DNA.fasta_example  --out out --database out.db  --lenThreshold 1200

__author__ = "Zhikun Wu"
__email__ = "5984662082@qq.com"
__date__ = "2018.05.28"

def is_exist_str(st, list_str):
    """
    a = ["Lachnospiraceae", "Lachnoclostridium 5"]
    b = is_exist_str("\d+", a)
    print(b)
    True
    """
    exist = False
    for s in list_str:
        match = re.findall(st, s)
        if match:
            exist = True
            break
    return exist

def RNA2DNA(seq):
    """
    Convert the RNA sequence to DNA sequence, U -> T, u -> T.
    argv:
    'GAGUUUGAUCAUGGCUCA'

    return:
    'GAGTTTGATCATGGCTCA'
    """
    seq = seq.upper()
    seq = seq.replace("U", "T")
    return seq


def delete_given_str(st, list_str):
    """
    a = "Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;Prevotellaceae;Prevotella 1;Prevotella bryantii"
    b= is_exist_str("\d+", a.split(";"))
    print(b)
    Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;Prevotellaceae;Prevotella;Prevotella bryantii
    """
    new_strs = []
    for s in list_str:
        match = re.findall('(.*)\s+%s' % st, s)
        if match:
            new_strs.append(match[0])
        else:
            new_strs.append(s)
    return new_strs


def silva_species(silva_fa, lenThreshold, out_file):
    """
    silva_fa example:

    >AY560638.1.507 Bacteria;Firmicutes;Bacilli;Bacillales;Bacillaceae;Bacillus;Nostoc sp. 'Peltigera collina cyanobiont'
    >AY316687.1.521 Bacteria;Proteobacteria;Alphaproteobacteria;Rhizobiales;Rhizobiaceae;Ensifer;Reichenowia ornatae
    >AY584513.1.2387 Bacteria;Cyanobacteria;Oxyphotobacteria;Nostocales;Nostocaceae;Nodularia PCC-9350;Nodularia harveyana CCAP 1452/1
    >AY691371.1.390 Bacteria;Proteobacteria;Gammaproteobacteria;Betaproteobacteriales;Burkholderiaceae;Burkholderia-Caballeronia-Paraburkholderia;Burkholderia sp. tpig5.3
    >U87828.1.1484  Bacteria;Firmicutes;Bacilli;Lactobacillales;Streptococcaceae;Streptococcus;unidentified eubacterium
    >U87831.1.1433  Bacteria;Proteobacteria;Alphaproteobacteria;Rhizobiales;Rhizobiaceae;Bartonella;unidentified proteobacterium
    >JQ811215.1.1434        Bacteria;Tenericutes;Mollicutes;Mollicutes Incertae Sedis;Unknown Family;Candidatus Phytoplasma;Tomato big    
    """
    lenThreshold = int(lenThreshold)
    GenusLevel = set()
    out_h = open(out_file, "w")
    for record in FastaParser(silva_fa):
        desc = str(record.description)
        seq = str(record.sequence)
        seqLen = len(seq)
        accession = desc.split()[0]
        descs = ' '.join(desc.split()[1:]).split(";")
        descLen = len(descs)
        kindom = descs[0]
        ### select Bacterial with seven levels.   ## Archaea
        if kindom == "Bacteria" and descLen == 7:
            genus_length = [len(s.split()) for s in descs[:-1]]
            genus_length_sum = sum(genus_length)
            non_normal = is_exist_str("-", descs[:-1])
            if genus_length_sum == 6 and not non_normal:
            # if not non_normal:
                genus_species = descs[-1].split()
                if len(genus_species) >= 2:
                    genus, spe_string = genus_species[:2]
                    spe_match1 = re.findall("sp\.$", spe_string)
                    # spe_match2 = re.findall("\d+", spe_string)
                    spe_match2 = re.findall("^bacterium", genus)
                    spe_match3 = re.findall("uncultured", genus)
                    spe_match4 = re.findall("unidentified", genus)
                    if genus == descs[5] and not (spe_match2 or spe_match3 or spe_match4):
                        seven_level = descs[:6] + ['%s %s' % (descs[5], spe_string)]
                        seven_level = delete_given_str("\d+", seven_level)
                        ### filte out the short length for sequence
                        GenusLevel.add(tuple(seven_level))
                        if seqLen >= lenThreshold:
                            seq = RNA2DNA(seq)
                            if spe_match1:
                                print(genus_species)
                                spe_string = '_'.join(genus_species[1:]).replace('sp._', '')
                                out_h.write("%s\t%s;%s %s\n%s\n" % (accession, ";".join(descs[:6]), genus, spe_string, seq))
                            else:
                                out_h.write("%s\t%s;%s %s\n%s\n" % (accession, ";".join(descs[:6]), genus, spe_string, seq))
                    else:
                        print(desc)
            # else:
            #     print(desc)
    return GenusLevel

def silva_species_to_db(silva_fa, lenThreshold, out_file, out_db):
    """
    out_db example:

    sqlite> select * from SpeciesTaxonomy limit 5;
    Kindom|Phylum|Class|Orders|Family|Genus|Species
    Bacteria|Firmicutes|Bacilli|Bacillales|Bacillaceae|Bacillus|circulans
    Bacteria|Proteobacteria|Gammaproteobacteria|Pseudomonadales|Moraxellaceae|Acinetobacter|nectaris
    Bacteria|Actinobacteria|Actinobacteria|Streptomycetales|Streptomycetaceae|Streptomyces|costaricanus
    Bacteria|Actinobacteria|Actinobacteria|Actinomycetales|Actinomycetaceae|Varibaculum|anthropi
    Bacteria|Proteobacteria|Alphaproteobacteria|Rhodospirillales|Magnetospirillaceae|Magnetospirillum|gryphiswaldense
    """
    GenusLevel = silva_species(silva_fa, lenThreshold, out_file)
    conn = sqlite3.connect(out_db)
    cursor= conn.cursor()
    cursor.execute("""
        CREATE TABLE SpeciesTaxonomy
        (Kindom text,
        Phylum text,
        Class text,
        Orders text,
        Family text,
        Genus text,
        Species text);
    """)
    for s in GenusLevel:
        # print(s)
        cursor.execute("""
            INSERT INTO SpeciesTaxonomy
            (Kindom, Phylum,
            Class, Orders,
            Family, Genus,
            Species)
            VALUES
            (?, ?, ?, ?, ?, ?, ?);
        """, (s))
    cursor.execute("""
        CREATE INDEX species2taxo
        on SpeciesTaxonomy (Species);
    """)
    cursor.close()
    conn.commit()
    conn.close()




def main():
    parser = argparse.ArgumentParser(description="Get the species with different level names for acterial.")
    parser.add_argument("-i", "--input", help="The input file with SILVA taxonomy name and sequence.")
    parser.add_argument("-t", "--lenThreshold", help="The threshold for sequence length.")
    parser.add_argument("-o", "--out", help="The output file.")
    parser.add_argument("-d", "--database", help="The output database.")
    args = parser.parse_args()
    silva_species_to_db(args.input, args.lenThreshold, args.out, args.database)

if __name__ == "__main__":
    main()



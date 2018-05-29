#/usr/bin/env python
import sqlite3
import argparse
import re
import os

#usage: python ~/github/Bacterial_rRNA_db/src/Bacterial_rRNA_db/assembly2db.py --input assembly_summary.txt --database assembly_summary.db

__author__ = "Zhikun Wu"
__email__ = "598466208@qq.com"
__date__ = "2018.05.28"

def parse_assembly_summary(assembly_summary):
    """
    assembly_summary example:
    GCF_000007725.1 PRJNA224116 SAMN02604289        representative genome   224915  9   Buchnera aphidicola str. Bp (Baizongia pistaciae)   strain=Bp (Baizongia pistaciae)     latest  Complete Genome Major   Full    2003/01/29  ASM772v1    Valencia Univ.  GCA_000007725.1 identical   ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/007/725/GCF_000007725.1_ASM772v1
    """
    in_h = open(assembly_summary, "r")
    Records = []
    for line in in_h:
        line = line.strip()
        if line.startswith("#") or line == "":
            continue
        else:
            lines =line.split("\t")
            accession, bioproject, biosample, wgs_master, refseq_category, taxid, species_taxid, organism_name, infraspecific_name, isolate, version_status, assembly_level, release_type, genome_rep, seq_rel_date, asm_name, submitter, gbrs_paired_asm, paired_asm_comp, ftp_path = lines[:20]
            base_name = os.path.basename(ftp_path)
            genome_name = base_name + "_genomic.fna"
            names = organism_name.split()
            ### species name contain genus and species, if it just contain genus name, the specie name may be "sp.", "genomosp.", "genosp.".
            spe_string = names[1]
            spe_match1 = re.findall("sp\.$", spe_string)
            spe_match2 = re.findall("\d+", spe_string)
            # spe_string[0] == spe_string[0].upper()
            # if names[1] in  ["sp.", "genomosp.", "genosp."]:
            if spe_match1 or spe_match2 or spe_string[0] == spe_string[0].upper():
                name_level = "genus"
                species = names[0]
            else:
                name_level = "species"
                species = " ".join(names[:2])
            ### strain and isolate may be empty
            if infraspecific_name:
                strain = infraspecific_name.split("=")[1]
            else:
                strain = isolate
            ### get the record
            record = (accession, asm_name, species, strain, name_level, bioproject, biosample, taxid, assembly_level,  gbrs_paired_asm, ftp_path, genome_name)
            Records.append(record)
    in_h.close()
    return Records

def record_to_db(assembly_summary, out_db):
    conn = sqlite3.connect(out_db)
    cursor = conn.cursor()
    cursor.execute("""
        CREATE TABLE AccessionSpecies
        (accession text,
        asm_name text,
        species text,
        strain text,
        taxa_level text,
        bioproject text,
        biosample text,
        taxid text,
        assembly_level text,
        gbrs_paired_asm text,
        ftp_path text,
        genome_name text);
    """)
    Records = parse_assembly_summary(assembly_summary)
    for r in Records:
        cursor.execute("""
            INSERT INTO AccessionSpecies
            (accession, asm_name,
            species, strain,
            taxa_level,
            bioproject, biosample,
            taxid, assembly_level,
            gbrs_paired_asm, ftp_path, genome_name)
            VALUES
            (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?);
        """, (r))
    ### Create index
    cursor.execute("""
        CREATE INDEX species2infor
        on AccessionSpecies (species);
    """)
    cursor.execute("""
        CREATE INDEX gname2infor
        on AccessionSpecies (genome_name);
    """)
    cursor.close()
    conn.commit()
    conn.close()

def main():
    parser = argparse.ArgumentParser(description="Get the information of accession and species of bacterial and them convert to database.")
    parser.add_argument("-i", "--input", help="Input file containing the assembly_summary information.")
    parser.add_argument("-d", "--database", help="The output database.")
    args = parser.parse_args()
    record_to_db(args.input, args.database)

if __name__ == "__main__":
    main()
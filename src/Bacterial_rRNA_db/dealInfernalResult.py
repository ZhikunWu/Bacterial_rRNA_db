#!/usr/bin/env python
import re
import argparse
from tinyfasta import FastaParser
import os
import sys
import sqlite3


#usage: python ~/github/Bacterial_rRNA_db/src/Bacterial_rRNA_db/dealInfernalResult.py --indir /home/wzk/database/Bacterial_genome/test_201805 --evalue 1e-20 --length 1400 --database Bacterial_rRNA.db

__author__ = "Zhikun Wu"
__email__ = "598466208@qq.com"
__date__ = "2018.05.25"

def parse_infernal_result(infernal_file, evalue_thr, length_thr):
    """
    infernal_file example:
    # query CM file:                         /home/wzk/database/Bacterial_genome/infernal_test/RF00177.seed.cm
    # target sequence database:              /home/wzk/database/Bacterial_genome/infernal_parallel/GCF_900217795.1_IMG-taxon_2740891832_annotated_assembly_genomic.fna
    # output directed to file:               /home/wzk/database/Bacterial_genome/infernal_parallel/GCF_900217795.1_IMG-taxon_2740891832_annotated_assembly_genomic.txt
    # number of worker threads:              40
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    Query:       SSU_rRNA_bacteria  [CLEN=1533]
    Accession:   RF00177
    Description: Bacterial small subunit ribosomal RNA
    Hit scores:
     rank     E-value  score  bias  sequence           start    end   mdl trunc   gc  description
     ----   --------- ------ -----  ----------------- ------ ------   --- ----- ----  -----------
      (1) !         0 1666.8   6.6  NZ_OBMQ01000025.1     36   1582 +  cm    no 0.54  Lysinibacillus xyleni strain JC22, whole genome shotgun sequ
      (2) !   2.3e-14   53.4   0.0  NZ_OBMQ01000019.1     60      1 -  cm    3' 0.52  Lysinibacillus xyleni strain JC22, whole genome shotgun sequ


    Hit alignments:
    >> NZ_OBMQ01000025.1  Lysinibacillus xyleni strain JC22, whole genome shotgun sequence
     rank     E-value  score  bias mdl mdl from   mdl to       seq from      seq to       acc trunc   gc
     ----   --------- ------ ----- --- -------- --------    ----------- -----------      ---- ----- ----
      (1) !         0 1666.8   6.6  cm        1     1533 []          36        1582 + .. 0.99    no 0.54

    """
    evalue_thr = float(evalue_thr)
    length_thr = int(length_thr)
    if os.path.exists(infernal_file):
        in_h = open(infernal_file, "r")
    else:
        gz_file = infernal_file + ".gz"
        if os.path.exists(gz_file):
            cmd = "gzip -dc %s > %s" % (gz_file, infernal_file)
            os.system(cmd)
            in_h = open(infernal_file, "r")

    record = []
    for line in in_h:
        line = line.strip()
        lines = line.split()
        if line.startswith("# target sequence database"):
            genome = lines[-1]
        elif line.startswith("Hit scores:"):
            while True:
                nextline = in_h.__next__()
                nextline = nextline.strip()
                if nextline.startswith("Hit alignments:"):
                    break
                rank = re.findall("^\(\d+\)", nextline)
                nextlines = nextline.split()
                if nextline and rank:
                    rank = nextlines[0]
                    evalue, score, bias, sequence, start, end, strand, mdl, trunc, gc = nextlines[2:12]
                    description = ' '.join(nextlines[12:])
                    evalue = float(evalue)
                    score = float(score)
                    start = int(start)
                    end = int(end)
                    gc = float(gc)
                    seqLen = abs(start - end) + 1
                    if seqLen >= length_thr and evalue <= evalue_thr:
                        record.append((genome, sequence, strand, start, end, gc, description))
    return record

def reverse_complement_seq(seq):
    """
    argv:
    seq_srting:
        GCTGAAATCGTATGGAAACGAAT
    return:
    seq_string:
        ATTCGTTTCCATACGATTTCAGC
    """
    seq = seq.lower()
    seq = seq.replace("a", "T")
    seq = seq.replace("t", "A")
    seq = seq.replace("g", "C")
    seq = seq.replace("c", "G")
    new_seq = seq[::-1]
    return new_seq

def get_seq(genome, target_contig, strand, start, end):
    """
    get target seq based on the given contig, strand, start and end information

    ### if strand is "-", the start and end position may be opposite.
    (1) !         0 1666.8   6.6  cm        1     1533 []          36        1582 + .. 0.99    no 0.54
    (2) !         0 1600.7  37.7  NZ_LT907982.1 3539179 3537661 -  cm    no 0.59  

    ### some file may contain two or more copies:
    rank     E-value  score  bias  sequence        start     end   mdl trunc   gc  description
 ----   --------- ------ -----  ------------- ------- -------   --- ----- ----  -----------
    (1) !         0 1600.7  37.7  NZ_LT907982.1 3207275 3208793 +  cm    no 0.59  Jatrophihabitans sp. GAS493 genome assembly, chromosome: I
    (2) !         0 1600.7  37.7  NZ_LT907982.1 3539179 3537661 -  cm    no 0.59  Jatrophihabitans sp. GAS493 genome  assembly, chromosome: I
    """
    start = int(start)
    end = int(end)
    target_seq = "-"
    for record in FastaParser(genome):
        desc = str(record.description)
        desc = desc.lstrip(">").split()[0]
        if desc == target_contig:
            seq = str(record.sequence)
            if strand == "+":
                target_seq = seq[start-1 : end]
            elif strand == "-":
                target_seq = seq[end-1 : start]
                target_seq = reverse_complement_seq(target_seq)
            elif strand != "-" or strand != "+":
                print("Please check the and make sure the strand %s in the target contig %s to be '+' of '-'." % (strand, target_contig))
            break
    return target_seq


def get_infernal_search_seq(infernal_file, evalue_thr, length_thr):
    record = parse_infernal_result(infernal_file, evalue_thr, length_thr)
    recordLen = len(record)
    seqResult = []
    for r in record:
        genome, seq_desc, strand, start, end, gc, description = r
        target_seq = get_seq(genome, seq_desc, strand, start, end)
        genome_base = os.path.basename(genome)
        seqResult.append((genome_base, seq_desc, strand, start, end, gc, target_seq))
    GeneCopy = (genome_base, recordLen)
    return seqResult, GeneCopy

def result_to_db(indir, evalue_thr, length_thr, out_db):
    """
    output database example:
    sqlite> select * from BacterialSeq limit 1;
    genome_name|contig_id|contig_strand|contig_start|contig_end|gc_ratio|sequence
    GCF_900241005.1_PRJEB22714_genomic.fna|NZ_OEST01000016.1|-|5062|3536|0.52|ACAATGAAGAGTTTGATCCTGGCTCAGGATGAACGCTAGCTACAGGCTTAACACATGCAAGTCGAGGGGCATCAGGGCATTAGCTTGCTAATGCCGCTGGCGACCGGCGCACGGGTGAGTAACACGTATCCAACCTGCCTTCAACTCAGGGATAGCCTTTCGAAAGAAAGATTAATACCTGATAGTATGTATTTTCCGCATGGCCTTTACATTAAAGGGTTCCGGTTGAAGATGGGGATGCGTTCCATTAGATAGTTGGCGGGGTAACGGCCCACCAAGTCCACGATGGATAGGGGTTCTGAGAGGAAGGTCCCCCACATTGGAACTGAGACACGGTCCAAACTCCTACGGGAGGCAGCAGTGAGGAATATTGGTCAATGGGCGCTAGCCTGAACCAGCCAAGTAGCGTGAAGGATGACTGCCCTATGGGTTGTAAACTTCTTTTATACGGGAATAAAGTTTGGCACGTGTGCCTTTTTGCATGTACCGTATGAATAAGGATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATCCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGAGCGTAGGCGGATGCTTAAGTCAGTTGTGAAAGTTTGCGGCTCAACCGTAAAATTGCAGTTGATACTGGGTGTCTTGAGTGCAGTATAGGCAGGCGGAATTCGTGGTGTAGCGGTGAAATGCTTAGATATCACGAAGAACTCCGATTGCGAAGGCAGCTTGCTGGACTGTAACTGACGCTGATGCTCGAAAGTGTGGGTATCAAACAGGATTAGATACCCTGGTAGTCCACACAGTAAACGATGAATACTCGCTGTTTGCGATATACGGTAAGCGGCCAAGCGAAAGCGTTAAGTATTCCACCTGGGGAGTACGCCGGCAACGGTGAAACTCAAAGGAATTGACGGGGGCCCGCACAAGCGGAGGAACATGTGGTTTAATTCGATGATACGCGAGGAACCTTACCCGGGCTTGAATTGCAAATGAATATAGTGGAGACATTATAGCCGCAAGGCATTTGTGAAGGTGCTGCATGGTTGTCGTCAGCTCGTGCCGTGAGGTGTCGGCTTAAGTGCCATAACGAGCGCAACCCTTATCTATAGTTACTATCAGGTAATGCTGAGGACTCTATGGAGACTGCCGTCGTAAGATGTGAGGAAGGTGGGGATGACGTCAAATCAGCACGGCCCTTACGTCCGGGGCTACACACGTGTTACAATGGGGGGTACAGAAGGCAGCTACACGGCGACGTGATGCTAATCCCGAAAACCTCTCTCAGTTCGGATTGGAGTCTGCAACCCGACTCCATGAAGCTGGATTCGCTAGTAATCGCGCATCAGCCACGGCGCGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCAAGCCATGAAAGCCGGGGGTACCTGAAGTGCGTAACCGCAAGGAGCGCCCTAGGGTAAAACTGGTGATTGGGGCTAAGTCGTAACAAGGTAGCCGTACCGGAAGGTGCGGCTGGAACACCTCCTTT

    sqlite> select * from GeneCopy  limit 1;
    genome_name|copy
    GCF_900241005.1_PRJEB22714_genomic.fna|1

    """
    ### Create the database
    conn = sqlite3.connect(out_db)
    cursor = conn.cursor()
    cursor.execute("""
        CREATE TABLE BacterialSeq
        (genome_name text,
        contig_id text,
        contig_strand text,
        contig_start integer,
        contig_end integer,
        gc_ratio real,
        sequence text)
    """)

    cursor.execute("""
        CREATE TABLE GeneCopy
        (genome_name text,
        copy text)
    """)

    if not os.path.isdir(indir):
        print("Please check whether the input %s is a directory." % indir)
        sys.exit(1)
    else:
        if not indir.endswith("/"):
            indir += "/"
        files = os.listdir(indir)
        for f in files:
            f = indir + f
            seqResult, GeneCopy = get_infernal_search_seq(f, evalue_thr, length_thr)
            ### Insert genome name and rRNA copies
            cursor.execute("""
                INSERT INTO GeneCopy
                (genome_name, copy)
                VALUES
                (?, ?)
            """, (GeneCopy))
            ### Insert all copies of target rRNA in the given genome
            for s in seqResult:
                cursor.execute("""
                    INSERT INTO BacterialSeq
                    (genome_name, 
                    contig_id, contig_strand,
                    contig_start, contig_end, 
                    gc_ratio, sequence)
                    VALUES
                    (?, ?, ?, ?, ?, ?, ?);
                """, (s))
    cursor.close()
    conn.commit()
    conn.close()



def main():
    parser = argparse.ArgumentParser(description="Get the records for the infernal result.")
    parser.add_argument("-i", "--indir", help="The input file or directory.")
    parser.add_argument("-l", "--length", help="The length threshold.")
    parser.add_argument("-e", "--evalue", help="The evalue threshold.")
    parser.add_argument("-d", "--database", help="The out database.")
    args = parser.parse_args()
    result_to_db(args.indir, args.evalue, args.length, args.database)

if __name__ == "__main__":
    main()
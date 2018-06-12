#!/usr/bin/env python
import collections
import argparse
import re

#usage: python ~/github/Bacterial_rRNA_db/src/Bacterial_rRNA_db/FiltSAM.py --sam fastqjoin.join.fastq.sam --lenThreshold 0.95 --matchThreshold 0.9 --out  fastqjoin.join.fastq.sam_out

__author__ = "Zhikun Wu"
__email__ = "598466208@qq.com"
__date__ = "2018.06.12"

def given_tag(tags, query):
    """
    a = ["AS:i:-35", "XS:i:-39", "XN:i:0", "XM:i:6", "XO:i:0", "XG:i:0", "NM:i:6", "MD:Z:75G0C1G37C14C152G8", "YT:Z:UU"]
    b = given_tag(a, "MD:Z")
    print(b)
    75G0C1G37C14C152G8
    """
    value = ""
    for t in tags:
        ts = t.split(":")
        t_name = ":".join(ts[:-1])
        t_value = ts[-1]
        if t_name == query:
            value = t_value
            break
    return value

def match_cigar_length(match_tag):
    """
    Get the match number for the cigar.
    match_number("110M4D3M4I176M")
    289
    """
    count = 0
    match = re.findall("(\d+\w)", match_tag)
    for m in match:
        if m.endswith("M"):
            number = int(m.rstrip("M"))
            count += number
    return count

def MD_match_number(MD_value):
    """
    Get the match number based on the MD value
    a = "75G0C1G37C14C152G8"
    MD_match_number(a)
    287
    """
    count = 0
    match = re.findall("(\d+)", MD_value)
    match = map(int, match)
    count = count + sum(match)
    return count

def filt_for_sam(sam_file, lenThreshold, matchThreshold, out_file):
    """
    sam_file: 293866 records with 0 flag
    $ grep -v '^@' fastqjoin.join.fastq.sam | cut -f 2 | sort | uniq -c 
     293866 0
         13 16
     109922 4
    
    out_file: 289410 records pass the threshold 0.9
    $ wc -l  fastqjoin.join.fastq.sam_out
    289410 fastqjoin.join.fastq.sam_out
    """
    matchThreshold = float(matchThreshold)
    lenThreshold = float(lenThreshold)
    in_h = open(sam_file, 'r')
    out_h = open(out_file, 'w')
    for line in in_h:
        line = line.strip()
        if line.startswith("@"):
            continue
        else:
            lines = line.split("\t")
            query, flag, subject, position, mapq, cigar, qmate, matePos, insertLen, seq, quality = lines[:11]
            flag = int(flag)
            seqLen = len(seq)
            Tags = lines[11:]
            MDValue = given_tag(Tags, "MD:Z")
            if flag == 0:
                ### Get the match length, it contain seq with match and mismatch
                MatchLen = match_cigar_length(cigar)
                MatchNumber = MD_match_number(MDValue)
                ### Filt with the match threshold
                if MatchLen >= seqLen * lenThreshold and MatchNumber >=  MatchLen *matchThreshold:
                    out_h.write("%s\t%s\n" % (query, subject))
    in_h.close()
    out_h.close()


def main():
    parser = argparse.ArgumentParser(description="Get the effective mapping record based on the SAM format file with given match length threshold.")
    parser.add_argument("-s", "--sam", help="The input sam file.")
    parser.add_argument("-l", "--lenThreshold", help="The length threshold for sequence, such as 0.95.")
    parser.add_argument("-t", "--matchThreshold", help="The threshold for match number of the sequence length, such as 0.9.")
    parser.add_argument("-o", "--out", help="The output file.")
    args = parser.parse_args()
    filt_for_sam(args.sam, args.lenThreshold, args.matchThreshold, args.out)

if __name__ == "__main__":
    main()


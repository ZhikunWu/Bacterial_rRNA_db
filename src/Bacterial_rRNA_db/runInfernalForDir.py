#!/usr/bin/env python
import os
import sys
import argparse

#usage: python  ~/github/Bacterial_rRNA_db/src/Bacterial_rRNA_db/runInfernalForDir.py --indir test_201805 --temp_dir test_temp --seed /home/wzk/database/Bacterial_genome/infernal_test/RF00177.seed.cm  --outdir test_out

__author__ = "ZhikunWu"
__email__ = "598466208@qq.com"
__date__ = "2018.05.25"

def check_out_dir(out_dir):
    if os.path.exists(out_dir):
        if not os.path.isdir(out_dir):
            os.makedirs(out_dir)
    else:
        os.makedirs(out_dir)
    if not out_dir.endswith("/"):
        out_dir += "/"
    return out_dir

def get_genome_file(in_dir, outdir):
    ### Check the in directory and out directory
    if not os.path.isdir(in_dir):
        print("Please check whether input paramater %s is a directory." % in_dir)
    else:
        if not in_dir.endswith("/"):
            in_dir += "/"

    out_dir = check_out_dir(outdir)

    ### copy the files in the in directory.
    files = os.listdir(in_dir)
    for f in files:
        f = in_dir + f
        if f.endswith(".gz") or f.endswith(".gzip"):
            f_base = os.path.basename(f)
            new_f = '.'.join(f_base.split(".")[:-1])
            new_ff = out_dir + new_f
            cmd = "gzip -dc %s > %s" % (f, new_ff)
        else:
            f_base = os.path.basename(f)
            new_f = out_dir + f_base
            cmd = "cp %s  %s" % (f, new_f)
        os.system(cmd)


def run_infernal_for_dir(out_dir, seed_file, resultDir):
    result_dir = check_out_dir(resultDir)

    files = os.listdir(out_dir)
    if not out_dir.endswith("/"):
        out_dir += "/"
    for f in files:
        f = out_dir + f
        f_base = os.path.basename(f)
        out_f = result_dir + f_base
        cmd = "cmsearch %s %s > %s" % (seed_file, f, out_f)
        try:
            os.system(cmd)
        except FileNotFoundError:
            print("Please check whether the tool cmsearch (of infernal) exists in your PATH.")

def main():
    parser = argparse.ArgumentParser(description="Get the rRNA information based on software infernal nad the genome sequence.")
    parser.add_argument("-i", "--indir", help="The input directory containing the bacterial genome sequences.")
    parser.add_argument("-t", "--temp_dir", help="The input temp directory containing uncompressed genomes. ")
    parser.add_argument("-s", "--seed", help="The seed file for infernal.")
    parser.add_argument("-o", "--outdir", help="The out directory.")
    args = parser.parse_args()
    get_genome_file(args.indir, args.temp_dir)
    run_infernal_for_dir(args.temp_dir, args.seed, args.outdir)

if __name__ == "__main__":
    main()


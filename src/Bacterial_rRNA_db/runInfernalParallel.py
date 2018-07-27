#!/usr/bin/env python
import os
import sys
import argparse
from multiprocessing import Pool


#usage: python ~/github/zkwu/kcMeta/pipeline/runInfernal.py --indir /home/wzk/database/Bacterial_genome/test --seed seeds/RF00177.seed.cm --threads 3 --outdir /home/wzk/database/Bacterial_genome/test1

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



def run_infernal(seed, in_file, out_file):
    cmd = "cmsearch %s %s > %s" % (seed, in_file, out_file)
    os.system(cmd)

def gunzip_file(file):
    cmd = "gunzip %s" % file
    if file.strip().endswith(".gz"):
        os.system(cmd)

def run_infernal_for_dir(in_dir, seed_file, resultDir, threads):
    result_dir = check_out_dir(resultDir)
    in_dir = check_out_dir(in_dir)
    files = os.listdir(in_dir)
    inFiles = [in_dir + f for f in files]
    ### run with mutiple threads
    threads = int(threads)
    print("@Start to uncompress the genome files which end with '.gz'.")    
    pool = Pool(threads)
    for i in inFiles:
        result = pool.apply_async(gunzip_file, (i,))
    pool.close()
    pool.join()
    if result.successful():
        print("Successfully uncompress the target files!\n")

    newFiles = os.listdir(in_dir)
    outFiles = [result_dir + f for f in newFiles]
    newFiles = [in_dir + f for f in newFiles]   
    print("@Start to find the taget gene using cmsearch based on the seed of %s." % seed_file)
    pool = Pool(threads)
    for i, j in zip(newFiles, outFiles):
        result = pool.apply_async(run_infernal, (seed_file, i, j))
    pool.close()
    pool.join()
    if result.successful():
        print("Successfully find the target genes!\n")



def main():
    parser = argparse.ArgumentParser(description="Get the rRNA information based on software infernal nad the genome sequence.")
    parser.add_argument("-i", "--indir", help="The input directory containing the bacterial genome sequences.")
    parser.add_argument("-t", "--threads", help="The threads. ")
    parser.add_argument("-s", "--seed", help="The seed file for infernal.")
    parser.add_argument("-o", "--outdir", help="The out directory.")
    args = parser.parse_args()
    run_infernal_for_dir(args.indir, args.seed,  args.outdir, args.threads)

if __name__ == "__main__":
    main()

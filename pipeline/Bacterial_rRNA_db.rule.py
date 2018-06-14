#!/usr/bin/env python
import yaml
import os
#usage: snakemake -s ~/github/Bacterial_rRNA_db/pipeline/Bacterial_rRNA_db.rule.py -j 30 --configfile ~/github/Bacterial_rRNA_db/pipeline/Bacterial_rRNA_db.rule.yaml


__author__ = 'Zhikun Wu'
__email__ = '598466208@qq.com'
__date__ = '2018.06.11'


IN_PATH = '/home/wzk/database/Bacterial_genome'
SRC_DIR = '/home/wzk/github/Bacterial_rRNA_db/src/Bacterial_rRNA_db'

SAMPLES = config["SAMPLES"]
THREADS = config["threads"]
SILVA_DB_DIR = config["SILVA_DB_DIR"]
# NCBI_DB_DIR = config["NCBI_DB_DIR"]



rule all:
    input:
        expand(IN_PATH + "/{sample}/species_infernal/temp", sample=SAMPLES),
        expand(IN_PATH + "/{sample}/database/Bacterial_16SrRNA.db", sample=SAMPLES),
        IN_PATH + "/SILVA_132_SSUParc_tax_silva_DNA_species.db",
        expand(IN_PATH + "/{sample}/database/NCBI_accession_16SrRNA.fasta", sample=SAMPLES),
        IN_PATH + "/fastqjoin.join_sample.sam",
        expand(IN_PATH + "/{sample}/temp1", sample=SAMPLES),
        expand(IN_PATH + "/{sample}/temp2", sample=SAMPLES),
        expand(IN_PATH + "/{sample}/database/NCBI_accession_16SrRNA.1.bt2", sample=SAMPLES),
        expand(IN_PATH + "/{sample}/mapping/fastqjoin.join.sam", sample=SAMPLES),
        expand(IN_PATH + "/{sample}/mapping/fastqjoin.join_record.txt", sample=SAMPLES),
        expand(IN_PATH + "/{sample}/mapping/fastqjoin.join_texonomy.txt", sample=SAMPLES),
        "/home/wzk/database/SILVA/SILVA_132_SSUParc_tax_silva_genus.fasta",
        SILVA_DB_DIR +  "/SILVA_132_SSUParc_tax_silva_genus.1.bt2",
        SILVA_DB_DIR + "/SILVA_132_SSUParc_tax_silva_species.1.bt2",


################################ For NCBI database  #####################################
#########################################################################################

########################## Extract the 16S rRNA gene sequence from genome ###################
rule runInfernal:
    input:
        indir = IN_PATH + "/species",
    output:
        out = IN_PATH + "/{sample}/species_infernal/temp",
    params:
        runInfernalForDir = SRC_DIR + "/runInfernalForDir.py",
        seed = config["seed"],
        temp_dir = IN_PATH + "/{sample}/species_temp",
        outdir = IN_PATH + "/{sample}/species_infernal",
    log:
        IN_PATH + "/log/runInfernal_{sample}.log"
    run:
        makedir(params.outdir)
        makedir(params.temp_dir)
        shell("python  {params.runInfernalForDir} --indir {input.indir} --temp_dir {params.temp_dir} --seed {params.seed}  --outdir {params.outdir} >{log} 2>&1")
        shell("touch {output.out}")


rule dealInfernalResult:
    input:
        indir = rules.runInfernal.params.outdir,
        temp = rules.runInfernal.output.out,
    output:
        database = IN_PATH + "/{sample}/database/Bacterial_16SrRNA.db",
    params:
        dealInfernalResult = SRC_DIR + "/dealInfernalResult.py",
        evalue = config["evalue"],
        length = config["length"], #The length of 16S rRNA is normally 1500
    log:
        IN_PATH + "/log/dealInfernalResult_{sample}.log"
    run:
        shell("python {params.dealInfernalResult} --indir {input.indir} --evalue {params.evalue} --length {params.length} --database {output.database} >{log} 2>&1")

rule assembly2db:
    input:
        assembly = IN_PATH + "/assembly_summary_uniq.txt",
    output:
        assembly_db = IN_PATH + "/assembly_summary_uniq.db",
    params:
        assembly2db = SRC_DIR + "/assembly2db.py",
    log:
        IN_PATH + "/log/assembly2db.log"
    run:
        shell("python {params.assembly2db} --input {input.assembly} --database {output.assembly_db} >{log} 2>&1")


rule SILVASpecies:
    input:
        taxo_record = SILVA_DB_DIR + "/SILVA_132_SSUParc_tax_silva.fasta",
    output:
        taxo_seq = SILVA_DB_DIR + "/SILVA_132_SSUParc_tax_silva_species.fasta",
        taxo_db = SILVA_DB_DIR + "/SILVA_132_SSUParc_tax_silva_species.db",
    params:
        getSilvaSpecies = SRC_DIR + "/getSilvaSpecies.py",
        lenThreshold = config["length"]
    log:
        IN_PATH + "/log/SILVASpecies.log"        
    run:
        shell("python {params.getSilvaSpecies} --input {input.taxo_record}  --out {output.taxo_seq} --database {output.taxo_db} --lenThreshold {params.lenThreshold} >{log} 2>&1")


rule NCBITaxoSeq:
    input:
        accession = rules.assembly2db.output.assembly_db,
        taxonomy = rules.SILVASpecies.output.taxo_db,
        sequence = rules.dealInfernalResult.output.database,
    output:
        fa = IN_PATH + "/{sample}/database/NCBI_accession_16SrRNA.fasta",
        copy = IN_PATH + "/{sample}/database/NCBI_accession_16SrRNA_geneCopy.xls",
    params:
        NCBI16STaxoSeq = SRC_DIR + "/NCBI16STaxoSeq.py",
    log:
        IN_PATH + "/log/NCBITaxoSeq_{sample}.log"
    run:
        shell("python {params.NCBI16STaxoSeq} --accession {input.accession} --taxonomy {input.taxonomy} --sequence {input.sequence} --out {output.fa} --geneCopy {output.copy} >{log} 2>&1")
########################################################################################



############################### Mapping reads using bowtie2 ##############################
rule bowtie2build:
    input:
        fa = rules.NCBITaxoSeq.output.fa,
    output:
        bt2 = IN_PATH + "/{sample}/database/NCBI_accession_16SrRNA.1.bt2",
    threads:
        THREADS
    params:
        ref_base = IN_PATH + "/{sample}/database/NCBI_accession_16SrRNA",
    log:
        IN_PATH + "/log/bowtie2build_{sample}.log"
    run:
        shell("bowtie2-build --threads {threads} {input.fa}  {params.ref_base} >{log} 2>&1")

rule bowtie2Align:
    input:
        bt2 = rules.bowtie2build.output.bt2,
        ref_base = rules.bowtie2build.params.ref_base,
        fastq = IN_PATH + "/fastqjoin.join.fastq",
    output:
        sam = IN_PATH + "/{sample}/mapping/fastqjoin.join.sam",
    threads:
        THREADS
    log:
        IN_PATH + "/log/bowtie2Align_{sample}.log"
    run:
        shell("bowtie2 -p {threads} -x {input.ref_base} -U {input.fastq} -S {output.sam} >{log} 2>&1")

rule FiltSAM:
    input:
        sam = rules.bowtie2Align.output.sam,
    output:
        record = IN_PATH + "/{sample}/mapping/fastqjoin.join_record.txt",
    params:
        FiltSAM = SRC_DIR + "/FiltSAM.py",
        lenThreshold = 0.95,
        matchThreshold = 0.9,
    log:
        IN_PATH + "/log/FiltSAM_{sample}.log"
    run:
        shell("python {params.FiltSAM} --sam {input.sam} --lenThreshold {params.lenThreshold} --matchThreshold {params.matchThreshold} --out {output.record} >{log} 2>&1")

rule TaxoCount:
    input:
        fasta = rules.NCBITaxoSeq.output.fa,
        record = rules.FiltSAM.output.record,
    output:
        taxo = IN_PATH + "/{sample}/mapping/fastqjoin.join_texonomy.txt",
    params:
        SpeciesNum = SRC_DIR + "/SpeciesNum.py",
    log:
        IN_PATH + "/log/TaxoCount_{sample}.log"       
    run:
        shell("python {params.SpeciesNum} --fasta {input.fasta} --mapping {input.record} --out {output.taxo} >{log} 2>&1") 

#############################################################################################


################## Compare the sequence difference among species within genus ############
rule GenusSeq:
    input:
        fa = rules.NCBITaxoSeq.output.fa,
    output:
        fa = IN_PATH + "/{sample}/temp1",
    params:
        SameGenusSeq = SRC_DIR + "/SameGenusSeq.py",
        outdir = IN_PATH + "/{sample}/database/NCBI_accession_16SrRNA",
    log:
        IN_PATH + "/log/GenusSeq_{sample}.log"        
    run:
        shell("python {params.SameGenusSeq} --input {input.fa} --outdir {params.outdir} >{log} 2>&1")
        shell("touch {output.fa}")

rule GenusDiffStats:
    input:
        fa = rules.GenusSeq.output.fa,
    output:
        temp = IN_PATH + "/{sample}/temp2",
    params:
        GenusSeqDiffStats = SRC_DIR + "/GenusSeqDiffStats.py",
        indir = rules.GenusSeq.params.outdir,
        outdir = IN_PATH + "/{sample}/database/NCBI_accession_16SrRNA_compare",
    log:
        IN_PATH + "/log/GenusDiffStats_{sample}.log"  
    run:
        shell("python {params.GenusSeqDiffStats} --indir {params.indir} --outdir {params.outdir} >{log} 2>&1")
        shell("touch {output.temp}")
##################################################################################################


################################ Mapping reads using bwa #############################
# rule BWAindex:
#     input:
#         taxo_seq = rules.SILVASpecies.output.taxo_seq,
#     output:
#         bwt = IN_PATH + "/SILVA_132_SSUParc_tax_silva_DNA_species.fasta.bwt",
#     log:
#         IN_PATH + "/log/BWAindex.log"
#     run:
#         shell("bwa index {input.taxo_seq} >{log} 2>&1")

# rule BWAalign:
#     input:
#         fq = IN_PATH + "/fastqjoin.join_sample.fastq",
#         ref = rules.SILVASpecies.output.taxo_seq,
#         bwt = rules.BWAindex.output.bwt,
#     output:
#         sam = IN_PATH + "/fastqjoin.join_sample.sam",
#     threads:
#         THREADS
#     log:
#         IN_PATH + "/log/BWAalign.log"
#     run:
#         shell("bwa mem -t {threads} {input.ref} {input.fq} > {output.sam} 2>{log}")

#########################################################################################








################################ For SILVA database ##################################
######################################################################################

################################# Get 16S reference with genus level #################
rule SILVAGenus:
    input:
        fasta = SILVA_DB_DIR + "/SILVA_132_SSUParc_tax_silva.fasta",
    output:
        genus = SILVA_DB_DIR + "/SILVA_132_SSUParc_tax_silva_genus.fasta",
        genus_db =  SILVA_DB_DIR + "/SILVA_132_SSUParc_tax_silva_genus.db",
    params:
        getSilvaGenus = SRC_DIR + "/getSilvaGenus.py",
        lenThreshold = config["length"],
    log: 
        IN_PATH + "/log/SILVAGenus.log"
    run:
        shell("python {params.getSilvaGenus} --input {input.fasta} --out {output.genus} --database {output.genus_db} --lenThreshold {params.lenThreshold} >{log} 2>&1")

rule GenusBowtie2Build:
    input:
        fa = rules.SILVAGenus.output.genus,
    output:
        bt2 = SILVA_DB_DIR +  "/SILVA_132_SSUParc_tax_silva_genus.1.bt2",
    threads:
        THREADS
    params:
        ref_base = SILVA_DB_DIR +  "/SILVA_132_SSUParc_tax_silva_genus",
    log:
        IN_PATH + "/log/GenusBowtie2Build.log"
    run:
        shell("bowtie2-build --threads {threads} {input.fa}  {params.ref_base} >{log} 2>&1")



rule SpeciesBowtie2Build:
    input:
        fa = SILVA_DB_DIR + "/SILVA_132_SSUParc_tax_silva_species.fasta",
    output:
        bt2 = SILVA_DB_DIR + "/SILVA_132_SSUParc_tax_silva_species.1.bt2",
    threads:
        THREADS 
    params:
        ref_base = SILVA_DB_DIR + "/SILVA_132_SSUParc_tax_silva_species",
    log:
        IN_PATH + "/log/SpeciesBowtie2Build.log"        
    run:
        shell("bowtie2-build --threads {threads} {input.fa}  {params.ref_base} >{log} 2>&1")        


########################################################################################

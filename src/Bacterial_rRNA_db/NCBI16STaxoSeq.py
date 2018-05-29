#!/usr/bin/env python
import sqlite3
import argparse

#usage: python /home/wzk/github/Bacterial_rRNA_db/src/Bacterial_rRNA_db/NCBI16STaxoSeq.py --accession /home/wzk/database/Bacterial_genomeassembly_summary.db --taxonomy /home/wzk/database/Bacterial_genome/SILVA_132_SSUParc_tax_silva_DNA_species.db --sequence /home/wzk/database/Bacterial_genome/Bacterial_16S/database/Bacterial_16SrRNA.db --out /home/wzk/database/Bacterial_genome/Bacterial_16S/database/NCBI_accession_16SrRNA.fasta  --geneCopy /home/wzk/database/Bacterial_genome/Bacterial_16S/database/NCBI_accession_16SrRNA_geneCopy.xls

__author__ = "Zhikun Wu"
__email__ = "598466208@qq.com"
__date__ = "2018.05.29"



def combine_assembly_silva_16S(AccessionSpecies_db, assembly_db, BacterialSeq_db, out_file, gene_copy):
    """
    AccessionSpecies:
    accession|asm_name|species|strain|taxa_level|bioproject|biosample|taxid|assembly_level|gbrs_paired_asm|ftp_path|genome_name
    GCF_000010525.1|ASM1052v1|Azorhizobium caulinodans|ORS 571|species|PRJNA224116|SAMD00060925|438753|Complete Genome|GCA_000010525.1|ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/010/525/GCF_000010525.1_ASM1052v1|GCF_000010525.1_ASM1052v1_genomic.fna

    SpeciesTaxonomy:
    Kindom|Phylum|Class|Orders|Family|Genus|Species
    Bacteria|Firmicutes|Bacilli|Bacillales|Bacillaceae|Thalassobacillus|Thalassobacillus acidophilus

    BacterialSeq:
    genome_name|contig_id|contig_strand|contig_start|contig_end|gc_ratio|seq
    GCF_000988235.1_ASM98823v1_genomic.fna|NZ_LAYS01000008.1|-|4078|2561|0.58|ATGT

    GeneCopy
    genome_name|copy
    GCF_000988235.1_ASM98823v1_genomic.fna|1
    """
    conn = sqlite3.connect(AccessionSpecies_db)
    cursor = conn.cursor()
    cursor.execute("""
        ATTACH DATABASE ? as SpeciesTaxonomy;
    """, (assembly_db,))

    cursor.execute("""
        ATTACH DATABASE ? as BacterialSeq;
    """, (BacterialSeq_db,))
    ### get the accession, taxonomy and sequence
    cursor.execute("""
        SELECT AccessionSpecies.accession,
        BacterialSeq.contig_id,
        BacterialSeq.contig_start, BacterialSeq.contig_end,
        SpeciesTaxonomy.Kindom, SpeciesTaxonomy.Phylum,
        SpeciesTaxonomy.Class, SpeciesTaxonomy.Orders,
        SpeciesTaxonomy.Family, SpeciesTaxonomy.Genus,
        SpeciesTaxonomy.Species, AccessionSpecies.strain, 
        BacterialSeq.genome_name,  
        BacterialSeq.sequence
        from AccessionSpecies
        left join SpeciesTaxonomy.SpeciesTaxonomy SpeciesTaxonomy
        on SpeciesTaxonomy.Species =  AccessionSpecies.species
        left join BacterialSeq.BacterialSeq BacterialSeq
        on BacterialSeq.genome_name = AccessionSpecies.genome_name;
    """)

    Records = cursor.fetchall()

    ### output the records
    out_h = open(out_file,"w")
    for r in Records:
        accession_desc = ','.join(map(str, r[:4]))
        taxonomy = ';'.join(map(str, r[4:11]))
        seq = r[-1]
        ######################################################
        ### sequence maybe None bacause the alignment result with infernal may filted out based on the length threshold
        if seq != None:
            out_h.write(">%s\t%s\n%s\n" % (accession_desc, taxonomy, seq))
    out_h.close()

    ### get the taxonomy and gene coppy
    cursor.execute("""
        SELECT 
        SpeciesTaxonomy.Kindom, SpeciesTaxonomy.Phylum,
        SpeciesTaxonomy.Class, SpeciesTaxonomy.Orders,
        SpeciesTaxonomy.Family, SpeciesTaxonomy.Genus,
        SpeciesTaxonomy.Species, AccessionSpecies.strain, 
        GeneCopy.copy
        from AccessionSpecies
        left join SpeciesTaxonomy.SpeciesTaxonomy SpeciesTaxonomy
        on SpeciesTaxonomy.Species =  AccessionSpecies.species
        left join BacterialSeq.GeneCopy GeneCopy
        on GeneCopy.genome_name = AccessionSpecies.genome_name;
    """)
    copy_h = open(gene_copy, "w")
    copy_h.write("Taxonomy\tGene_copy\n")
    Records = cursor.fetchall()
    for r in Records:
        #############################################################
        # the taxonomy in record maybe empty
        # (None, None, None, None, None, None, None, 'ATCC 15551', '2')
        copy = r[-1]
        if copy != None and r[0] != None:
            taxonomy = ';'.join(r[:7])
            copy_h.write("%s\t%s\n" % (taxonomy, copy))
    copy_h.close()


def main():
    parser = argparse.ArgumentParser(description="Combine the accession, taxonomy and sequences, then output the records.")
    parser.add_argument("-a", "--accession", help="The input accessionSpecies database.")
    parser.add_argument("-t", "--taxonomy", help="The input taxonomy database with species level.")
    parser.add_argument("-s", "--sequence", help="The input bacterial 16S rRNA sequence database.")
    parser.add_argument("-c", "--geneCopy", help="The output file with 16S rRNA gene copy.")
    parser.add_argument("-o", "--out", help="The output file.")
    args = parser.parse_args()
    combine_assembly_silva_16S(args.accession, args.taxonomy, args.sequence, args.out, args.geneCopy)

if __name__ == "__main__":
    main()






# Bacterial_rRNA_db

## Download Bacterial genomes from NCBI

Assembly summary
```
$ curl 'ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt' > assembly_summary_20180625.txt
```

Download genomes
```
$ curl 'ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt' |  awk '{FS="\t"} !/^#/ {print $20} ' | sed -r 's|(ftp://ftp.ncbi.nlm.nih.gov/genomes/all/.+/)(GCF_.+)|\1\2/\2_genomic.fna.gz|'  > genomic_file_20180625

$ wget --input genomic_file_20180625
```

## [SILVA database](https://www.arb-silva.de)

download file [ SILVA_132_SSUParc_tax_silva.fasta.gz](https://www.arb-silva.de/fileadmin/silva_databases/release_132/Exports/SILVA_132_SSUParc_tax_silva.fasta.gz)



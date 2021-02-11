# this script prepares the reference genome and annotation 
# for GRCh37 (hg19) from Gencode

# download reference genome
wget --directory-prefix=data/rawdata/ \
   ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/GRCh37_mapping/GRCh37.primary_assembly.genome.fa.gz
samtools faidx data/rawdata/GRCh37.primary_assembly.genome.fa.gz
# download the gencode annotation gtf 
wget --directory-prefix=data/rawdata/ \
   ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz
# filter for protein-coding genes 
zcat data/rawdata/gencode.v19.annotation.gtf.gz | 
   awk '$3=="gene"{print $0}' | 
   grep 'protein_coding' \
   > data/streamLined/gencode.v19.annotation.codingonly.gtf
# reformats it 
cat data/streamLined/gencode.v19.annotation.codingonly.gtf |
   awk 'BEGIN{OFS="\t"}; {print $1,$4-1,$5,$10" "$1" "$4-1" "$5,$18,$7}' | 
   tr -d '";' \
   > data/streamLined/gencode.v19.annotation.codingonly.bed
# add the sequences of each gene (from reference genome)
bedtools getfasta \
   -fi data/rawdata/GRCh37.primary_assembly.genome.fa.gz \
   -bed data/streamLined/gencode.v19.annotation.codingonly.bed \
   -name -tab \
   > data/streamLined/gencode.v19.annotation.codingonly.withsequences.bed

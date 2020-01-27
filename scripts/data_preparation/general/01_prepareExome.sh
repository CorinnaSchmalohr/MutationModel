zcat data/rawdata/gencode.v19.annotation.gtf.gz | awk '$3=="gene"{print $0}' | grep 'protein_coding' > data/procData/gencode.v19.annotation.codingonly.gtf
cat data/procData/gencode.v19.annotation.codingonly.gtf | awk 'BEGIN{OFS="\t"}; {print $1,$4-1,$5,$10" "$1" "$4-1" "$5,$18,$7}' | tr -d '";' > data/procData/gencode.v19.annotation.codingonly.bed
bedtools getfasta -fi data/rawdata/GRCh37.p11.genome.fa -bed data/procData/gencode.v19.annotation.codingonly.bed -name -tab > data/procData/gencode.v19.annotation.codingonly.withsequences.bed

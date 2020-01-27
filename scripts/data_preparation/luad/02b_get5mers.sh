mkdir data/procData/5mers
cat data/procData/luad/Muts.bed | awk 'BEGIN{OFS="\t"}; {print $1,$2-2,$3+2,$1"_"$2"_"$3}' > data/procData/5mers/Muts_5merPositions.bed

cat data/procData/luad/Muts.bed | awk 'BEGIN{OFS="\t"}; {print $1,$2-3,$3+3,$1"_"$2"_"$3}' > data/procData/5mers/Muts_7merPositions.bed

bedtools getfasta -fi data/rawdata/GRCh37.p11.genome.fa -bed data/procData/5mers/Muts_5merPositions.bed -name -tab > data/procData/5mers/5mers.txt
bedtools getfasta -fi data/rawdata/GRCh37.p11.genome.fa -bed data/procData/5mers/Muts_7merPositions.bed -name -tab > data/procData/5mers/7mers.txt

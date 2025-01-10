# translate gtf into bed
cat data/procData/codExons.gtf |
   awk 'BEGIN{OFS="\t"}; {$4=$4-3; $5=$5+2; print $1,$4,$5, $1" "$4" "$5" "$9}' | 
   bedtools getfasta -fi data/rawdata/GRCh37.p11.genome.fa \
   -bed - \
   -name -tab > data/procData/codExons.withsequences.bed


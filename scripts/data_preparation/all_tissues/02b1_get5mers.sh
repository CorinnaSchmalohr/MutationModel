
for t in luad breast skin colon ovary kidney prostate
do
   echo $t
   rm -r data/procData/${t}/context
   mkdir data/procData/${t}/context
   cat data/procData/${t}/Muts.bed | 
      awk 'BEGIN{OFS="\t"}; {print $1,$2-3,$3+3,$1"_"$2"_"$3}' \
      > data/procData/${t}/context/context.bed
   bedtools getfasta -fi data/rawdata/GRCh37.p11.genome.fa \
      -bed data/procData/${t}/context/context.bed \
      -name -tab > data/procData/${t}/context/context.txt
done

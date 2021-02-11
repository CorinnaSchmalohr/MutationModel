for t in luad breast skin colon ovary kidney prostate
do
(  echo $t
   mkdir -p data/procData/${t}/methylation
   for i in $(ls data/rawdata/methylation/${t}/ | grep '.bed.gz'); do
      echo $i
      id=$(echo $i | cut -d'.' -f1)
      if [ ! -f data/procData/${t}/methylation/${id}_sorted.bed ]
      then
         gunzip -k data/rawdata/methylation/${t}/${id}.bed.gz
         sort -k1,1 -k2,2n -S 8G --parallel 4 -T data/ \
            data/rawdata/methylation/${t}/${id}.bed \
            > data/procData/${t}/methylation/${id}_sorted.bed
         rm data/rawdata/methylation/${t}/${id}.bed
      fi
      # bedtools map -a data/procData/${t}/Muts.bed \
      #    -b data/procData/${t}/methylation/${id}_sorted.bed \
      #    -c 5,10,11 -o mean,mean,mean \
      #    > data/procData/${t}/methylation/${id}_local.out
      bedtools slop -i  data/procData/${t}/Muts.bed \
         -g data/rawdata/GRCh37.p11.genome.fa.fai -b 50 | \
         bedtools map -a - -b data/procData/${t}/methylation/${id}_sorted.bed \
         -c 5,10,11 -o mean,mean,mean \
         > data/procData/${t}/methylation/${id}_100bp.out
      # # bedtools slop -i data/procData/${t}/Muts.bed   \
      # #    -g data/rawdata/GRCh37.p11.genome.fa.fai -b 500000 | \
      # #    bedtools map -a - -b data/procData/${t}/methylation/${id}_sorted.bed \
      # #    -c 5,10,11 -o mean,mean,mean \
      # #    > data/procData/${t}/methylation/${id}_1Mbp.out
   done
) &
# 5: score 1-1000, 10:Number of reads or coverage, 11:Percentage of reads that show methylation at this position in the genome
done

# breast: column 10 and 11


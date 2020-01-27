for t in luad breast skin colon  kidney prostate #ovary
do
   echo ${t}
   mkdir data/procData/${t}/methylation_methbank
   case $t in
      luad)
         T=Lung
         ;;
      breast)
         T=Breast
         ;;
      skin)
         T=Skin
         ;;
      colon)
         T=Colon
         ;;
      kidney)
         T=Kidney
         ;;
      prostate)
         T=Prostate
         ;;
   esac
   echo $T
   gunzip data/rawdata/methylation_methbank/${t}/${T}.wig.gz
   convert2bed --input=wig --output=bed --do-not-sort \
         < data/rawdata/methylation_methbank/${t}/${T}.wig \
         > data/procData/${t}/methylation_methbank/${T}.bed
   lib/liftOver data/procData/${t}/methylation_methbank/${T}.bed  \
      data/rawdata/hg38ToHg19.over.chain \
      data/procData/${t}/methylation_methbank/${T}_liftover.bed  \
      data/procData/${t}/methylation_methbank/${T}_liftover.unmapped
   sort -k1,1 -k2,2n -S 8G --parallel 4 -T data/ \
      data/procData/${t}/methylation_methbank/${T}_liftover.bed \
      > data/procData/${t}/methylation_methbank/${T}_liftover_sorted.bed
   bedtools map -a data/procData/${t}/Muts.bed -b \
      data/procData/${t}/methylation_methbank/${T}_liftover_sorted.bed \
      -c 5 -o mean > data/procData/${t}/methylation_methbank/${T}.out.bed
   bedtools slop -i  data/procData/${t}/Muts.bed \
      -g data/rawdata/GRCh37.p11.genome.fa.fai -b 500 | \
      bedtools map -a - -c 5 -o mean \
      -b data/procData/${t}/methylation_methbank/${T}_liftover_sorted.bed \
      > data/procData/${t}/methylation_methbank/${T}_1kbp.out.bed
   bedtools slop -i  data/procData/${t}/Muts.bed \
      -g data/rawdata/GRCh37.p11.genome.fa.fai -b 500000 | \
      bedtools map -a - -c 5 -o mean \
      -b data/procData/${t}/methylation_methbank/${T}_liftover_sorted.bed \
      > data/procData/${t}/methylation_methbank/${T}_1Mbp.out.bed
done

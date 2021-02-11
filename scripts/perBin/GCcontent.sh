for bin in bins10bp bins100bp bins1kb bins10kb bins100kb bins1Mb ; do
   echo ${bin}
   bedtools getfasta -fi data/rawdata/GRCh37.p11.genome.fa \
         -bed data/procData/bins/${bin}/${bin}.bed -tab |
         grep -o -n [GC] | cut -d':' -f1 |
         uniq -c > data/procData/bins/${bin}/GCcontent.out
   for t in luad breast skin colon ovary kidney prostate; do
      echo ${t}
      bedtools getfasta -fi data/rawdata/GRCh37.p11.genome.fa \
         -bed data/procData/bins/${bin}/${bin}_${t}_covered.bed -tab |
         grep -o -n [GC] | cut -d':' -f1 |
         uniq -c > data/procData/bins/${bin}/${t}_covered_GCcontent.out
   done
done

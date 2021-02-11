for bin in bins10kb bins100kb bins1Mb ; do
   echo ${bin}
   bedtools getfasta -fi data/rawdata/GRCh37.p11.genome.fa \
      -bed data/procData/${bin}_allChrs.bed -tab |
      grep -o -n [GC] | cut -d':' -f1 |
      uniq -c > data/procData/${bin}_GCcontent.out
done

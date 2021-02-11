for t in luad breast skin colon ovary kidney prostate
do
   echo $t
   for i in $(ls data/procData/${t}/GCcontent/ | grep '.bed$' ); do
      echo $i
      id=$(echo $i | cut -d'.' -f1)
      bedtools getfasta -fi data/rawdata/GRCh37.p11.genome.fa \
      -bed data/procData/${t}/GCcontent/${id}.bed -tab | \
      grep -o -n [GC] | cut -d':' -f1 | \
      uniq -c > data/procData/${t}/GCcontent/${id}.out
   done
done
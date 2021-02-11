for t in luad breast skin colon ovary kidney prostate; do
   echo ${t}
   for bin in bins10bp bins100bp bins1kb bins10kb bins100kb bins1Mb ; do
      echo ${bin}
      bedops --partition data/procData/${t}/covered_regions_withIDs.bed \
         data/procData/bins/${bin}/${bin}.bed | 
         bedops --element-of 1 - data/procData/${t}/covered_regions_withIDs.bed | 
         awk -v OFS='\t' '{print $0,NR}' \
         > data/procData/bins/${bin}/${bin}_${t}_covered.bed
   done
done
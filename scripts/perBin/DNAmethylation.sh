for t in luad breast skin colon ovary kidney prostate; do
   echo ${t}
   for bin in bins10bp bins100bp bins1kb bins10kb bins100kb bins1Mb ; do
      echo ${bin}
      mkdir data/procData/bins/${bin}/${t}/methylation/
      for i in $(ls data/rawdata/methylation/${t}/ | grep '.bed.gz'); do
         id=$(echo $i | cut -d'.' -f1)
         bedmap --echo --wmean --delim '\t' \
            data/procData/bins/${bin}/${bin}.bed \
            data/procData/${t}/methylation/${id}_sorted.bed \
            > data/procData/bins/${bin}/${t}/methylation/${id}.out.bed
         bedmap --echo --wmean --delim '\t' \
            data/procData/bins/${bin}/${bin}_${t}_covered.bed \
            data/procData/${t}/methylation/${id}_sorted.bed \
            > data/procData/bins/${bin}/${t}/methylation/covered_${id}.out.bed 
      done
   done  
done

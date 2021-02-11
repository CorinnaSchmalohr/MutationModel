for t in luad breast skin colon ovary kidney prostate; do
   echo ${t}
   for bin in bins10bp bins100bp bins1kb bins10kb bins100kb bins1Mb ; do
      echo ${bin}
      mkdir data/procData/bins/${bin}/${t}/DNAaccessibility/
      mkdir data/procData/bins/${bin}/${t}/UCSC_tracks/
      for i in $(ls data/rawdata/DNAaccessibility/${t}/ | grep bigWig); do
         id=$(echo $i | cut -d'.' -f1)
         genome=$(cat data/rawdata/DNAaccessibility/${t}/metadata.tsv | grep ^${id} | cut -f44)
         if [ $genome == "GRCh38" ]
         then
            # cat data/procData/${t}/DNAaccessibility/${id}_liftOver_sorted.bedGraph |
            #    awk -vOFS="\t" '{ print $1, $2, $3, ".", $4 }' \
            #    > data/procData/${t}/DNAaccessibility/${id}_liftOver_sorted_UCSC.bedGraph
            # bedmap --echo --wmean --delim '\t' \
            #    data/procData/bins/${bin}/${bin}.bed \
            #    data/procData/${t}/DNAaccessibility/${id}_liftOver_sorted_UCSC.bedGraph \
            #    > data/procData/bins/${bin}/${t}/DNAaccessibility/${id}.out.bed
            bedmap --echo --wmean --delim '\t' \
               data/procData/bins/${bin}/${bin}_${t}_covered.bed \
               data/procData/${t}/DNAaccessibility/${id}_liftOver_sorted_UCSC.bedGraph \
               > data/procData/bins/${bin}/${t}/DNAaccessibility/covered_${id}.out.bed 
         else
            # lib/bigWigAverageOverBed \
            #    -bedOut=data/procData/bins/${bin}/${t}/DNAaccessibility/${id}.out.bed \
            #    data/rawdata/DNAaccessibility/${t}/${id}.bigWig \
            #    data/procData/bins/${bin}/${bin}.bed \
            #    data/procData/temp.out.tab
            lib/bigWigAverageOverBed \
               -bedOut=data/procData/bins/${bin}/${t}/DNAaccessibility/covered_${id}.out.bed \
               data/rawdata/DNAaccessibility/${t}/${id}.bigWig \
               data/procData/bins/${bin}/${bin}_${t}_covered.bed \
               data/procData/temp.out.tab
         fi
      done
      # bedmap  --echo --wmean --delim '\t' \
      #    data/procData/bins/${bin}/${bin}.bed \
      #    <(cat data/rawdata/UCSC_tracks/wgEncodeRegDnaseClusteredV3.sorted.bed |
      #       awk -vOFS="\t" '{ print $1, $2, $3, ".", $4 }') \
      #    > data/procData/bins/${bin}/${t}/UCSC_tracks/wgEncodeRegDnaseClusteredV3.out
      bedmap  --echo --wmean --delim '\t' \
         data/procData/bins/${bin}/${bin}_${t}_covered.bed \
         <(cat data/rawdata/UCSC_tracks/wgEncodeRegDnaseClusteredV3.sorted.bed | 
            awk -vOFS="\t" '{ print $1, $2, $3, ".", $4 }') \
         > data/procData/bins/${bin}/${t}/UCSC_tracks/covered_wgEncodeRegDnaseClusteredV3.out
   done  
done

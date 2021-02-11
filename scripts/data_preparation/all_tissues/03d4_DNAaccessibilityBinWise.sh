for t in luad breast skin prostate ovary colon  kidney ; do
   echo ${t}
   for bin in bins10kb bins100kb bins1Mb ; do
      #bins10bp bins100bp bins1kb
      echo ${bin}
      # create missing folders
      if [ ! -d data/procData/${t}/DNAaccessibility/bins ]
      then
         mkdir data/procData/${t}/DNAaccessibility/bins
      fi
      if [ ! -d data/procData/${t}/UCSC_tracks/bins ]
      then
         mkdir data/procData/${t}/UCSC_tracks/bins
      fi
      # iterate through DNAaccessibility files and get values for bins
      for i in $(ls data/rawdata/DNAaccessibility/${t}/ | grep bigWig); do
         id=$(echo $i | cut -d'.' -f1)
         genome=$(cat data/rawdata/DNAaccessibility/${t}/metadata.tsv | grep ^${id} | cut -f44)
         echo $id $genome
         if [ $genome == "GRCh38" ]
         # if source file from GRCh38, we have to use liftOver
         then
            if [ ! -d data/procData/${t}/DNAaccessibility/${id}_liftOver_sorted_UCSC.bedGraph ]
            then
               cat data/procData/${t}/DNAaccessibility/${id}_liftOver_sorted.bedGraph |
                  awk -vOFS="\t" '{ print $1, $2, $3, ".", $4 }' \
                  > data/procData/${t}/DNAaccessibility/${id}_liftOver_sorted_UCSC.bedGraph
            fi
            bedmap --echo --wmean --delim '\t' \
               data/procData/${bin}_allChrs.bed \
               data/procData/${t}/DNAaccessibility/${id}_liftOver_sorted_UCSC.bedGraph \
               > data/procData/${t}/DNAaccessibility/bins/${bin}_${id}.out.bed
         else
            # it crashes on large files, for those we have to change to bedGraph first
            if [ $(stat -c%s data/rawdata/DNAaccessibility/${t}/${id}.bigWig) -gt 1000000000 ] && [ ${bin} == bins1Mb ]
            then
               lib/bigWigToBedGraph \
                  data/rawdata/DNAaccessibility/${t}/${id}.bigWig \
                  data/procData/${t}/DNAaccessibility/${id}.bedGraph
               cat data/procData/${t}/DNAaccessibility/${id}.bedGraph |
                  awk -vOFS="\t" '{ print $1, $2, $3, ".", $4 }' \
                  > data/procData/${t}/DNAaccessibility/${id}_UCSC.bedGraph
               bedmap --echo --wmean --delim '\t' \
                  data/procData/${bin}_allChrs.bed \
                  data/procData/${t}/DNAaccessibility/${id}_UCSC.bedGraph \
                  > data/procData/${t}/DNAaccessibility/bins/${bin}_${id}.out.bed
            else
               lib/bigWigAverageOverBed \
                  -bedOut=data/procData/${t}/DNAaccessibility/bins/${bin}_${id}.out.bed \
                  data/rawdata/DNAaccessibility/${t}/${id}.bigWig \
                  data/procData/${bin}_allChrs.bed \
                  data/procData/temp.out.tab
            fi
         fi
      done
      # DNAaccessibility from UCSC
      echo 'UCSC'
      bedmap  --echo --wmean --delim '\t' \
         data/procData/${bin}_allChrs.bed \
         <(cat data/rawdata/UCSC_tracks/wgEncodeRegDnaseClusteredV3.sorted.bed |
            awk -vOFS="\t" '{ print $1, $2, $3, ".", $4 }') \
         > data/procData/${t}/UCSC_tracks/bins/${bin}_wgEncodeRegDnaseClusteredV3.out.bed
   done  
   rm data/procData/temp.out.tab
done

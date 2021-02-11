for t in luad breast skin colon ovary kidney prostate
do
   echo $t
   for bin in bins10bp bins100bp bins1kb bins10kb bins100kb bins1Mb ; do
      echo ${bin}
      if [ ! -d data/procData/bins/${bin}/${t}/ ]
      then
         mkdir data/procData/bins/${bin}/${t}/
      fi
      if [ ! -d data/procData/bins/${bin}/${t}/DNAbinding/ ]
      then
      mkdir data/procData/bins/${bin}/${t}/DNAbinding/
      fi
      for i in $(ls data/rawdata/DNAbinding/${t}/ | grep 'bigWig'); do
         id=$(echo $i | cut -d'.' -f1)
         # lib/bigWigAverageOverBed \
         #    -bedOut=data/procData/bins/${bin}/${t}/DNAbinding/${id}.out.bed \
         #    data/rawdata/DNAbinding/${t}/${id}.bigWig \
         #    data/procData/bins/${bin}/${bin}.bed \
         #    data/procData/temp_DNAbinding.out.tab
         lib/bigWigAverageOverBed \
            -bedOut=data/procData/bins/${bin}/${t}/DNAbinding/covered_${id}.out.bed \
            data/rawdata/DNAbinding/${t}/${id}.bigWig \
            data/procData/bins/${bin}/${bin}_${t}_covered.bed \
            data/procData/temp_DNAbinding.out.tab
      done
      rm data/procData/temp_DNAbinding.out.tab
      mkdir data/procData/bins/${bin}/${t}/UCSC_tracks/
      # bedtools map -a data/procData/bins/${bin}/${bin}.bed \
      #    -b data/rawdata/UCSC_tracks/wgEncodeRegTfbsClusteredV3.bed.gz \
      #    -c 4 -o distinct \
      #    > data/procData/bins/${bin}/${t}/UCSC_tracks/wgEncodeRegTfbsClusteredV3.out
      bedtools map -a data/procData/bins/${bin}/${bin}_${t}_covered.bed \
         -b data/rawdata/UCSC_tracks/wgEncodeRegTfbsClusteredV3.bed.gz \
         -c 4 -o distinct \
         > data/procData/bins/${bin}/${t}/UCSC_tracks/covered_wgEncodeRegTfbsClusteredV3.out
      for i in $(ls data/rawdata/UCSC_tracks/ | grep 'wgEncodeBroadHistone'); do
         id=$(echo $i | cut -d'.' -f1)
         # lib/bigWigAverageOverBed \
         #    -bedOut=data/procData/bins/${bin}/${t}/UCSC_tracks/${id}.out.bed \
         #    data/rawdata/UCSC_tracks/${id}.bigWig \
         #    data/procData/bins/${bin}/${bin}.bed \
         #    data/procData/temp_DNAbinding.out.tab
         lib/bigWigAverageOverBed \
            -bedOut=data/procData/bins/${bin}/${t}/UCSC_tracks/covered_${id}.out.bed \
            data/rawdata/UCSC_tracks/${id}.bigWig \
            data/procData/bins/${bin}/${bin}_${t}_covered.bed \
            data/procData/temp_DNAbinding.out.tab
      done
      rm data/procData/temp_DNAbinding.out.tab
   done 
done
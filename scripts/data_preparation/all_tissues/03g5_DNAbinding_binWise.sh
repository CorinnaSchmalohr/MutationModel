for t in luad breast skin colon ovary kidney prostate
do
   echo $t
   for bin in bins10kb bins100kb bins1Mb ; do
      echo ${bin}
      # create folders
      if [ ! -d data/procData/${t}/DNAbinding/bins/ ]
      then
         mkdir data/procData/${t}/DNAbinding/bins/
      fi
      if [ ! -d data/procData/${t}/UCSC_tracks/bins/ ]
      then
         mkdir data/procData/${t}/UCSC_tracks/bins/
      fi
      # DNAbinding
      echo DNAbinding
      for i in $(ls data/rawdata/DNAbinding/${t}/ | grep 'bigWig'); do
         id=$(echo $i | cut -d'.' -f1)
         echo $id
         lib/bigWigAverageOverBed \
            -bedOut=data/procData/${t}/DNAbinding/bins/${bin}_${id}.out.bed \
            data/rawdata/DNAbinding/${t}/${id}.bigWig \
            data/procData/${bin}_allChrs.bed \
            data/procData/temp_DNAbinding.out.tab
      done
      #tfbs
      echo tfbs
      bedtools map -a data/procData/${bin}_allChrs.bed \
         -b data/rawdata/UCSC_tracks/wgEncodeRegTfbsClusteredV3.bed.gz \
         -c 4 -o distinct \
         > data/procData/${t}/UCSC_tracks/bins/${bin}_wgEncodeRegTfbsClusteredV3.out
      # UCSC histones
      echo UCSChistones
      for i in $(ls data/rawdata/UCSC_tracks/ | grep 'wgEncodeBroadHistone'); do
         id=$(echo $i | cut -d'.' -f1)
         lib/bigWigAverageOverBed \
            -bedOut=data/procData/${t}/UCSC_tracks/bins/${bin}_${id}.out.bed \
            data/rawdata/UCSC_tracks/${id}.bigWig \
            data/procData/${bin}_allChrs.bed \
            data/procData/temp_DNAbinding.out.tab
      done
      rm data/procData/temp_DNAbinding.out.tab
   done 
done

for t in luad breast skin colon ovary kidney prostate; do
   echo ${t}
   if [ ! -d data/procData/${t}/UCSC_tracks/bins/ ]
      then
         mkdir data/procData/${t}/UCSC_tracks/bins/
   fi
   for bin in bins10kb bins100kb bins1Mb ; do
      echo ${bin}
      lib/bigWigAverageOverBed \
         -bedOut=data/procData/${t}/UCSC_tracks/bins/${bin}_wgEncodeSydhNsomeGm12878Sig.out.bed \
         data/rawdata/UCSC_tracks/wgEncodeSydhNsomeGm12878Sig.bigWig \
         data/procData/${bin}_allChrs.bed \
         data/procData/temp_nucleosome.out.tab
      lib/bigWigAverageOverBed \
         -bedOut=data/procData/${t}/UCSC_tracks/bins/${bin}_wgEncodeSydhNsomeK562Sig.out.bed \
         data/rawdata/UCSC_tracks/wgEncodeSydhNsomeK562Sig.bigWig \
         data/procData/${bin}_allChrs.bed \
         data/procData/temp_nucleosome.out.tab   
   done  
done

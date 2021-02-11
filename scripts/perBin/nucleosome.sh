for t in luad breast skin colon ovary kidney prostate; do
   echo ${t}
   for bin in bins10bp bins100bp bins1kb bins10kb bins100kb bins1Mb ; do
      echo ${bin}
      lib/bigWigAverageOverBed \
         -bedOut=data/procData/bins/${bin}/${t}/wgEncodeSydhNsomeGm12878Sig.out.bed \
         data/rawdata/UCSC_tracks/wgEncodeSydhNsomeGm12878Sig.bigWig \
         data/procData/bins/${bin}/${bin}.bed \
         data/procData/temp_nucleosome.out.tab
      lib/bigWigAverageOverBed \
         -bedOut=data/procData/bins/${bin}/${t}/wgEncodeSydhNsomeK562Sig.out.bed \
         data/rawdata/UCSC_tracks/wgEncodeSydhNsomeK562Sig.bigWig \
         data/procData/bins/${bin}/${bin}.bed \
         data/procData/temp_nucleosome.out.tab   
      lib/bigWigAverageOverBed \
         -bedOut=data/procData/bins/${bin}/${t}/wgEncodeSydhNsomeGm12878Sig_covered.out.bed \
         data/rawdata/UCSC_tracks/wgEncodeSydhNsomeGm12878Sig.bigWig \
         data/procData/bins/${bin}/${bin}_${t}_covered.bed \
         data/procData/temp_nucleosome.out.tab
      lib/bigWigAverageOverBed \
         -bedOut=data/procData/bins/${bin}/${t}/wgEncodeSydhNsomeK562Sig_covered.out.bed \
         data/rawdata/UCSC_tracks/wgEncodeSydhNsomeK562Sig.bigWig \
         data/procData/bins/${bin}/${bin}_${t}_covered.bed \
         data/procData/temp_nucleosome.out.tab         
   done  
done

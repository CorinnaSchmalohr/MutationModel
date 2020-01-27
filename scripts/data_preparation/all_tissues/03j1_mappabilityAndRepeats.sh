for t in luad breast skin colon ovary kidney prostate
do
   echo $t
   for i in $(ls data/rawdata/UCSC_tracks | grep 'wgEncodeCrgMapabilityAlign'); do
      echo $i
      id=$(echo $i | cut -d'.' -f1)
      lib/bigWigAverageOverBed -bedOut=data/procData/${t}/UCSC_tracks/${id}.out.bed \
         data/rawdata/UCSC_tracks/${id}.bigWig \
         data/procData/${t}/MutsWithIDs.bed \
         data/procData/temp.out.tab
   done   
   echo 'repeatMasker'
   bedtools map -a data/procData/${t}/Muts.bed \
      -b data/rawdata/UCSC_tracks/hg19.repeatMasker.UCSC.bed \
      -c 4 -o distinct >  data/procData/${t}/UCSC_tracks/hg19.repeatMasker.UCSC.out
   
   echo 'Trf'
   bedtools map -a data/procData/${t}/Muts.bed \
      -b data/rawdata/UCSC_tracks/hg19.TandemRepeatsFinder_woutextrascaffolds.UCSC.bed \
      -c 4 -o distinct >  data/procData/${t}/UCSC_tracks/hg19.TandemRepeatsFinder.UCSC.out
done
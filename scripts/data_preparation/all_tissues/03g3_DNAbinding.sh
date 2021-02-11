for t in luad breast skin colon ovary kidney prostate
do
   mkdir data/procData/${t}/DNAbinding/
   echo 'foldchange'
   for i in $(ls data/rawdata/DNAbinding/${t}/ | grep 'bigWig'); do
      echo $i
      id=$(echo $i | cut -d'.' -f1)
      lib/bigWigAverageOverBed -bedOut=data/procData/${t}/DNAbinding/${id}.out.bed \
         data/rawdata/DNAbinding/${t}/${id}.bigWig \
         data/procData/${t}/MutsWithIDs.bed data/procData/temp_DNAbinding.out.tab
      lib/bigWigAverageOverBed -sampleAroundCenter=100\
         -bedOut=data/procData/${t}/DNAbinding/${id}_100bp.out.bed \
         data/rawdata/DNAbinding/${t}/${id}.bigWig \
         data/procData/${t}/MutsWithIDs.bed data/procData/temp_DNAbinding.out.tab
   done
   rm data/procData/temp_DNAbinding.out.tab

   echo 'Tfbs'
   mkdir data/procData/${t}/UCSC_tracks/
   bedtools map -a data/procData/${t}/Muts.bed \
      -b data/rawdata/UCSC_tracks/wgEncodeRegTfbsClusteredV3.bed.gz \
      -c 4 -o distinct >  data/procData/${t}/UCSC_tracks/wgEncodeRegTfbsClusteredV3.out
   bedtools slop -i  data/procData/${t}/Muts.bed \
      -g data/rawdata/GRCh37.p11.genome.fa.fai -b 50 | \
      bedtools map -a - -b data/rawdata/UCSC_tracks/wgEncodeRegTfbsClusteredV3.bed.gz \
      -c 4 -o distinct \
      > data/procData/${t}/UCSC_tracks/wgEncodeRegTfbsClusteredV3_100bp.out


   for i in $(ls data/rawdata/UCSC_tracks/ | grep 'wgEncodeBroadHistone'); do
      id=$(echo $i | cut -d'.' -f1)
      echo $id
      lib/bigWigAverageOverBed -bedOut=data/procData/${t}/UCSC_tracks/${id}.out.bed \
         data/rawdata/UCSC_tracks/${id}.bigWig data/procData/${t}/MutsWithIDs.bed \
         data/procData/temp_DNAbinding.out.tab
      lib/bigWigAverageOverBed -sampleAroundCenter=100 \
         -bedOut=data/procData/${t}/UCSC_tracks/${id}_100bp.out.bed \
         data/rawdata/UCSC_tracks/${id}.bigWig data/procData/${t}/MutsWithIDs.bed \
         data/procData/temp_DNAbinding.out.tab   
      
   done
   rm data/procData/temp_DNAbinding.out.tab
done
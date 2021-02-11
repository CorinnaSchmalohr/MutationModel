for t in luad breast skin colon ovary kidney prostate
do
   mkdir data/procData/${t}/conservation
   lib/bigWigAverageOverBed \
      -bedOut=data/procData/${t}/conservation/hg19.100way.phyloP100way.out.bed \
      data/rawdata/UCSC_tracks/hg19.100way.phyloP100way.bigWig \
      data/procData/${t}/MutsWithIDs.bed data/procData/temp.out.tab
   lib/bigWigAverageOverBed -sampleAroundCenter=100 \
      -bedOut=data/procData/${t}/conservation/hg19.100way.phyloP100way_100bp.out.bed \
      data/rawdata/UCSC_tracks/hg19.100way.phyloP100way.bigWig \
      data/procData/${t}/MutsWithIDs.bed data/procData/temp.out.tab
done
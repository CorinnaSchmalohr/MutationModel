for t in luad breast skin colon ovary kidney prostate
do
   lib/bigWigAverageOverBed \
      -bedOut=data/procData/${t}/UCSC_tracks/wgEncodeSydhNsomeGm12878Sig.out.bed \
      data/rawdata/UCSC_tracks/wgEncodeSydhNsomeGm12878Sig.bigWig \
      data/procData/${t}/MutsWithIDs.bed \
      data/procData/temp_nucleosome.out.tab
   lib/bigWigAverageOverBed \
      -bedOut=data/procData/${t}/UCSC_tracks/wgEncodeSydhNsomeK562Sig.out.bed \
      data/rawdata/UCSC_tracks/wgEncodeSydhNsomeK562Sig.bigWig \
      data/procData/${t}/MutsWithIDs.bed \
      data/procData/temp_nucleosome.out.tab
done
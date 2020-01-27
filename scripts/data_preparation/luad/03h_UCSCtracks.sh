mkdir data/procData/luad/UCSC_tracks
for i in $(ls data/rawdata/UCSC_tracks | grep '.bigWig'); do
echo $i
lib/bigWigAverageOverBed -bedOut=data/procData/luad/UCSC_tracks/${i}.out.bed data/rawdata/UCSC_tracks/${i} data/procData/luad/MutsWithIDs.bed data/procData/temp.out.tab
#/lib/bigWigAverageOverBed -bedOut=/cellnet/somaticMutations/data/Mutmodel/Lung/UCSC_tracks/${i}_within20bp.out.bed -sampleAroundCenter=10 ${i} /cellnet/somaticMutations/data/Mutmodel/Lung/LungMutsForUCSCtracks.bed /cellnet/somaticMutations/data/Mutmodel/Lung/UCSC_tracks/${i}_within20bp.out.tab
done

echo 'gc5base'
bedtools map -a data/procData/luad/Muts.bed -b data/rawdata/UCSC_tracks/hg19.gc5Base.UCSC.sorted.bed -c 4 -o distinct > data/procData/luad/UCSC_tracks/hg19.gc5Base.UCSC.out

echo 'repeatMasker'
bedtools map -a data/procData/luad/Muts.bed -b data/rawdata/UCSC_tracks/hg19.repeatMasker.UCSC.bed -c 4 -o distinct >  data/procData/luad/UCSC_tracks/hg19.repeatMasker.UCSC.out

echo 'Trf'
bedtools map -a data/procData/luad/Muts.bed -b data/rawdata/UCSC_tracks/hg19.TandemRepeatsFinder_woutextrascaffolds.UCSC.bed -c 4 -o distinct >  data/procData/luad/UCSC_tracks/hg19.TandemRepeatsFinder.UCSC.out

echo 'Tfbs'
bedtools map -a data/procData/luad/Muts.bed -b data/rawdata/UCSC_tracks/wgEncodeRegTfbsClusteredV3.bed.gz -c 4 -o distinct >  data/procData/luad/UCSC_tracks/wgEncodeRegTfbsClusteredV3.out

echo 'DNAse'
bedtools map -a data/procData/luad/Muts.bed -b data/rawdata/UCSC_tracks/wgEncodeRegDnaseClusteredV3.sorted.bed -c 4 -o distinct >  data/procData/luad/UCSC_tracks/wgEncodeRegDnaseClusteredV3.out


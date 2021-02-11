# prepare UCSC tracks #####
cat data/rawdata/UCSC_tracks/hg19.repeatMasker.UCSC.out | awk 'BEGIN{OFS="\t"};NR>3{print $5,$6-1,$7,$11}' | bedtools sort -i - > data/rawdata/UCSC_tracks/hg19.repeatMasker.UCSC.bed

lib/varStepToBedGraph.pl data/rawdata/UCSC_tracks/hg19.gc5Base.UCSC.wigVarStep | grep -v '_' > data/rawdata/UCSC_tracks/hg19.gc5Base.UCSC.bed 
sort -k1,1 -k2,2n data/rawdata/UCSC_tracks/hg19.gc5Base.UCSC.bed > data/rawdata/UCSC_tracks/hg19.gc5Base.UCSC.sorted.bed
cat data/rawdata/UCSC_tracks/hg19.TandemRepeatsFinder.UCSC.bed |  grep -v '_' > data/rawdata/UCSC_tracks/hg19.TandemRepeatsFinder_woutextrascaffolds.UCSC.bed
bedtools sort -i data/rawdata/UCSC_tracks/wgEncodeRegDnaseClusteredV3.bed.gz > data/rawdata/UCSC_tracks/wgEncodeRegDnaseClusteredV3.sorted.bed
#####
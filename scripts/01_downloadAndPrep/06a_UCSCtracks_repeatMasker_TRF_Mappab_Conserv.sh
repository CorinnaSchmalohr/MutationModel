# repeatMasker
wget -O data/rawdata/UCSC_tracks/hg19.repeatMasker.UCSC.out.gz \
   http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.out.gz 
gunzip data/rawdata/UCSC_tracks/hg19.repeatMasker.UCSC.out.gz
convert2bed  --input=rmsk --output=bed --do-not-sort \
   < data/rawdata/UCSC_tracks/hg19.repeatMasker.UCSC.out \
   > data/rawdata/UCSC_tracks/hg19.repeatMasker.UCSC.bed
# Tandem Repeats Finder
wget --directory-prefix=data/rawdata/UCSC_tracks/ \
   http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.trf.bed.gz
# mappability:    
wget --directory-prefix=data/rawdata/UCSC_tracks/ \
   http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeCrgMapabilityAlign100mer.bigWig
wget --directory-prefix=data/rawdata/UCSC_tracks/ \
   http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeCrgMapabilityAlign24mer.bigWig 
wget --directory-prefix=data/rawdata/UCSC_tracks/ \
   http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeCrgMapabilityAlign40mer.bigWig
# consensus excludable
wget --directory-prefix=data/rawdata/UCSC_tracks/ \
   http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeDacMapabilityConsensusExcludable.bed.gz 
# conservation
wget --directory-prefix=data/rawdata/UCSC_tracks/ \
   http://hgdownload.soe.ucsc.edu/goldenPath/hg19/phyloP100way/hg19.100way.phyloP100way.bw
# TFBS 
wget --directory-prefix=data/rawdata/UCSC_tracks/ \
   http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeRegTfbsClustered/wgEncodeRegTfbsClusteredWithCellsV3.bed.gz
# UCSC histones was removed


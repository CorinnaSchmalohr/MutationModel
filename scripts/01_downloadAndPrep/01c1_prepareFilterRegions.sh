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


# convert mappability files to regions with mappability=1 (unique regions)
lib/bigWigToBedGraph \
  data/rawdata/UCSC_tracks/wgEncodeCrgMapabilityAlign100mer.bigWig \
  temp/wgEncodeCrgMapabilityAlign100mer.bedGraph
cat temp/wgEncodeCrgMapabilityAlign100mer.bedGraph | awk '$4 == 1 {print $0}' > \
  data/processedData/MappabilityAlign100mer_regionsToKeep.bedGraph
rm temp/wgEncodeCrgMapabilityAlign100mer.bedGraph

lib/bigWigToBedGraph \
data/rawdata/UCSC_tracks/wgEncodeCrgMapabilityAlign24mer.bigWig \
temp/wgEncodeCrgMapabilityAlign24mer.bedGraph
cat temp/wgEncodeCrgMapabilityAlign24mer.bedGraph | awk '$4 == 1 {print $0}' > \
data/processedData/MappabilityAlign24mer_regionsToKeep.bedGraph
rm temp/wgEncodeCrgMapabilityAlign24mer.bedGraph

lib/bigWigToBedGraph \
data/rawdata/UCSC_tracks/wgEncodeCrgMapabilityAlign40mer.bigWig \
temp/wgEncodeCrgMapabilityAlign40mer.bedGraph
cat temp/wgEncodeCrgMapabilityAlign40mer.bedGraph | awk '$4 == 1 {print $0}' > \
data/processedData/MappabilityAlign40mer_regionsToKeep.bedGraph
rm temp/wgEncodeCrgMapabilityAlign40mer.bedGraph
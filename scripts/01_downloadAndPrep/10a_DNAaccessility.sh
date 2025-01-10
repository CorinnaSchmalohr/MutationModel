# Nucleosome from UCSC
# wget --directory-prefix=data/rawdata/UCSC_tracks/ \
#    http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeSydhNsome/wgEncodeSydhNsomeGm12878Sig.bigWig
# wget --directory-prefix=data/rawdata/UCSC_tracks/ \
#    http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeSydhNsome/wgEncodeSydhNsomeK562Sig.bigWig

# DNAse from UCSC (use? it is a summary of ENCODE)
# wget --directory-prefix=data/rawdata/UCSC_tracks/ \
#    http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeRegDnaseClustered/wgEncodeRegDnaseClusteredV3.bed.gz    

# ENCODE
for t in brain breast colon esophagus kidney liver lung ovary prostate skin
  do
  echo $t
  cd data/rawdata/DNAseENCODE/${t}
  head -n1 files.txt | tr -d '"' > metafile.txt
  wget -O metadata.tsv  -i metafile.txt
  cd ../../../..
done




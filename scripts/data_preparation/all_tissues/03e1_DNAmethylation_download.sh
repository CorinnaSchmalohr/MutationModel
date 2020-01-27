for t in luad breast skin colon ovary kidney prostate
do
   echo ${t}
   cd data/rawdata/methylation/${t}
   head -n1 files.txt > metafile.txt
   xargs -L 1 curl -O -L < metafile.txt
   cat metadata.tsv | grep "methylation state at CpG" | grep "bed bedMethyl" | cut -f43 > filesToDownload.txt
   xargs -L 1 curl -O -L < filesToDownload.txt
   cd ../../../..
done
for t in luad breast skin colon ovary kidney prostate
do
   cd data/rawdata/DNAaccessibility/${t}
   head -n1 files.txt > metafile.txt
   xargs -L 1 curl -O -L < metafile.txt
   cat metadata.tsv | grep -e 'fold change over control' \
   -e 'read-depth normalized signal' -e 'base overlap signal' |
   cut -f43 > filesToDownload.txt
   xargs -L 1 curl -O -L < filesToDownload.txt
   cd ../../../..
done

for t in  breast skin colon kidney 
do
   echo ${t}
   cd data/rawdata/HiC_Encode/${t}
   head -n1 files.txt > metafile.txt
   xargs -L 1 curl -O -L < metafile.txt
   cat metadata.tsv | 
      grep -e 'chromatin interactions' \
      -e 'topologically associated domains' \
      -e 'genome compartments' |
      grep  'released' | grep 'hg19' | 
      cut -f43 > filesToDownload.txt
   xargs -L 1 curl -O -L < filesToDownload.txt
   cd ../../../..
done

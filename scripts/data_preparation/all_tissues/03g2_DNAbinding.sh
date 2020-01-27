for t in luad breast skin colon ovary kidney prostate
do
   cd data/rawdata/DNAbinding/${t}
   head -n1 files.txt > metafile.txt
   xargs -L 1 curl -O -L < metafile.txt
   cd ../../../..
   Rscript --vanilla --verbose scripts/data_preparation/all_tissues/03g1_DNAbinding.R $t
   cd data/rawdata/DNAbinding/${t}
   xargs -L 1 curl -O -L < filesToDownload.txt
   cd ../../../..
done
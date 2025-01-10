for t in luad breast skin colon ovary kidney prostate liver brain esophagus
do
   echo $t
   cd data/rawdata/DNAbinding/${t}
   xargs -L 1 curl -O -J -L < filesToDownload.txt
   cd ../../../..
done

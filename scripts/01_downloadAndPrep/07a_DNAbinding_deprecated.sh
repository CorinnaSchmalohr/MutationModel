for t in luad breast skin colon ovary kidney prostate liver brain esophagus
do
   echo $t
   cd data/rawdata/DNAbinding/${t}
   head -n1 files.txt > metafile.txt
   xargs -L 1 curl -O -L -J < metafile.txt
   cd ../../../..
done

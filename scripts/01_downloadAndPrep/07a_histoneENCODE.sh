for t in brain breast colon esophagus kidney liver lung ovary prostate skin
do
   echo $t
   cd data/rawdata/histoneENCODE/${t}
   head -n1 files.txt | tr -d '"' > metafile.txt
   wget -O metadata.tsv  -i metafile.txt
   cd ../../../..
done

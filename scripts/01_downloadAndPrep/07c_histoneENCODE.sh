for t in brain breast colon esophagus kidney liver lung ovary prostate skin
do
   echo $t
   # cd data/rawdata/histoneENCODE/${t}
   # xargs -L 1 curl -O -J -L < filesToDownload.txt
   # cd ../../../..
   flist=$(cat data/rawdata/histoneENCODE/${t}/filesToDownload.txt)
   for f in ${flist}; do
      wget --directory-prefix=data/rawdata/histoneENCODE/${t} $f
   done
done

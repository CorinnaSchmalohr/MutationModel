for t in brain breast colon esophagus kidney liver lung ovary prostate skin 
do
   echo $t
   flist=$(cat data/rawdata/HiCENCODE/${t}/filesToDownload.txt)
   for f in ${flist}; do
      wget --directory-prefix=data/rawdata/HiCENCODE/${t} $f
   done
done


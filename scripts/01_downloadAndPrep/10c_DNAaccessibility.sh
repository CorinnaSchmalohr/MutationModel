for t in brain breast colon esophagus kidney liver lung ovary prostate
  do
  echo $t
  flist=$(cat data/rawdata/DNAseENCODE/${t}/filesToDownload.txt)
  for f in ${flist}; do
    wget --directory-prefix=data/rawdata/DNAseENCODE/${t} $f
  done
done

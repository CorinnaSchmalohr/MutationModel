for t in brain breast colon esophagus kidney liver lung ovary prostate
  do
  echo $t
  flist=$(cat data/rawdata/ATACseqENCODE/${t}/filesToDownload.txt)
  for f in ${flist}; do
    wget --directory-prefix=data/rawdata/ATACseqENCODE/${t} $f
  done
done

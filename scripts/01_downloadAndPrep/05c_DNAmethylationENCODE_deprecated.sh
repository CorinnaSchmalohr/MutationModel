flist=$(cat data/rawdata/methylationENCODE/files2Download.csv | cut -f2 -d,)
for f in ${flist}; do
  wget --directory-prefix=data/rawdata/methylationENCODE $f
done


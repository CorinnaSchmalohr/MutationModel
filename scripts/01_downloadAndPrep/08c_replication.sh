
cd data/rawdata/replication
xargs -L 1 curl -O -J -L < filesToDownload.txt
cd ../../..

mkdir data/procData/luad/DNAbinding/
echo 'foldchange'
for i in $(cat data/procData/lung/DNAbinding/fc_files.txt); do 
echo $i
lib/bigWigAverageOverBed -bedOut=data/procData/luad/DNAbinding/${i}.out.bed data/rawdata/DNAbinding/lung/${i}.bigWig data/procData/luad/MutsWithIDs.bed data/procData/temp.out.tab
done


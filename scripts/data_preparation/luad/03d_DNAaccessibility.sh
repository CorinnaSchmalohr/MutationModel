cat data/procData/luad/Muts.bed | awk 'BEGIN{OFS="\t"};{print $0,NR}' > data/procData/luad/MutsWithIDs.bed # bigwigaverageoverbed needs a "name" for each bed entry

mkdir data/procData/luad/DNAaccessibility
for i in $(ls data/rawdata/DNAaccessibility/luad/ | grep bigWig	); do
echo $i
lib/bigWigAverageOverBed -bedOut=data/procData/luad/DNAaccessibility/${i}.out.bed data/rawdata/DNAaccessibility/luad/${i} data/procData/luad/MutsWithIDs.bed data/procData/temp.out.tab
done


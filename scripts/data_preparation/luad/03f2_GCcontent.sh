
for i in $(ls data/procData/luad/GCcontent/ | grep '.bed$'); do
echo $i
bedtools getfasta -fi data/rawdata/GRCh37.p11.genome.fa -bed data/procData/luad/GCcontent/${i} -tab |  grep -o -n [GC] | cut -d':' -f1 | uniq -c > data/procData/luad/GCcontent/${i}.out
done

mkdir data/procData/luad/replicationdomain
for i in $(ls data/procData/lung/replicationdomain/ | grep scaled); do
echo $i
bedtools sort -i data/procData/lung/replicationdomain/${i} | bedtools map -a data/procData/luad/Muts.bed -b - -c 4 -o sum > data/procData/luad/replicationdomain/${i}.out
done

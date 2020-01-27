mkdir data/procData/luad/methylation
for i in $(ls data/rawdata/methylation/luad/ | grep '.bed.gz'); do
echo $i
bedtools sort -i data/rawdata/methylation/luad/${i} > data/procData/luad/methylation/${i}_sorted.bed
bedtools map -a data/procData/luad/Muts.bed -b data/procData/luad/methylation/${i}_sorted.bed -c 5,10,11 -o mean,mean,mean > data/procData/luad/methylation/${i}_local.out
bedtools slop -i data/procData/luad/methylation/${i}_sorted.bed -g data/rawdata/GRCh37.p11.genome.fa.fai -b 500 | bedtools map -a data/procData/luad/Muts.bed -b - -c 5,10,11 -o mean,mean,mean > data/procData/luad/methylation/${i}_1kbp.out
bedtools slop -i data/procData/luad/methylation/${i}_sorted.bed -g data/rawdata/GRCh37.p11.genome.fa.fai -b 500000 | bedtools map -a data/procData/luad/Muts.bed -b - -c 5,10,11 -o mean,mean,mean > data/procData/luad/methylation/${i}_1Mbp.out
done

# 5: score 1-1000, 10:Number of reads or coverage, 11:Percentage of reads that show methylation at this position in the genome

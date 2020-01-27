for t in luad breast skin colon ovary kidney prostate
do
   mkdir data/procData/${t}/replication
   
   bedtools sort -i data/procData/Koren_hg19lifted_withSlopes.bed |
      bedtools closest -a data/procData/${t}/Muts.bed \
      -t first -b - > data/procData/${t}/replication/Koren.out
done
mkdir data/procData/luad/replication

bedtools sort -i data/procData/AsymTools/per_base_territories_20kb.bed | bedtools map -a data/procData/luad/Muts.bed -b - -c 8,9,10,11,12,13 -o distinct,distinct,distinct,distinct,distinct,distinct > data/procData/luad/replication/Asymtools.out

bedtools sort -i data/procData/Koren_hg19lifted_withSlopes.bed | bedtools closest -a data/procData/luad/Muts.bed -t first -b  -   > data/procData/luad/replication/Koren.out


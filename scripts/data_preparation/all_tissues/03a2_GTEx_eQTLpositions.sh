for t in luad breast skin colon ovary prostate
do
   echo $t
   bedtools map -a data/procData/${t}/Muts.bed \
   -b data/procData/${t}/GTEx_eqtls.bed \
   -c 4,5 -o distinct,min > data/procData/${t}/GTEx_eqtls.out
# 5: score 1-1000, 10:Number of reads or coverage, 11:Percentage of reads that show methylation at this position in the genome
done
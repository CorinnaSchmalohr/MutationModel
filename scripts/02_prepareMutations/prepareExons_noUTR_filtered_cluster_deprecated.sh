#!/bin/bash -l
#SBATCH --cpus-per-task=1
#SBATCH --account=cschmal1
#SBATCH --mail-type=END,FAIL,REQUEUE
#SBATCH --mail-use=cschmal1@uni-koeln.de
#SBATCH --output=/cellnet/MutationModel/output/codExons_noUTR_filtered_o-%j
#SBATCH --error=/cellnet/MutationModel/output/codExons_noUTR_filtered_e-%j

# translate gtf into bed
cat data/procData/codExons_noUTR_filtered.gtf |
  awk 'BEGIN{OFS="\t"}; {$4=$4-3; $5=$5+2; print $1,$4,$5, $1" "$4" "$5" "$8}' | 
  bedtools getfasta -fi data/rawdata/GRCh37.p11.genome.fa \
-bed - \
-name -tab > data/procData/codExons_noUTR_filtered.withsequences.bed

echo 'done'

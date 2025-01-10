for tissue in luad brca skcm coadread ov kirc prad 
do
(  filelist=$(ls data/rawdata/coverage/${tissue}/beds/ | grep bed)
   foldertissue=$(echo -e 'luad luad\nbrca breast\nskcm skin\ncoadread colon\nov ovary\nkirc kidney\nprad prostate' | 
   grep $tissue | cut -d' ' -f2)
   echo $tissue
   for file in $filelist 
   do
      id=$(echo $file | cut -d'.' -f1)
      cat data/rawdata/coverage/${tissue}/beds/${file} | 
         sort -k 1,1 -k2,2n |
         bedtools merge |
         awk -v OFS='\t' '{print $0,"1"}' > data/rawdata/coverage/${tissue}/beds/${id}.bg
   done
   echo $tissue 'unionbedg'
   bglist=$(ls data/rawdata/coverage/${tissue}/beds/*.bg)
   bedtools unionbedg -i $bglist > data/procData/${foldertissue}/coverage.bed
   cat data/procData/${foldertissue}/coverage.bed | 
      awk -v OFS='\t' '{ for(i=4; i<=NF;i++) j+=$i; 
            if ($1 ~ /^[:0-9:]/ || $1 ~ /[X,Y]/ ) {$1="chr"$1}
            if (j > NF/2) {print $1,$2,$3};
            j=0 
      }' | bedtools merge \
      > data/procData/${foldertissue}/covered_regions.bed
   echo $tissue 'getfasta'
   cat data/procData/${foldertissue}/covered_regions.bed | 
      awk -v OFS='\t' '{$2=$2-2; $3=$3+2; print $1,$2,$3,$1" "$2" "$3}' |
      bedtools getfasta -fi data/rawdata/GRCh37.p11.genome.fa \
      -bed - \
      -name -tab > data/procData/${foldertissue}/covered_regions_withsequences.bed
) &
done
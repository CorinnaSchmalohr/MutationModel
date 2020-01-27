for t in breast prostate #luad skin colon ovary kidney 
do
(  mkdir data/procData/${t}/DNAaccessibility
   for i in $(ls data/rawdata/DNAaccessibility/${t}/ | grep bigWig); do
      id=$(echo $i | cut -d'.' -f1)
      genome=$(cat data/rawdata/DNAaccessibility/${t}/metadata.tsv | grep ^${id} | cut -f44)
      echo $id $genome
      if [ $genome == "GRCh38" ]
      then
         lib/bigWigToBedGraph data/rawdata/DNAaccessibility/${t}/${id}.bigWig \
            data/procData/${t}/DNAaccessibility/${id}.bedGraph
         lib/liftOver  data/procData/${t}/DNAaccessibility/${id}.bedGraph \
            data/rawdata/hg38ToHg19.over.chain \
            data/procData/${t}/DNAaccessibility/${id}_liftOver.bedGraph \
            data/procData/${t}/DNAaccessibility/${id}_liftOver.unmapped
         sort -k1,1 -k2,2n -S 8G --parallel 8 -T data/ \
            data/procData/${t}/DNAaccessibility/${id}_liftOver.bedGraph \
            > data/procData/${t}/DNAaccessibility/${id}_liftOver_sorted.bedGraph
         bedtools map -a data/procData/${t}/Muts.bed -b \
            data/procData/${t}/DNAaccessibility/${id}_liftOver_sorted.bedGraph \
            -c 4 -o mean > data/procData/${t}/DNAaccessibility/${id}.out.bed
         bedtools slop -i  data/procData/${t}/Muts.bed \
            -g data/rawdata/GRCh37.p11.genome.fa.fai -b 500 | \
            bedtools map -a - -c 4 -o mean \
            -b data/procData/${t}/DNAaccessibility/${id}_liftOver_sorted.bedGraph \
            > data/procData/${t}/DNAaccessibility/${id}_1kbp.out.bed 
      else
         lib/bigWigAverageOverBed -bedOut=data/procData/${t}/DNAaccessibility/${id}.out.bed \
            data/rawdata/DNAaccessibility/${t}/${id}.bigWig \
            data/procData/${t}/MutsWithIDs.bed data/procData/temp.out.tab
         lib/bigWigAverageOverBed -sampleAroundCenter=1000 \
            -bedOut=data/procData/${t}/DNAaccessibility/${id}_1kbp.out.bed \
            data/rawdata/DNAaccessibility/${t}/${id}.bigWig \
            data/procData/${t}/MutsWithIDs.bed data/procData/temp.out.tab
      fi
   done
   bedtools map -a data/procData/${t}/Muts.bed \
      -b data/rawdata/UCSC_tracks/wgEncodeRegDnaseClusteredV3.sorted.bed \
      -c 4 -o distinct >  data/procData/${t}/UCSC_tracks/wgEncodeRegDnaseClusteredV3.out
) &
done
# When given two separate files for each sex, we took the female version 
for t in lung breast_F skin_F colon kidney Prostate esophagus brain_F liver_F
do
   echo $t
   wget --directory-prefix=data/rawdata/methylation_methbank \
      ftp://download.big.ac.cn/methbank/Human/CRM/${t}/${t}.wig.gz
   gunzip data/rawdata/methylation_methbank/${t}.wig.gz
   convert2bed --input=wig --output=bed --do-not-sort \
         < data/rawdata/methylation_methbank/${t}.wig \
         > data/rawdata/methylation_methbank/${t}.bed
   lib/liftOver data/rawdata/methylation_methbank/${t}.bed  \
      data/rawdata/hg38ToHg19.over.chain \
      data/rawdata/methylation_methbank/${t}_liftover.bed  \
      data/rawdata/methylation_methbank/${t}_liftover.unmapped
done
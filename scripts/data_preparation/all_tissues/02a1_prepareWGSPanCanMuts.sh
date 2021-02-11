for t in luad breast skin colon ovary kidney prostate; do
   name=$(cat data/rawdata/pancan/tissue2CancerType.txt | grep ${t} | cut -d',' -f2)
   echo $name
   head -n1 data/rawdata/pancan/final_consensus_passonly.snv_mnv_indel.icgc.public.maf \
      > data/rawdata/pancan/pcawg_icgc_WGS_${t}.maf 
   cat data/rawdata/pancan/final_consensus_passonly.snv_mnv_indel.icgc.public.maf | 
      grep ${name} >> data/rawdata/pancan/pcawg_icgc_WGS_${t}.maf
done
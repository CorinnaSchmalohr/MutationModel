# TCGA Pancan  mutations
wget -O data/rawdata/pancan/mc3.v0.2.8.PUBLIC.maf.gz \
  http://api.gdc.cancer.gov/data/1c8cfe5f-e52d-41ba-94da-f15ea1337efc


# ICGC PCAWG
wget -O data/rawdata/pancan/final_consensus_passonly.snv_mnv_indel.icgc.public.maf.gz \
https://dcc.icgc.org/api/v1/download?fn=/PCAWG/consensus_snv_indel/final_consensus_passonly.snv_mnv_indel.icgc.public.maf.gz


# SomaMutDB V1.4 (January 2024)
# manually selected tissues from https://vijglab.einsteinmed.org/SomaMutDB/browse/ and 
# downloaded into data/rawdata/SomaMutDB/

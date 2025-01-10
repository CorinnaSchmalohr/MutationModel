function validate_url()
{
    wget -O/dev/null -q $1
    return $?
}

for t in brain breast colon esophagus kidney liver lung ovary skin
do
  echo $t
  for f in ".txt.gz"  "_F.txt.gz" "_M.txt.gz" "_850K.txt.gz" "_850K_F.txt.gz" "_850K_M.txt.gz" 
  do
    url=ftp://download.cncb.ac.cn/methbank/MethBank4.0/Human_Healthy_Consensus_Reference_Methylomes/${t}${f}
    echo $url
    if validate_url $url; then
      echo "downloading"
      wget --directory-prefix=data/rawdata/methbank $url
    else
      echo "skipping"
    fi
  done
done

# get data for prostate separately (somehow not available in methbank 4.0)
wget --directory-prefix=data/rawdata/methbank \
ftp://download.cncb.ac.cn/methbank/MethBank3.0/Human/CRM_update/From_EWAS_Data_202203116/prostate_850K_M/prostate_850K_M.txt

# get probe ID sheets
wget -O data/rawdata/methbank/MethylationEPIC_v2.0_Files.zip \
  https://support.illumina.com/content/dam/illumina-support/documents/downloads/productfiles/methylationepic/MethylationEPIC%20v2%20Files.zip
unzip -d data/rawdata/methbank/ data/rawdata/methbank/MethylationEPIC_v2.0_Files.zip
wget --directory-prefix=data/rawdata/methbank \
  https://webdata.illumina.com/downloads/productfiles/humanmethylation450/humanmethylation450_15017482_v1-2.csv
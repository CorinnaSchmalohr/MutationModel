cd data/rawdata
wget -O non-b_db_human_hg19_tsv.tar.gz \
   https://ncifrederick.cancer.gov/bids/ftp/actions/download/?resource=/bioinfo/static/nonb_dwnld/human_hg19/human_hg19.tsv.tar.gz
tar -xz --one-top-level -f non-b_db_human_hg19_tsv.tar.gz
cd ../..

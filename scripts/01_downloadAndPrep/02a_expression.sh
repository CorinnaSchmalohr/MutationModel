######### cancer expression  #####################
cd data/rawdata/cancer_expr
for t in brca coadread esca gbm kirc lgg lihc luad ov prad skcm ucec 
do
echo $t
wget http://download.cbioportal.org/${t}_tcga_pan_can_atlas_2018.tar.gz
mkdir ${t}_tcga_pan_can_atlas_2018
tar -xzf ${t}_tcga_pan_can_atlas_2018.tar.gz -C ${t}_tcga_pan_can_atlas_2018
done
cd ..
# cd data/rawdata/cancer_expr
# declare -a tissues="SKCM BRCA COAD READ KIRC OV PRAD LUAD LGG GBM LIHC UCEC ESCA"
# for t in $(echo $tissues)
# do
#    echo $t
#    if [ ! -f  TCGA-${t}.htseq_fpkm.tsv.gz ]
#    then
#       echo "downloading " $t
#       wget https://gdc-hub.s3.us-east-1.amazonaws.com/latest/TCGA-${t}.htseq_fpkm.tsv.gz
#    fi
#    # average all the samples because it is too big to open in R
#    zcat TCGA-${t}.htseq_fpkm.tsv.gz | 
#       awk 'BEGIN{ OFS="\t";
#                   print "Ensembl_ID", "meanExpression";}
#            (NR>1){ sum=0;
#                    for(i=2; i<=NF; i++){sum+=log($i+1)};
#                    sum/=NF; 
#                    print $1, sum }' \
#       > ${t}_averaged.tsv
# done
# # file used for RNAseq quantification
# wget -O gencode.v22.annotation.gtf.gz \
#    https://api.gdc.cancer.gov/data/b011ee3e-14d8-4a97-aed4-e0b10f6bbe82 




################ GTEx healthy tissues###############
# have to use v7 of GTEx since it is still hg19
wget https://storage.googleapis.com/gtex_analysis_v7/rna_seq_data/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_median_tpm.gct.gz
# definitions of gene boundaries used by v7 GTEx
wget https://storage.googleapis.com/gtex_analysis_v7/reference/gencode.v19.genes.v7.patched_contigs.gtf 

cd ../..


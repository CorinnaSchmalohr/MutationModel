cd data/rawdata/cancer_expr
wget http://download.cbioportal.org/skcm_tcga_pan_can_atlas_2018.tar.gz
wget http://download.cbioportal.org/brca_tcga_pan_can_atlas_2018.tar.gz
wget http://download.cbioportal.org/coadread_tcga_pan_can_atlas_2018.tar.gz
wget http://download.cbioportal.org/kirc_tcga_pan_can_atlas_2018.tar.gz
wget http://download.cbioportal.org/ov_tcga_pan_can_atlas_2018.tar.gz
wget http://download.cbioportal.org/prad_tcga_pan_can_atlas_2018.tar.gz
wget http://download.cbioportal.org/luad_tcga_pan_can_atlas_2018.tar.gz

mkdir skcm_tcga_pan_can_atlas_2018
mkdir brca_tcga_pan_can_atlas_2018
mkdir coadread_tcga_pan_can_atlas_2018
mkdir kirc_tcga_pan_can_atlas_2018
mkdir ov_tcga_pan_can_atlas_2018
mkdir prad_tcga_pan_can_atlas_2018
mkdir luad_tcga_pan_can_atlas_2018

tar -xzf skcm_tcga_pan_can_atlas_2018.tar.gz     -C skcm_tcga_pan_can_atlas_2018
tar -xzf brca_tcga_pan_can_atlas_2018.tar.gz     -C brca_tcga_pan_can_atlas_2018
tar -xzf coadread_tcga_pan_can_atlas_2018.tar.gz -C coadread_tcga_pan_can_atlas_2018
tar -xzf kirc_tcga_pan_can_atlas_2018.tar.gz     -C kirc_tcga_pan_can_atlas_2018
tar -xzf ov_tcga_pan_can_atlas_2018.tar.gz       -C ov_tcga_pan_can_atlas_2018
tar -xzf prad_tcga_pan_can_atlas_2018.tar.gz     -C prad_tcga_pan_can_atlas_2018
tar -xzf luad_tcga_pan_can_atlas_2018.tar.gz     -C luad_tcga_pan_can_atlas_2018
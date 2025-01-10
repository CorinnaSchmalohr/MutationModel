# this script downloads and prepares the reference genome and annotation 
# for GRCh37 (hg19) from Gencode


# download reference genome
wget --directory-prefix=data/rawdata/ \
   https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_17/GRCh37.p11.genome.fa.gz
zcat data/rawdata/GRCh37.p11.genome.fa.gz \
   > data/rawdata/GRCh37.p11.genome.fa
samtools faidx data/rawdata/GRCh37.p11.genome.fa
# get chr sizes for bigWig creation
wget --directory-prefix=data/rawdata/ \
   http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes
# download the gencode annotation gtf 
# gtfs are 1-based
wget --directory-prefix=data/rawdata/ \
   ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz
# Chain files for liftOver
wget --directory-prefix=data/rawdata/ https://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz
wget --directory-prefix=data/rawdata/ http://hgdownload.cse.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz
# get positions of centromeres and telomeres
wget --directory-prefix=data/rawdata/\
   https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/gap.txt.gz


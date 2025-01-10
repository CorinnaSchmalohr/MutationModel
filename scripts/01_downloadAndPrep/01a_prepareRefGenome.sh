# this script downloads and prepares the reference genome and annotation 
# for GRCh37 (hg19) from Gencode



# download reference genome
wget --directory-prefix=data/rawdata/genome/ \
   https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/GRCh37_mapping/GRCh37.primary_assembly.genome.fa.gz
zcat data/rawdata/genome/GRCh37.primary_assembly.genome.fa.gz \
   > data/rawdata/genome/GRCh37.primary_assembly.genome.fa
samtools faidx data/rawdata/genome/GRCh37.primary_assembly.genome.fa
# get chr sizes for bigWig creation
wget --directory-prefix=data/rawdata/genome/ \
   http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes
wget  --directory-prefix=data/rawdata/genome/ \
https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.chrom.sizes
# download the gencode annotation gtf 
# gtfs are 1-based
wget --directory-prefix=data/rawdata/genome/ \
    https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/GRCh37_mapping/gencode.v43lift37.basic.annotation.gtf.gz
   
# Chain files for liftOver
wget --directory-prefix=data/rawdata/genome/ https://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz
wget --directory-prefix=data/rawdata/genome/ http://hgdownload.cse.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz
# get positions of centromeres and telomeres
wget --directory-prefix=data/rawdata/\
   https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/gap.txt.gz

# get problematic regions
wget -O data/rawdata/genome/UCSC_unusualRegions.bigBed https://hgdownload.soe.ucsc.edu/gbdb/hg19/bbi/problematic/comments.bb
lib/bigBedToBed data/rawdata/genome/UCSC_unusualRegions.bigBed data/rawdata/genome/UCSC_unusualRegions.bed
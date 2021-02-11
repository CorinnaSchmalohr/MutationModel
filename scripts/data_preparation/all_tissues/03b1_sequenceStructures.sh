for t in luad breast skin colon ovary kidney prostate
do
   mkdir data/procData/${t}/structures
   # a phased repeats
   bedtools intersect -a data/procData/${t}/Muts.bed \
      -b data/rawdata/non-b_db_human_hg19/a-phased_repeats_combined.gff \
      -c > data/procData/${t}/structures/aPhasedRepeats.tsv
   bedtools slop -i data/rawdata/non-b_db_human_hg19/a-phased_repeats_combined.gff \
      -g data/rawdata/GRCh37.p11.genome.fa.fai -b 100 | \
      bedtools intersect -a data/procData/${t}/Muts.bed \
      -b - -c  \
      > data/procData/${t}/structures/aPhasedRepeats_100bp.tsv
   # direct repeats
   bedtools intersect -a data/procData/${t}/Muts.bed -b \
      data/rawdata/non-b_db_human_hg19/direct_repeats_combined.gff \
      -c > data/procData/${t}/structures/directRepeats.tsv
   bedtools slop -i data/rawdata/non-b_db_human_hg19/direct_repeats_combined.gff \
      -g data/rawdata/GRCh37.p11.genome.fa.fai -b 100 | \
      bedtools intersect -a data/procData/${t}/Muts.bed \
      -b - -c  \
      > data/procData/${t}/structures/directRepeats_100bp.tsv
   # g quadruplex
   # bedtools intersect -a data/procData/${t}/Muts.bed -b \
   #    data/rawdata/non-b_db_human_hg19/g-quadruplex_forming_repeats_combined.gff \
   #    -c > data/procData/${t}/structures/gQuadruplex.tsv
   # bedtools slop -i data/rawdata/non-b_db_human_hg19/g-quadruplex_forming_repeats_combined.gff \
   #    -g data/rawdata/GRCh37.p11.genome.fa.fai -b 100 | \
   #    bedtools intersect -a data/procData/${t}/Muts.bed \
   #    -b - -c  \
   #    > data/procData/${t}/structures/gQuadruplex_100bp.tsv
   # # inverted repeats
   # bedtools intersect -a data/procData/${t}/Muts.bed -b \
   #    data/rawdata/non-b_db_human_hg19/inverted_repeats_combined.gff \
   #    -c > data/procData/${t}/structures/invertedRepeats.tsv
   # bedtools slop -i data/rawdata/non-b_db_human_hg19/inverted_repeats_combined.gff \
   #    -g data/rawdata/GRCh37.p11.genome.fa.fai -b 100 | \
   #    bedtools intersect -a data/procData/${t}/Muts.bed \
   #    -b - -c  \
   #    > data/procData/${t}/structures/invertedRepeats_100bp.tsv
   # # mirror repeats
   # bedtools intersect -a data/procData/${t}/Muts.bed -b \
   #    data/rawdata/non-b_db_human_hg19/mirror_repeats_combined.gff \
   #    -c > data/procData/${t}/structures/mirrorRepeats.tsv
   # bedtools slop -i data/rawdata/non-b_db_human_hg19/mirror_repeats_combined.gff \
   #    -g data/rawdata/GRCh37.p11.genome.fa.fai -b 100 | \
   #    bedtools intersect -a data/procData/${t}/Muts.bed \
   #    -b - -c  \
   #    > data/procData/${t}/structures/mirrorRepeats_100bp.tsv
   # # str
   # bedtools intersect -a data/procData/${t}/Muts.bed -b \
   #    data/rawdata/non-b_db_human_hg19/short_tandem_repeats_combined.gff \
   #    -c > data/procData/${t}/structures/shortTandemRepeats.tsv
   # bedtools slop -i data/rawdata/non-b_db_human_hg19/short_tandem_repeats_combined.gff \
   #    -g data/rawdata/GRCh37.p11.genome.fa.fai -b 100 | \
   #    bedtools intersect -a data/procData/${t}/Muts.bed \
   #    -b - -c  \
   #    > data/procData/${t}/structures/shortTandemRepeats_100bp.tsv
   # # z-DNA
   # bedtools intersect -a data/procData/${t}/Muts.bed -b \
   #    data/rawdata/non-b_db_human_hg19/z-dna_motifs_combined.gff \
   #    -c > data/procData/${t}/structures/zDNAmotifs.tsv
   # bedtools slop -i data/rawdata/non-b_db_human_hg19/z-dna_motifs_combined.gff \
   #    -g data/rawdata/GRCh37.p11.genome.fa.fai -b 100 | \
   #    bedtools intersect -a data/procData/${t}/Muts.bed \
   #    -b - -c  \
   #    > data/procData/${t}/structures/zDNAmotifs_100bp.tsv
   done
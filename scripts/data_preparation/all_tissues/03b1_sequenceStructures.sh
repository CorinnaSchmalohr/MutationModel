for t in luad breast skin colon ovary kidney prostate
do
   mkdir data/procData/${t}/structures
   bedtools intersect -a data/procData/${t}/Muts.bed \
      -b data/rawdata/non-b_db_human_hg19/a-phased_repeats_combined.gff \
      -c > data/procData/${t}/structures/aPhasedRepeats.tsv
   bedtools intersect -a data/procData/${t}/Muts.bed -b \
      data/rawdata/non-b_db_human_hg19/direct_repeats_combined.gff \
      -c > data/procData/${t}/structures/directRepeats.tsv
   bedtools intersect -a data/procData/${t}/Muts.bed -b \
      data/rawdata/non-b_db_human_hg19/g-quadruplex_forming_repeats_combined.gff \
      -c > data/procData/${t}/structures/gQuadruplex.tsv
   bedtools intersect -a data/procData/${t}/Muts.bed -b \
      data/rawdata/non-b_db_human_hg19/inverted_repeats_combined.gff \
      -c > data/procData/${t}/structures/invertedRepeats.tsv
   bedtools intersect -a data/procData/${t}/Muts.bed -b \
      data/rawdata/non-b_db_human_hg19/mirror_repeats_combined.gff \
      -c > data/procData/${t}/structures/mirrorRepeats.tsv
   bedtools intersect -a data/procData/${t}/Muts.bed -b \
      data/rawdata/non-b_db_human_hg19/short_tandem_repeats_combined.gff \
      -c > data/procData/${t}/structures/shortTandemRepeats.tsv
   bedtools intersect -a data/procData/${t}/Muts.bed -b \
      data/rawdata/non-b_db_human_hg19/z-dna_motifs_combined.gff \
      -c > data/procData/${t}/structures/zDNAmotifs.tsv
done
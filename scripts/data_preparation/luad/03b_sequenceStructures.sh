mkdir data/procData/luad/structures
bedtools intersect -a data/procData/luad/Muts.bed -b data/rawdata/non-b_db_human_hg19/a-phased_repeats_combined.gff -c > data/procData/luad/structures/aPhasedRepeats.tsv
bedtools intersect -a data/procData/luad/Muts.bed -b data/rawdata/non-b_db_human_hg19/direct_repeats_combined.gff -c > data/procData/luad/structures/directRepeats.tsv
bedtools intersect -a data/procData/luad/Muts.bed -b data/rawdata/non-b_db_human_hg19/g-quadruplex_forming_repeats_combined.gff -c > data/procData/luad/structures/gQuadruplex.tsv
bedtools intersect -a data/procData/luad/Muts.bed -b data/rawdata/non-b_db_human_hg19/inverted_repeats_combined.gff -c > data/procData/luad/structures/invertedRepeats.tsv
bedtools intersect -a data/procData/luad/Muts.bed -b data/rawdata/non-b_db_human_hg19/mirror_repeats_combined.gff -c > data/procData/luad/structures/mirrorRepeats.tsv
bedtools intersect -a data/procData/luad/Muts.bed -b data/rawdata/non-b_db_human_hg19/short_tandem_repeats_combined.gff -c > data/procData/luad/structures/shortTandemRepeats.tsv
bedtools intersect -a data/procData/luad/Muts.bed -b data/rawdata/non-b_db_human_hg19/z-dna_motifs_combined.gff -c > data/procData/luad/structures/zDNAmotifs.tsv

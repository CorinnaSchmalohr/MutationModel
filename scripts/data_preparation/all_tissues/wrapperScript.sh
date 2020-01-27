scripts/data_preparation/all_tissues/01_createFolders.sh

echo 'prepare Muts'
scripts/data_preparation/all_tissues/02a_preparePanCanCoverage.sh
Rscript --vanilla --verbose  \
   scripts/data_preparation/all_tissues/02a1_preparePanCanMuts.R
for n in $(seq 1 7); do Rscript --vanilla --verbose  \
   scripts/data_preparation/all_tissues/02a2_preparePanCanMuts.R $n; done
for n in $(seq 1 7); do Rscript --vanilla --verbose  \
   scripts/data_preparation/all_tissues/02a3_preparePanCanMuts.R $n; done
   
echo 'GTEx eQTL positions'
for n in $(seq 1 6); do Rscript --vanilla --verbose  \
   scripts/data_preparation/all_tissues/03a1_GTEx_eQTLpositions.R $n; done
scripts/data_preparation/all_tissues/03a2_GTEx_eQTLpositions.sh
for n in $(seq 1 7); do Rscript --vanilla --verbose  \
   scripts/data_preparation/all_tissues/03a3_GTEx_eQTLpositions.R $n; done

echo 'sequence structures'
scripts/data_preparation/all_tissues/03b1_sequenceStructures.sh
for n in $(seq 1 7); do Rscript --vanilla --verbose  \
   scripts/data_preparation/all_tissues/03b2_sequenceStructures.R $n; done

echo 'DNA accessibility'
scripts/data_preparation/all_tissues/03d1_DNAaccessibility.sh
scripts/data_preparation/all_tissues/03d2_DNAaccessibility.sh
for n in $(seq 1 7); do Rscript --vanilla --verbose  \
   scripts/data_preparation/all_tissues/03d3_DNAaccessibility.R $n; done

echo 'methylation'
scripts/data_preparation/all_tissues/03e1_DNAmethylation_download.sh
scripts/data_preparation/all_tissues/03e1_DNAmethylation.sh
scripts/data_preparation/all_tissues/03e1_DNAmethylation_methbank.sh
for n in $(seq 1 7); do Rscript --vanilla --verbose \
    scripts/data_preparation/all_tissues/03e2_DNAmethylation.R $n; done

echo 'GC content'
for n in $(seq 1 7); do Rscript --vanilla --verbose \
   scripts/data_preparation/all_tissues/03f1_GCcontent.R $n; done
scripts/data_preparation/all_tissues/03f2_GCcontent.sh
for n in $(seq 1 7); do Rscript --vanilla --verbose \
   scripts/data_preparation/all_tissues/03f3_GCcontent.R $n; done

echo 'exon/intron'
for n in $(seq 1 7); do Rscript --vanilla --verbose \
   scripts/data_preparation/all_tissues/03k_inexon.R $n; done
   
echo 'expression'
scripts/data_preparation/all_tissues/03c1_expression.sh
Rscript --vanilla --verbose \
   scripts/data_preparation/all_tissues/03c2_expression.R
   
echo 'DNA binding'
scripts/data_preparation/all_tissues/03g2_DNAbinding.sh
scripts/data_preparation/all_tissues/03g3_DNAbinding.sh
for n in $(seq 1 7); do Rscript --vanilla --verbose \
   scripts/data_preparation/all_tissues/03g4_DNAbinding.R $n; done

echo 'replication Direction and timing'
Rscript --vanilla --verbose \
   scripts/data_preparation/all_tissues/03i1_replication.R
scripts/data_preparation/all_tissues/03i2_replication.sh
Rscript --vanilla --verbose \
   scripts/data_preparation/all_tissues/03i3_replication.R
scripts/data_preparation/all_tissues/03i4_replication.sh
for n in $(seq 1 7); do Rscript --vanilla --verbose \
   scripts/data_preparation/all_tissues/03i5_replication.R $n; done

   
echo 'conservation'
scripts/data_preparation/all_tissues/03h1_conservation.sh
for n in $(seq 1 7); do Rscript --vanilla --verbose \
   scripts/data_preparation/all_tissues/03h2_conservation.R $n; done
   
   
echo 'mappability and repeats'
scripts/data_preparation/all_tissues/03j1_mappabilityAndRepeats.sh
for n in $(seq 1 7); do Rscript --vanilla --verbose \
   scripts/data_preparation/all_tissues/03j2_mappabilityAndRepeats.R $n; done

   
echo 'nucleosome'
scripts/data_preparation/all_tissues/03l1_nucleosome.sh
for n in $(seq 1 7); do Rscript --vanilla --verbose \
   scripts/data_preparation/all_tissues/03l2_nucleosome.R $n; done

# echo 'HiC'
scripts/data_preparation/all_tissues/03m1_HiC.sh
for n in $(seq 1 7); do Rscript --vanilla --verbose \
   scripts/data_preparation/all_tissues/03m2_HiC.R $n; done

echo '5mers'
scripts/data_preparation/all_tissues/02b1_get5mers.sh
for n in $(seq 1 7); do Rscript --vanilla --verbose  \
   scripts/data_preparation/all_tissues/02b2_get5mers.R $n; done

echo 'prepare dataframe'
for n in $(seq 1 7); do Rscript --vanilla --verbose \
   scripts/data_preparation/all_tissues/04_prepareDataframe.R $n; done
   
   
echo 'run RF'
for n in $(seq 1 7); do Rscript --vanilla --verbose \
   scripts/RFapplication/RF_application_alltissues.R $n; done
   
   
echo 'make summary plots'
for n in $(seq 1 7); do Rscript --vanilla --verbose \
   scripts/RFapplication/RF_application_alltissues_summaryStats.R $n; done
      
   
echo 'run RF per trimer'
for n in $(seq 1 7); do Rscript --vanilla --verbose \
   scripts/RFapplication/RF_application_alltissues_modelPerSignature.R $n; done
   
   
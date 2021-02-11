tissues = c("skin","ovary", "kidney", "prostate")
chrs = as.character(1:22)
for(tissue in tissues){
   print(tissue)
   muts = read.table(paste0("data/rawdata/pancan/pcawg_icgc_WGS_",
                            tissue, ".maf"),
                     sep = "\t", header = T, as.is = T, quote="")
   if(nrow(muts)<1){
      print(tissue)
      next
   }
   muts = muts[muts$Variant_Type == "SNP" & 
                  muts$Chromosome %in% chrs,]
   mafs = as.numeric(muts$i_1000genomes_AF)
   muts = muts[mafs<0.01 | is.na(mafs),]
   muts = muts[muts$i_snv_near_indel == "False",]
   muts = muts[muts$t_alt_count+muts$t_ref_count>14,]
   muts$Chromosome = paste0("chr", muts$Chromosome)
   # positions = paste(muts$Chromosome, muts$Start_position, sep = "_")
   # multiple = duplicated(positions) | duplicated(positions, fromLast=T)
   # submuts = submuts[!multiple,]
   muts = muts[, c("Chromosome", "Start_position", "End_position",
                           "Reference_Allele",  "Tumor_Seq_Allele2", "ref_context"),]
   muts$ref_context = substr(muts$ref_context, 9,13)
   
   save(muts, file = paste0("data/rdata/pancanWGS_icgc_muts_", 
                            tissue, ".RData"))
}

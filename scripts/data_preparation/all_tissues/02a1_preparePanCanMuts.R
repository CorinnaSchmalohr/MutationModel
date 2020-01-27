tissue2Cancer = list("luad" = "Lung adenocarcinoma",
                     "breast" = "Breast invasive carcinoma",
                     "skin" = "Skin Cutaneous Melanoma",
                     "colon" = c("Colon adenocarcinoma", "Rectum adenocarcinoma"),
                     "ovary" = "Ovarian serous cystadenocarcinoma",
                     "kidney" = "Kidney renal clear cell carcinoma",
                     "prostate" = "Prostate adenocarcinoma")
muts = read.table("data/rawdata/pancan/mc3.v0.2.8.PUBLIC.maf.gz",
                  header = T, na.strings=".", quote="", as.is=T)
cancerIDs = read.table("data/rawdata/pancan/TCGA-CDR-SupplementalTableS1_TSScode_info.csv", sep = ",", header = T, as.is = T, skip=2, na.strings="")
cancerID2cancer = cancerIDs$Disease.type
names(cancerID2cancer) = cancerIDs$TSS.Code
cancerCode = substr(muts$Tumor_Sample_Barcode, 6,7)
muts$cancerType = cancerID2cancer[cancerCode]
muts = muts[muts$cancerType %in% unlist(tissue2Cancer) & muts$Variant_Type == "SNP",]
muts = muts[,sapply(muts, function(x){!all(is.na(x))})]
muts$Chromosome = paste0("chr", muts$Chromosome)
chrs = paste0("chr", c(1:22))
muts = muts[muts$Chromosome %in% chrs,]
save(muts, file = "data/rdata/TCGA_muts.RData")
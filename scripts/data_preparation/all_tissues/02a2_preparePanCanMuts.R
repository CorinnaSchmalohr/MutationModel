# give 1 to 7 to get the different tissues
args = as.numeric(commandArgs(trailingOnly=T))
tissues = c("luad", "breast", "skin", "colon", "ovary", "kidney", "prostate")
tissue = tissues[args]
print(tissue)

setwd("/cellnet/MutationModel")

tissue2Cancer = list("luad" = "Lung adenocarcinoma",
                     "breast" = "Breast invasive carcinoma",
                     "skin" = "Skin Cutaneous Melanoma",
                     "colon" = c("Colon adenocarcinoma", "Rectum adenocarcinoma"),
                     "ovary" = "Ovarian serous cystadenocarcinoma",
                     "kidney" = "Kidney renal clear cell carcinoma",
                     "prostate" = "Prostate adenocarcinoma")
chrs = paste0("chr", c(1:22))

load("data/rdata/TCGA_muts.RData")
covered = read.table(paste0("data/procData/", tissue,
                            "/covered_regions.bed"),as.is=T)
colnames(covered) = c("chr", "start", "end")
covered = covered[covered$chr %in% chrs,]

library(parallel)


# prepare mutations #####
submuts = muts[muts$cancerType %in% tissue2Cancer[[tissue]],]
rm(muts)
save(submuts, file = paste0("data/rdata/", tissue, "/submuts.RData"))
# filter out variants with >1% population frequency
mafs = submuts[,c("GMAF", "AFR_MAF", "AMR_MAF", "EAS_MAF",
                  "EUR_MAF", "SAS_MAF", "AA_MAF", "EA_MAF")]
mafs = apply(mafs,2,function(x){
   sapply(x, function(y){
      temp = strsplit(y,split=",")[[1]]
      mean(as.numeric(substr(temp, 3,20)))
   })
})
mafsToRemove = apply(mafs,1,function(x){!sum(x > 0.01,na.rm=T)>0})
submuts = submuts[mafsToRemove,]
# filter out variants with to little depth
submuts = submuts[submuts$n_depth>=8 & submuts$t_depth >= 14,]
# these cutoffs are based on what was used for the "sufficient coverage" files.
# exclude positions that were mutated more than once, to exclude possible selection
positions = paste(submuts$Chromosome, submuts$Start_Position, sep = "_")
multiple = duplicated(positions) | duplicated(positions, fromLast=T)
submuts = submuts[!multiple,]
exomemuts = submuts[, c("Chromosome", "Start_Position", "End_Position",
                        "Reference_Allele",  "Tumor_Seq_Allele2", "CONTEXT"),]
exomemuts$CONTEXT = substr(exomemuts$CONTEXT, 4,8)
print("filtering for coverage")
cl = makeCluster(8,type = "SOCK")
filtered = parApply(cl,exomemuts,1,function(x, covered = covered){
   cr = x["Chromosome"]
   p = as.integer(x["Start_Position"])
   temp = covered[covered$chr == cr,]
   any(temp$start < p & temp$end >= p)
}, covered = covered)
stopCluster(cl)
print(table(filtered))
exomemuts = exomemuts[filtered,]
save(exomemuts, file = paste0("data/rdata/", tissue, "/exomemuts.RData"))
#####



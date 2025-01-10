setwd("/cellnet/MutationModel/")

library(vcfR)
library(stringr)

tissues = c("luad",  "skin", "colon", "ovary", "kidney", "prostate", "breast")
t = "breast"

bases = c("A", "C", "G", "T")


load(paste0("/cellnet/MutationModel/data/rdata/", t,"/Muts.RData"))
  
# Generate VCF File #####
input_file <- lapply(1:nrow(Muts), function(i){
  CHROM = as.numeric(strsplit(x = Muts$chr[i], split = "chr")[[1]][2])
  POS = Muts$pos[i]
  REF = Muts$ref[i]
  ID = paste0(CHROM,"_",POS)
  ALT = Muts$alt[i]
  QUAL = 100
  FILTER = "PASS"
  INFO = paste0("ST=", Muts$strand[i])
  
  if(is.na(ALT)){
    alleles = bases[-which(bases == REF)] # all alternativ alleles 
    
    chr3 = rep(x = CHROM, 3)
    pos3 = rep(x = POS, 3)
    ref3 = rep(x = REF, 3)
    id3 = rep(x = ID, 3)
    qual3 = rep(x = QUAL, 3)
    filter3 = rep(x = FILTER, 3)
    info3 = rep(x = INFO, 3)
    
    return(cbind(CHROM = chr3, POS = pos3, ID = id3, REF = ref3, 
                 ALT = alleles, QUAL = qual3, FILTER = filter3, INFO = info3))
  } else {
    return(cbind(CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO))
  }
})
  
vcffile = as.data.frame(do.call(rbind, input_file))
vcffile = vcffile[order(as.numeric(vcffile[,1]), as.numeric(vcffile[,2])),] # Sort by chrom and position
  
colnames(vcffile)[1] = "#CHROM"
  
# Write Meta-info lines
fileConn<-file(paste0("./data/predictors/SIFT/SIFT_input_",t,".vcf"))
writeLines(c("##fileformat=VCFv4.0", 
             paste0("##source=/cellnet/MutationModel/data/rdata/",t,"/Muts.RData"), 
             "##reference=./lib/GRCh37.74/", 
             "##INFO=<ID=ST,Number=1,Type=String,Description=Strand>"), fileConn)
close(fileConn)
  
# Write VCF file columns 
write.table(vcffile, file = paste0("data/predictors/SIFT/SIFT_input_",t,".vcf"), append = T, 
            sep = "\t", quote = F, col.names = T, row.names = F)
######

# Function to extract Variant Type from SIFT output and return as score  #####
FUN_VarType = function(i){
  VarType = strsplit(x = i, split = "|", fixed = T)[[1]][6]
  
  if(VarType == "NONCODING" | is.na(VarType)){ # If NA or Noncoding
    VarType_Score = 0
    
  } else if(VarType == "SYNONYMOUS" & !is.na(VarType)){ # If synonymous
    VarType_Score = 1
    
  } else if(VarType == "NONSYNONYMOUS" & !is.na(VarType)){ # If nonsynonymous
    VarType_Score = 2
    
  } else { # If START-LOST, STOP-GAIN or STOP-LOSS
    VarType_Score = 3
  }
  
  return(VarType_Score)
}
#####


# Run SIFT in the console

cmd = paste0("java -jar ./lib/SIFT4G_Annotator.jar -c -i ./data/predictors/SIFT/SIFT_input_",t,".vcf -d ./lib/GRCh37.74/ -r ./data/predictors/SIFT/")
system(command = cmd) # Run SIFT
  
# Read results
res = read.vcfR(file = paste0("./data/predictors/SIFT/Sift_input_",t,"_SIFTpredictions.vcf")) # Open SIFT results 
res = as.data.frame(res@fix)
    
# Get Variant Type
res_wVarType = data.frame(res, VarType = apply(X = res[8], MARGIN = 1, FUN = FUN_VarType)) # Turn Variant type output in a score 
  
ID_multALT = unique(res_wVarType[(duplicated(res_wVarType$ID)),]$ID) # IDs of all Positions where all alt. alleles were checked
VarType = sapply(ID_multALT, function(id){ # Calculates the mean VarType score of each alternative allel set !Dauert sehr lange!
  i = res_wVarType[res_wVarType$ID == id,]
  return(c(as.numeric(i$POS[1]),as.numeric(mean(i$VarType))))
}) %>% t()

res_singleALT = res_wVarType[!res_wVarType$ID %in% ID_multALT, ]
res_singleALT = t(apply(res_singleALT[c(2,9)], 1, as.numeric)) # Position and VarType score of all pos with known alternativ allel

VarType = rbind(res_singleALT, VarType) # Combine all extracted scored 
VarType = VarType[order(as.numeric(VarType[,1]), decreasing = F),] # Sort by position 
colnames(VarType) = c("pos", "VarType")

save(VarType, file = paste0("./data/predictors/SIFT/SIFT_VarType_",t,".RData"))
  
# give 1 to 7 to get the different tissues
args = as.numeric(commandArgs(trailingOnly=T))

##### Select tissue #####
tissues = c("luad",  "skin", "colon", "ovary", "kidney", "prostate", "breast")
tissue = tissues[args]
#########################

chroms = paste0("chr", 1:22)

path = paste0("/cellnet/MutationModel/data/rdata/", tissue, "/logRegression")


# Drop1 command 
pValues <- sapply(chroms, function(chrom){
  
  load(paste0(path, "/", chrom,"/logR_", chrom,".RData")) 
  
  drop1_features <- drop1(object = logR, test = "LRT")$`Pr(>Chi)`
  #save(drop1_features, file =  paste0(path, "/", chrom,"/pValues_features_glm.RData"))
  return(drop1_features)
})

pValues <- pValues[-1,] # remove first empty row
pValues <- as.data.frame(pValues)
save(pValues, file = paste0(path, "/pValues_allChr.RData"))

setwd("/cellnet/MutationModel/")
### Packages ###

###############

##### Select tissue #####
tissues = c("luad",  "skin", "colon", "ovary", "kidney", "prostate", "breast")
#tissue = tissues[1]
#########################

set.seed(1)

lapply(tissues, function(tissue){
  print(tissue)

  ##### Load preprocessed data with bins and create paths for data #####
  dir.create(paste0("data/rdata/", tissue, "/logRegression"), showWarnings=F)
  path = paste0("data/rdata/", tissue, "/logRegression")
  
  load(paste0("data/rdata/", tissue, "/completeData_withBinwise.RData"))
  
  data$context = NULL
  data$pentamer = NULL
  data$trimer = NULL
  data$septamer = NULL
  data$inexon = NULL

  data_wChr <- cbind(removed$chr, data) # Dataset with corresponding chromosome assignment
  colnames(data_wChr)[1] <- "chr"
  
  chroms <- unique(removed$chr) # List of all chromosomes
  
  
  ##### Create glm with leave-one-chromosome-out CV #####
  lapply(chroms, function(chrom){
    print(chrom)
    dir.create(paste0(path, "/", chrom ), 
               showWarnings=F) # Folder for each chromosom
    
    
    # Data for training; n-1 chromosomes without the column chr 
    data_train <- subset(data_wChr, chr != chrom, select = -c(chr))
    
    # Create log glm; culumn mutated VS all the others
    logR <- glm(formula = mutated ~ ., data = data_train, 
                family = binomial(link = "logit"))
    save(logR, file = paste0(path, "/", chrom,"/logR_", chrom,".RData"))
  })
  
  
  
  ##### Test model on remaining chromosome --> Leave-one-chromosome CV #####
  lapply(chroms, function(chrom){
    print(chrom)
    load(paste0(path, "/", chrom,"/logR_", chrom,".RData")) # Load glm with the current chrom not trained on
    
    # Data for testing; current chromosome data without the columnes chr and mutated 
    data_test <- subset(data_wChr, chr == chrom, select = -c(chr, mutated))
    
    # Prediction on the current chromosome with the glm where this chromosome is not trained
    yhat <- predict(logR, newdata = data_test, type = "response") # type=respones to get the actual predicted values between 0 and 1
    save(yhat, file = paste0(path, "/", chrom,"/yhat_",chrom,".RData"))
  })
  
})  

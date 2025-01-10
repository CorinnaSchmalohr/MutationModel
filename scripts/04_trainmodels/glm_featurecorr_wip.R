source("scripts/05_analysis/00_NamesAndColors.R")
library(tidyr)
library(tibble)
predictorOrder = read.table("scripts/04_trainmodels/predictorOrder.txt")[,1]
rarePredictors = c("ZBTB33_100bp", "YY1_100bp", "TAF1_100bp", "SP1_100bp",
                   "RXRA_100bp", "REST_100bp", "RAD21_100bp",
                   "NR2F2_100bp", "MAX_100bp", "JUND_100bp", "HNF4G_100bp",
                   "HNF4A_100bp", "GABPA_100bp", "FOXA2_100bp", 
                   "FOXA1_100bp", "EGR1_100bp", "ATF3_100bp")

# load
load("data/rdata/RFmodel/predPerTissueRF.RData")
load("data/rdata/GLMmodel/predPerTissueGLM.RData")
load("data/rdata/dataInfos.RData")
load("data/rdata/RFmodel/ROC_PR_RF_perChr.RData")
load("data/rdata/RFmodel/ROC_PR_RF_concat.RData")
load("data/rdata/GLMmodel/ROC_PR_glm_perChr.RData")
load("data/rdata/GLMmodel/ROC_PR_glm_concat.RData")
load("data/rdata/RFmodel/RF_imps.RData")
load("data/rdata/GLMmodel/GLM_impsAndPvals.RData")
load("data/rdata/nPerTissue.RData")
predictorOrder = rownames(rf_gini[[1]])
#####

pValue_features_allChr = glm_pvals[[tissue]]

## Adjust p-Values for multiple testing with FDRcorrection
glm_pvals_adj = sapply(tissues, function(tissue){
  pvals = glm_pvals[[tissue]]
  
  temp = sapply(1:ncol(pvals), function(cols){
    return(p.adjust(p = pvals[,cols], method = "fdr"))
  }) %>% as.data.frame()
  colnames(temp) <- colnames(pvals)
  return(temp)
}, simplify = F)


## 
pCutOff <- 0.01

set.seed(1)

lapply(tissues, function(tissue){
  print(tissue)
  
  ##### Load preprocessed data with bins and create paths for data and figures #####
  load(paste0(paste0("data/procData/traindata/traindata_processed_", tissue, ".RData")))
  chroms = unique(datchroms)
  
  pvals = glm_pvals_adj[[tissue]]
  
  # Subset of all features with a specific pValue on at least one chromosome
  print("Creating data...")
  pvals_wRown <- rownames_to_column(pvals)
  
  feature_list <- sapply(chroms, function(chrom){
    temp_subset <- subset(x = pvals_wRown[c("rowname", chrom)], subset = pvals_wRown[,chrom] <= pCutOff) # All features with a certain p-Value for each chr left out
    return(temp_subset["rowname"])
  }) %>% unlist() %>% as.data.frame()
  
  unique_features <- unique(feature_list) # Names of all features which showed up at least one time in the feature list
  subset_dat <- dat[unique_features] # Data set containing only the significant features
  subset_dat <- cbind(mutated = dat$mutated, subset_data_features_padj) # Dataset with corresponding chromosome and mutation assignme
  save(subset_dat, file = paste0(path, "/data_subset_pCutOff",pCutOff,".RData"))
  
  
  ##### Create glm with leave-one-chromosome-out CV and selected features#####
  print("Train and test model...")
  
  temp = lapply(chroms, function(cr){
    cat(cr, ' ')
    trainData = subset_dat[datchroms != cr,]
    logR = glm(formula = mutated ~ ., data = trainData, 
               family = binomial(link = "logit"))
    save(logR, file = paste0("data/rdata/GLMmodel/", tissue, "_", cr, "_sigFeatures.RData"))
    return(NULL)
  })
  cat('\n')
  #####
  
  # get variable importances ####
  print("var Importances")
  imp = sapply(chroms, function(cr){
    cat(cr, ' ')
    load(paste0("data/rdata/GLMmodel/", tissue, "_", cr, "_sigFeatures.RData"))
    return(logR$coefficients)
  })
  save(imp, file = paste0("data/rdata/GLMmodel/", tissue,
                          "_importances_sigFeatures.RData"))
  cat('\n')
  #####
})
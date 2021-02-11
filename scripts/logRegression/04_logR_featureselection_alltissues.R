setwd("/cellnet/MutationModel/")

### Packages ###
library(ggplot2)
library(tidyr)
library(tibble)
library(ggforce)
###############

##### Select tissue #####
tissues = c("luad",  "skin", "colon", "ovary", "kidney", "prostate", "breast")
tissue = tissues[1]
#########################


#replace pValues_features_allChr_glm with pValues_allChr and pVlaues with pValue_features_allChr when drop1 finished !!!!!!

##### VARIABLE IMPORTANCES #####
lapply(tissues, function(tissue){
  print(tissue)
  
  ##### Load preprocessed data with bins and create paths for data and figures #####
  path = paste0("data/rdata/", tissue, "/logRegression")
  
  load(paste0("data/rdata/", tissue, "/completeData_withBinwise.RData"))
  
  data$context = NULL
  data$pentamer = NULL
  data$trimer = NULL
  data$septamer = NULL
  data$inexon = NULL
  
  chroms = paste0("chr", seq(1,22)) # List of all chromosomes
  
  
  ##### Extract coefficiens and p-values ####
  # Get all calculated coefficients for each leave-one-chromosome-out CV
  print("Get coefficients...")
  
  coefficients <- sapply(chroms, function(chrom){
    
    load(paste0(path, "/", chrom,"/logR_", chrom,".RData"))
    
    return(coef(logR))
  }) %>% as.data.frame()
  
  save(coefficients, file = paste0(path, "/coefficients_allChr.RData"))
  
  
  # Get p-values for each model feature - Done by using the cluster with the script 04a_logR_cluster_drop1
  
  #pValues <- sapply(chroms, function(chrom){
  #  
  #  load(paste0(path, "/", chrom,"/logR_", chrom,".RData")) 
  #  
  #  drop1_features <- drop1(object = logR, test = "LRT")$`Pr(>Chi)` # Drop one feature after another and calculates how much the model changed 
  #  return(drop1_features)
  #  
  #})
#  
#  pValues <- pValues[-1,] # remove first empty row
#  rownames(pValues) <- colnames(data)[-c(2)] # add names without mutated
#  pValues <- as.data.frame(pValues)
#  save(pValues, file = paste0(path, "/pValues_allChr.RData"))
  

  # Adjust p-Values for multiple testing with FDRcorrection  
  load(file = paste0(path, "/pValues_features_allChr_glm.RData"))
  rownames(pValue_features_allChr) <- colnames(data)[-c(2)] # add names without mutated
  

  pValues_adjust <- sapply(1:ncol(pValue_features_allChr), function(cols){
    return(p.adjust(p = pValue_features_allChr[,cols], method = "fdr"))
  }) %>% as.data.frame()
  
  rownames(pValues_adjust) <- rownames(pValue_features_allChr)
  colnames(pValues_adjust) <- colnames(pValue_features_allChr)
  save(pValues_adjust, file = paste0(path, "/pValues_allChr_adjust.RData"))
  
})



##### REDUCE FEATRUES BASED ON SIGNIFICANCE ##### 
pCutOff <- 0.05

set.seed(1)

lapply(tissues, function(tissue){
  print(tissue)
  
  ##### Load preprocessed data with bins and create paths for data and figures #####
  path = paste0("data/rdata/", tissue, "/logRegression")
  
  load(paste0("data/rdata/", tissue, "/completeData_withBinwise.RData"))
  load(paste0(path, "/pValues_allChr_adjust.RData"))
  
  chroms = unique(removed$chr) # List of all chromosomes



  # Subset of all features with a specific pValue on at least one chromosome
  print("Creating data...")
  pValues_adjust_wRown <- rownames_to_column(as.data.frame(pValues_adjust))
  
  feature_list <- sapply(chroms, function(chrom){
    temp_subset <- subset(x = pValues_adjust_wRown[c("rowname", chrom)], subset = pValues_adjust_wRown[,chrom] <= pCutOff) # All features with a certain p-Value for each chr left out
    return(temp_subset["rowname"])
  }) %>% unlist() %>% as.data.frame()
  
  unique_features <- unique(feature_list) # Names of all features which showed up at least one time in the feature list
  
  subset_data_features_padj <- data[as.vector(unique_features[[1]])] # Data set containing only the significant features
  
  subset_data_features_padj <- cbind(data$mutated, removed$chr, subset_data_features_padj) # Dataset with corresponding chromosome and mutation assignme

  colnames(subset_data_features_padj)[c(1,2)] <- c("mutated", "chr")
  save(subset_data_features_padj, file = paste0(path, "/data_subset_pCutOff",pCutOff,".RData"))

  
  ##### Create glm with leave-one-chromosome-out CV and selected features#####
  print("Train and test model...")

  lapply(chroms, function(chrom){
    # Data for training, n-1 chromosomes 
    data_train <- subset(subset_data_features_padj, chr != chrom, select = -c(chr))
    
    # glm
    logR <- glm(formula = mutated ~ ., data = data_train, 
                family = binomial(link = "logit"))
    save(logR, file = paste0(path, "/", chrom,"/logR_featuresp",pCutOff,"_", chrom,".RData"))
  })
  
  
  ##### Test model on remaining chromosome #####
  lapply(chroms, function(chrom){
    print(chrom)
    load(file = paste0(path, "/", chrom,"/logR_featuresp",pCutOff,"_", chrom,".RData"))
    
    # Data for testing
    data_test <- subset(subset_data_features_padj, chr == chrom, select = -c(chr, mutated))
    
    # Prediction
    yhat <- predict(logR, newdata = data_test, type = "response")
    save(yhat, file = paste0(path, "/", chrom,"/yhat_featuresp",pCutOff,"_",chrom, ".RData"))
  })
  
  
  ###### Get all variable coefficients for each leave-one-chromosome-out CV #####
  coefficients_filteredPval <- sapply(chroms, function(chrom){
    
    load(file = paste0(path, "/", chrom,"/logR_featuresp",pCutOff,"_", chrom,".RData"))
    
    return(coef(logR))
  }) %>% as.data.frame()
  
  save(coefficients_filteredPval, file = paste0(path, "/coefficients_allChr_p",pCutOff,".RData"))
  
  
  
  ##### Concatinate yhat #####
  yhat_concat <- sapply(chroms, function(chrom){
    print(chrom)
    load(paste0(path, "/", chrom,"/yhat_featuresp",pCutOff,"_",chrom, ".RData"))
    return(as.vector(yhat))
  }) %>% unlist() %>% as.data.frame() # All predictions combined
  
  save(yhat_concat, file = paste0(path, "/yhat_concat_pCutOff",pCutOff,".RData"))
})




##### PLOT ##### 
lapply(tissues, function(tissue){  
  print(tissue)
  
  ##### Paths & Data #####
  path = paste0("data/rdata/", tissue, "/logRegression")
  
  path_fig = paste0("fig/", tissue, "/logRegression")
  
  
  # Data
  chroms = paste0("chr", seq(1,22)) # List of all chromosomes
  colors_chr <- as.vector(rainbow(22))
  
  load(paste0(path, "/coefficients_allChr_p",pCutOff,".RData"))
  load(paste0(path, "/coefficients_allChr.RData"))
  load(paste0(path, "/pValues_allChr_adjust.RData"))
  
  
  
  # All coefs 
  significant <- c()
  
  for (row in 1:nrow(coefficients)){
    if (rownames(coefficients[row,]) %in% rownames(coefficients_filteredPval)){
      print(paste0(row, " - not filtered"))
      significant <- c(significant, paste0("yes, p-value <= ", pCutOff))
    } else if (!rownames(coefficients[row,]) %in% rownames(coefficients_filteredPval)) {
      print(paste0(row, " - filtered"))
      significant <- c(significant, paste0("no, p-value > ", pCutOff))
    }
  }
  
  coefficients$significant <- significant
  
  coefficients_gg <- cbind(predictors = as.factor(rownames(coefficients)), coefficients)
  rownames(coefficients_gg) <- NULL
  coefficients_gg <- gather(data = coefficients_gg, key = chromosome, value = coefficient, chr1:chr22)
  
  # All p-Values
  pValues_adjust_gg <- cbind(predictors = as.factor(rownames(pValues_adjust)), pValues_adjust)
  rownames(pValues_adjust_gg) <- NULL
  pValues_adjust_gg <- gather(data = pValues_adjust_gg, key = chromosome, value = pValue, chr1:chr9)
  
  
  
  ###### Plots #####
  # All coefficients
  ggplot(coefficients_gg, aes(predictors, coefficient)) + 
    geom_boxplot() + 
    geom_jitter(shape = 20, aes(colour = significant))+ 
    geom_hline(yintercept = 0) + 
    labs(x = "Model feature",y = "Feature coefficient") +
    theme_classic()+
    theme(panel.grid.major = element_line(size = 0.3, linetype = 'solid', colour = "lightgrey"), 
          panel.grid.minor = element_line(size = 0.3, linetype = 'solid', colour = "lightgrey"))+
    theme(text = element_text(size=30))+
    coord_flip()
  
  ggsave(filename = paste0("jitter_pValsFeatures",pCutOff,"_coefs_",tissue,".png"), path = paste0(path_fig),
         width = 55, height = 68, dpi = 400, units = "cm")
  
  
  
  # All p-Values
  ggplot(pValues_adjust_gg, aes(predictors, -log10(pValue))) + 
    geom_hline(yintercept = -log10(0.01), color = "blue") +
    geom_hline(yintercept = -log10(0.05), color = "orange") +
    geom_boxplot() + 
    geom_jitter(shape = 20) + 
    labs(x = "Model feature",y = "-log10(p-value)") +
    theme_light() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +    
    theme(axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme(text = element_text(size=25))+ 
    theme(legend.position="bottom")+
    facet_zoom(ylim = c(0, 4.5))
  
  ggsave(filename = paste0("jitter_pValsFeatures_-log10_zoom_",tissue,".png"), path = paste0(path_fig),
         width = 65, height = 40, dpi = 300, units = "cm")
  
  
  
  # P-Value plot seperated for more detail
  ggplot(pValues_adjust_gg, aes(predictors, -log10(pValue))) + 
    geom_hline(yintercept = -log10(0.01), color = "blue") +
    geom_hline(yintercept = -log10(0.05), color = "orange") +
    geom_boxplot() + 
    geom_jitter(shape = 20) + 
    labs(x = "Model feature",y = "-log10(p-value)") +
    theme_classic() +
    ylim(0,75) +
    theme(panel.grid.major = element_line(size = 0.3, linetype = 'solid', colour = "lightgrey"), 
          panel.grid.minor = element_line(size = 0.3, linetype = 'solid', colour = "lightgrey"))+
    theme(text = element_text(size=45))+
    coord_flip()
  
  ggsave(filename = paste0("jitter_pValsFeatures_-log10_",tissue,"1.png"), path = paste0(path_fig),
         width = 55, height = 70, dpi = 400, units = "cm")
  
  
  ggplot(pValues_adjust_gg, aes(predictors, -log10(pValue))) + 
    geom_boxplot() + 
    geom_jitter(shape = 20) + 
    theme_classic() +
    theme(panel.grid.major = element_line(size = 0.3, linetype = 'solid', colour = "lightgrey"), 
          panel.grid.minor = element_line(size = 0.3, linetype = 'solid', colour = "lightgrey"))+
    ylim(225,300) +
    theme(text = element_text(size=45))+
    coord_flip()
  
 ggsave(filename = paste0("jitter_pValsFeatures_-log10_",tissue,"2.png"), path = paste0(path_fig),
         width = 55, height = 70, dpi = 400, units = "cm")
})

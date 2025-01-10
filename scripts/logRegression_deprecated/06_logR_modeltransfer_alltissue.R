setwd("/cellnet/MutationModel/")

### Packages ###
library(dplyr)
library(ModelMetrics)
library(ROCR)
library(ggplot2)
library(gtools)
library(tibble)
library(magrittr)
###############

##### Select tissue #####
tissues = c("luad",  "skin", "colon", "ovary", "kidney", "prostate", "breast")
tissue_model = tissues[1]
#tissue_apply = tissues[2]
#tissue_comp = c(tissue_apply, paste0(tissue_model,"ON", tissue_apply))
#########################

lapply(tissues[-grep(tissue_model, tissues)], function(tissue_apply){
  print(tissue_apply)
  tissue_comp = c(tissue_model, tissue_apply)
  
  ##### Load data and create paths for processed data and figures #####
  path_model = paste0("data/rdata/", tissue_model, "/logRegression")
  path_apply = paste0("data/rdata/", tissue_apply, "/logRegression")
  
  dir.create(paste0("fig/", tissue_model, "/logRegression/ModelON", tissue_apply), showWarnings=F)
  path_fig = paste0("fig/", tissue_model, "/logRegression/ModelON", tissue_apply)
  
  pCutOff <- 0.01
  
  load(paste0(path_model, "/data_subset_pCutOff",pCutOff,".RData")) # Load subset of data containing only the significant features
  model_features <- colnames(subset_data_features_padj)       # Get all significant features the model was trained on
  
  load(paste0("data/rdata/", tissue_apply, "/completeData_withBinwise.RData")) # Load data to use the model on
  data$context = NULL
  data$pentamer = NULL
  data$trimer = NULL
  data$septamer = NULL
  data$inexon = NULL
  
  data_wChr <- cbind(removed$chr, data) # Data with chromosomes
  colnames(data_wChr)[1] <- "chr"
  
  chroms <- unique(removed$chr)# List of all chromosomes
  #################################################################################

  print("Creating data...")
  # Create data with only signiifcant features to apply the model on 
  data_applymodel <- sapply(model_features, function(feature){ # Select all features used in the model from the data set 
    if(feature %in% colnames(data_wChr)){
      return(data_wChr[feature])
    } else {
      return(rep(x = 0, nrow(data_wChr))) # If feature is not in the data, create a dummy variable 
    }
  }) %>% as.data.frame()
  
  colnames(data_applymodel) <- model_features

  
  
  print("Model testing...")
  ##### MODEL TESTING #####
  
  lapply(chroms, function(chrom){
    load(file = paste0(path_model, "/", chrom,"/logR_featuresp",pCutOff,"_", chrom,".RData")) # Load glm with the current chrom not trained 
    
    # Data for testing; current chromosome data without the columnes chr and mutated 
    data_test <- subset(data_applymodel, chr == chrom, select = -c(chr, mutated))
    
    # Prediction on the current chromosome with which the model was not trained on
    yhat <- predict(logR, newdata = data_test, type = "response") # type=respones to get the actual predicted values between 0 and 1
    
    save(yhat, file = paste0(path_model, "/", chrom,"/yhat_modelON",tissue_apply,"_",chrom,".RData"))
  })
  
  
  print("Model performance...")
  ##### MODEL PERFORMANCES #####
  
  # Calculate LogLoss error on testing chromosome prediction VS true labels
  logloss_CV <- sapply(chroms, function(chrom){
    load(file = paste0(path_model, "/", chrom,"/yhat_modelON",tissue_apply,"_",chrom,".RData"))
    yhat <- as.vector(yhat)
    data_test <- subset(data_applymodel, chr == chrom) # Subset data chromosome-wise
    
    return(logLoss(actual = as.numeric(as.character(data_test$mutated)), predicted = yhat, distribution = "binomial"))
  })
  save(logloss_CV, file = paste0(path_model, "/logLoss_CV_modelON",tissue_apply,".RData"))
  
  
  
  # ROC performances per leave-one-chromosome-out CV
  ROC_CV <- lapply(chroms, function(chrom){ 
    load(file = paste0(path_model, "/", chrom,"/yhat_modelON",tissue_apply,"_",chrom,".RData"))
    
    data_test <- subset(data_applymodel,  chr == chrom)
    
    predict <- prediction(yhat, as.numeric(as.character(data_test$mutated)))
    
    return(performance(predict, "tpr", "fpr"))
  })
  names(ROC_CV) <- chroms
  save(ROC_CV, file = paste0(path_model, "/ROC_CV_modelON",tissue_apply,".RData"))
  
  ## Precision/Recall (PR) for leave-one-chromosome-out CV
  PR_CV <- lapply(chroms, function(chrom){ 
    load(file = paste0(path_model, "/", chrom,"/yhat_modelON",tissue_apply,"_",chrom,".RData"))
    
    data_test <- subset(data_applymodel,  chr == chrom)
    
    predict <- prediction(yhat, as.numeric(as.character(data_test$mutated)))
    
    return(performance(predict, "prec", "rec"))
  })
  names(PR_CV) <- chroms
  save(PR_CV, file = paste0(path_model, "/PR_CV_modelON",tissue_apply,".RData"))
  
  
  # Plot all ROC and PR curves
  load(paste0(path_model, "/ROC_CV_modelON",tissue_apply,".RData"))
  load(paste0(path_model, "/PR_CV_modelON",tissue_apply,".RData"))
  
  colors_chr <- as.vector(rainbow(22))
  
  png(filename = paste0(path_fig,"/ROC_PR_CV_modelON",tissue_apply,".png"), width = 720*2, height = 720, pointsize = 22)
  par(mfrow = c(1,2), cex.lab = 1.4)
  
  # ROC
  plot(ROC_CV$chr1, col = colors_chr[1])
  sapply(seq(2,length(ROC_CV)), function(i){
    plot(ROC_CV[[i]], add = TRUE, col = colors_chr[i])
  })
  abline(coef = c(0,1), lty = 2)
  legend("bottomright", legend = chroms,bg="transparent", col = colors_chr, ncol = 2, lwd = 5, lty=1, box.lty=0, seg.len = 0.8, 
         y.intersp = 1, x.intersp = 0.8, text.width = c(0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1))
  
  plot(PR_CV$chr1, col = colors_chr[1], xlim = c(0,1), ylim = c(0,1))
  sapply(seq(2,length(PR_CV)), function(i){
    plot(PR_CV[[i]], add = TRUE, col = colors_chr[i], xlim = c(0,1), ylim = c(0,1))
  })
  dev.off()
  
  print("Concatenated results...")
  ##### CONCATENATED RESULTS #####
  yhat_concat <- sapply(chroms, function(chrom){
    load(paste0(path_model, "/", chrom,"/yhat_modelON",tissue_apply,"_",chrom,".RData"))
    return(as.vector(yhat))
  }) %>% unlist() %>% as.data.frame() # All predictions combined
  
  save(yhat_concat, file = paste0(path_model, "/yhat_concat_modelON",tissue_apply,".RData"))
  
  
  
  # Evaluation of the concatenated results 
  predict_concat <- prediction(yhat_concat, as.numeric(as.character(data_applymodel$mutated))) # Prediction data for ROC, PR and AUC
  ## ROC 
  ROC_concat <- performance(predict_concat, "tpr", "fpr")
  save(ROC_concat, file = paste0(path_model, "/ROC_concat_modelON",tissue_apply,".RData"))
  ## PR 
  PR_concat <- performance(predict_concat, "prec", "rec")
  save(PR_concat, file = paste0(path_model, "/PR_concat_modelON",tissue_apply,".RData"))
  ## AUROC 
  AUC_concat <- performance(predict_concat, measure = "auc")
  save(AUC_concat, file = paste0(path_model, "/AUC_concat_modelON",tissue_apply,".RData"))
  

  
  # Compare ROC
  allROC_concat <- sapply(tissue_comp, function(tissue){
    if (tissue == tissue_apply){
      load(paste0(path_apply, "/ROC_concat_pCutOff",pCutOff,".RData"))     
      return(ROC_concat)
    } else {
      load(paste0(path_model, "/ROC_concat_modelON",tissue_apply,".RData"))
      return(ROC_concat)
    }
  })
  
  
  
  # Compare PR
  allPR_concat <- sapply(tissue_comp, function(tissue){
    if (tissue == tissue_apply){
      load(paste0(path_apply, "/PR_concat_pCutOff",pCutOff,".RData"))     
      return(PR_concat)
    } else {
      load(paste0(path_model, "/PR_concat_modelON",tissue_apply,".RData"))
      return(PR_concat)
    }
  })
  
  print("Compare all...")
  # Plot ROC and PR curve for the same-tissue model and the other-tissue model performance 
  png(filename = paste0(path_fig, "/ROC_PR_concat_",tissue_comp[1],"VS",tissue_comp[2],".png"), width = 720*2, height = 720, pointsize = 22)
  par(mfrow = c(1,2), cex.lab = 1.4)
  
  # ROC
  plot(allROC_concat[[1]], col = "blue", lwd = 4)
  plot(allROC_concat[[2]], add = TRUE, col = "orange", lwd = 4, lty = 2)
  abline(coef = c(0,1), lty = 2)
  legend(0.65, 0.2, legend = tissue_comp,bg="transparent", col = c("blue", "orange"), lwd = 5, lty=c(1,2), box.lty=0, seg.len = 0.8, 
         y.intersp = 1, x.intersp = 0.8, text.width = c(0.1,0.1))
  
  # PR
  plot(allPR_concat[[1]], col = "blue", lwd = 4, ylim = c(0,1))
  plot(allPR_concat[[2]], add = TRUE, col = "orange", lwd = 4, ylim = c(0,1), lty = 2)
  dev.off()

  
  
  ### Compare LogLoss ###
  all_LogLoss_cv_error <- sapply(tissue_comp, function(tissue){
    if (tissue == tissue_apply){
      load(paste0(path_apply, "/logloss_CV_pCutOff",pCutOff,".RData"))
      return(logloss_CV)
    } else {
      load(file = paste0(path_model, "/logLoss_CV_modelON",tissue_apply,".RData"))
      return(logloss_CV)
    }
  })
  
  
  png(paste0(path_fig, "/boxplot_logLoss_",tissue_comp[1],"VS",tissue_comp[2],".png"), width = 2200, height = 1200, pointsize = 30)
  boxplot(all_LogLoss_cv_error, ylab = "LogLoss", xlab = "model")
  lapply(c(1:length(tissue_comp)), function(i){
    text(y = boxplot.stats(all_LogLoss_cv_error[,i])$stats[3]+0.001, 
         labels = round(boxplot.stats(all_LogLoss_cv_error[,i])$stats[3], digits=4), x = i)
  })
  dev.off()
  
})



###### Compares the results of the one model on all different tissues #####
lapply(tissues, function(tissue_model){
  print(tissue_model)
  path_model = paste0("data/rdata/", tissue_model, "/logRegression")
  
  
  # Compare ROC
  allROC_concat_allTissues <- sapply(tissues, function(tissue){
    if (tissue != tissue_model){
      load(file = paste0(path_model, "/ROC_concat_modelON",tissue,".RData"))
      return(ROC_concat) # Return different model on data performance
    }
  })
  allROC_concat_allTissues <- allROC_concat_allTissues[!sapply(allROC_concat_allTissues,is.null)]
  
  allROC_concat_allTissues_sameModel <- sapply(tissues, function(tissue){
    path = paste0("data/rdata/", tissue, "/logRegression")
    load(file = paste0(path, "/ROC_concat_pCutOff0.01.RData"))
    return(ROC_concat) # Return ROC for model and data from the same tissue
  })
  
  
  # Compare PR
  allPR_concat_allTissues <- sapply(tissues, function(tissue){
    if (tissue != tissue_model){
      load(file = paste0(path_model, "/PR_concat_modelON",tissue,".RData"))
      return(PR_concat) # Return different model on data performance
    }
  })
  allPR_concat_allTissues <- allPR_concat_allTissues[!sapply(allPR_concat_allTissues,is.null)]
  
  allPR_concat_allTissues_sameModel <- sapply(tissues, function(tissue){
    path = paste0("data/rdata/", tissue, "/logRegression")
    load(file = paste0(path, "/PR_concat_pCutOff0.01.RData"))
    return(PR_concat) # Return PR for model and data from the same tissue
  })
  
  colors_t <- as.vector(rainbow(length(tissues)))
  names(colors_t) <- tissues
  
  png(filename = paste0("/cellnet/MutationModel/fig/ROC_PR_allTissues_", tissue_model,"ONall_comparison.png"), width = 720*2, height = 720, pointsize = 22)
  par(mfrow = c(1,2), cex.lab = 1.3)
  
  #ROC
  plot(allROC_concat_allTissues[[1]], col = colors_t[names(allROC_concat_allTissues[1])], lwd = 2)
  sapply(seq(2,length(allROC_concat_allTissues)), function(i){
    plot(allROC_concat_allTissues[[i]], add = TRUE, col = colors_t[names(allROC_concat_allTissues[i])], lwd = 2)
  })
  sapply(setdiff(tissues, tissue_model), function(i){
    plot(allROC_concat_allTissues_sameModel[[i]], add = TRUE, col = colors_t[names(allROC_concat_allTissues[i])], lwd = 2, lty = 3)
  })
  plot(allROC_concat_allTissues_sameModel[[tissue_model]], add = TRUE, lwd = 2, col = colors_t[tissue_model])
  abline(coef = c(0,1), lty = 2)
  legend(0.7, 0.5, legend = tissues, bg="transparent", col = colors_t, lwd = 5, lty=1, box.lty=0, seg.len = 0.8, 
         y.intersp = 1, x.intersp = 0.8, text.width = c(0.1,0.1,0.1,0.1,0.1,0.1, 0.1))
  legend(0.15, 0.18, legend = c(paste0(tissue_model, " model"), "tissue-specific model"), bg="transparent", col = "darkgrey", lwd = 2, lty=c(1,3), box.lty=0, seg.len = 0.8, 
         y.intersp = 1, x.intersp = 0.8, text.width = c(0.1,0.1))
  
  
  #PR
  plot(allPR_concat_allTissues[[1]], col = colors_t[names(allROC_concat_allTissues[1])], lwd = 2, ylim = c(0,1))
  sapply(seq(2,length(allPR_concat_allTissues)), function(i){
    plot(allPR_concat_allTissues[[i]], add = TRUE, col = colors_t[names(allROC_concat_allTissues[i])], lwd = 2, ylim = c(0,1))
  })
  sapply(setdiff(tissues, tissue_model), function(i){
    plot(allPR_concat_allTissues_sameModel[[i]], add = TRUE, col = colors_t[names(allROC_concat_allTissues[i])], lwd = 2, lty = 3, ylim = c(0,1))
  })
  plot(allPR_concat_allTissues_sameModel[[tissue_model]], add = TRUE, lwd = 2, col = colors_t[tissue_model], ylim = c(0,1))
  
  dev.off()
  
})




###### COMAPRE ALL MODEL-TISSUE-APPLICATIONS AT ONCE - HEATMAP #####
fileorder = c()
for(tissue in tissues){
  print(tissue)
  tissue_model = tissue
  
  fileorder = c(fileorder, tissue_model)
  
  for(tissue_apply in tissues[!tissues %in% tissue_model]){
    fileorder = c(fileorder, paste0(tissue_model, "ON", tissue_apply))
  }
}


AUROC = c()
for(tissue in tissues){
  tissue_model = tissue
  path_model = paste0("data/rdata/", tissue_model, "/logRegression")
    
  load(paste0(path_model, "/AUC_concat_pCutOff0.01.RData"))
  AUROC = c(AUROC, AUC_concat@y.values[[1]])
    
  for(tissue_apply in tissues[!tissues %in% tissue_model]){
    load(paste0(path_model, "/AUC_concat_modelON",tissue_apply,".RData"))
    AUROC = c(AUROC, AUC_concat@y.values[[1]])
  }
} 
allAUC_concat <- cbind(fileorder, AUROC) %>% as.data.frame()
allAUC_concat$AUROC <- as.numeric(as.character(allAUC_concat$AUROC))
allAUC_concat$Model = NA
allAUC_concat$Tissue = NA
  
for(i in 1:nrow(allAUC_concat)){ # Scan names to determine the model used on the tissue for each AUROC; ModelONTissue or tissue/data_glm for same model and tissue
  tempName = as.character(allAUC_concat$fileorder[i])
  
  for(tissue in tissues){
    nCharTissue = nchar(tissue)
    
    if(tempName == tissue){ # Same model and tissue 
      allAUC_concat[[3]][i] = tissue
      allAUC_concat[[4]][i] = tissue
      
    } else if(substr(tempName, start = 1, stop = nCharTissue) == tissue){ # Looking for the models used  
      allAUC_concat[[3]][i] = tissue
      
    } else if(substr(tempName, start = nchar(tempName) - nCharTissue+1, stop = nchar(tempName)) == tissue){ # Looking for the tissues the models were used on
      allAUC_concat[[4]][i] = tissue
    }
  }
}  
  
  
  
ggplot(allAUC_concat, aes(Model, Tissue, fill = AUROC)) +
  geom_tile() +
  theme_minimal()+ 
  scale_fill_viridis_c(limits =c(0.4, 0.65))

ggsave(filename = "heatmap_AUROC_allModelsONallTissues_adjScale.png", path = "/cellnet/MutationModel/fig", 
       width = 12, height = 12, units = "cm")



##### BASE-WISE EVALUATION #####

for(tissue_model in tissues){
  print(tissue_model)
  
  lapply(tissues[-grep(tissue_model, tissues)], function(tissue_apply){
  print(tissue_apply)

  ##### Load data and create paths for processed data and figures #####
  path_model = paste0("data/rdata/", tissue_model, "/logRegression")

  bases <- c("A", "C", "T", "G")
  colors_bases <- as.vector(rainbow(4))
  
  load(paste0("data/rdata/", tissue_apply, "/completeData_withBinwise.RData")) # Load data to use the model on
  
  chroms <- unique(removed$chr)# List of all chromosomes
  
  load(paste0(path_model, "/yhat_concat_modelON",tissue_apply,".RData"))
  #################################################################################
  
  # AUROC
  AUROC_bases <- lapply(bases, function(base){ 
    
    data_yhat <- as.data.frame(cbind(yhat_concat, data))
    data_yhat_subset <- subset(data_yhat, ref == base)
    
    predict <- prediction(data_yhat_subset$., as.numeric(as.character(data_yhat_subset$mutated)))
    
    return(performance(predict, measure = "auc"))
  })
  names(AUROC_bases) <- bases
  save(AUROC_bases, file = paste0(path_model, "/AUROC_concat_modelON",tissue_apply,"_basewise.RData"))
  })  

}

##### Plot heatmap - basewise evaluation of all tissue-model-combinations #####
bases <- c("A", "C", "T", "G")

fileorder = c()
for(tissue in tissues){
  print(tissue)
  tissue_model = tissue

  fileorder = c(fileorder, tissue_model)
  
  for(tissue_apply in tissues[!tissues %in% tissue_model]){
    fileorder = c(fileorder, paste0(tissue_model, "ON", tissue_apply))
  }
}


lapply(bases, function(base){
  print(base)
  
  AUROC = c()
  for(tissue in tissues){
    tissue_model = tissue
    path_model = paste0("data/rdata/", tissue_model, "/logRegression")
    
    load(paste0(path_model, "/AUROC_concat_basewise.RData"))
    AUROC = c(AUROC, AUROC_bases[[base]]@y.values[[1]])
      
    for(tissue_apply in tissues[!tissues %in% tissue_model]){
      load(paste0(path_model, "/AUROC_concat_modelON",tissue_apply,"_basewise.RData"))
      AUROC = c(AUROC, AUROC_bases[[base]]@y.values[[1]])
    }
  } 
  allAUC_concat <- cbind(fileorder, AUROC) %>% as.data.frame()
  allAUC_concat$AUROC <- as.numeric(as.character(allAUC_concat$AUROC))
  allAUC_concat$Model = NA
  allAUC_concat$Tissue = NA
  allAUC_concat$count = NA
  
  
  for(i in 1:nrow(allAUC_concat)){ # Scan names to determine the model used on the tissue for each AUROC; ModelONTissue or tissue/data_glm for same model and tissue
    tempName = as.character(allAUC_concat$fileorder[i])
    
    for(tissue in tissues){
      nCharTissue = nchar(tissue)
      
      if(tempName == tissue){ # Same model and tissue 
        allAUC_concat[[3]][i] = tissue
        allAUC_concat[[4]][i] = tissue
        
      } else if(substr(tempName, start = 1, stop = nCharTissue) == tissue){ # Looking for the models used  
        allAUC_concat[[3]][i] = tissue
        
      } else if(substr(tempName, start = nchar(tempName) - nCharTissue+1, stop = nchar(tempName)) == tissue){ # Looking for the tissues the models were used on
        allAUC_concat[[4]][i] = tissue
      }
    }
  }

  #countB <- sapply(tissues, function(tissue){
  #  load(paste0("data/rdata/", tissue, "/completeData_withBinwise.RData")) # Load data to use the model on
  #  return(table(data$ref)[base])
  #})
  #names(countB) <- tissues
#  
#    
#  for(i in 1:nrow(allAUC_concat)){
#    tissuedata = as.character(allAUC_concat$Tissue[i])
#    pos <- match(tissuedata, names(countB))
#    allAUC_concat[[5]][i] = countB[pos]
#  }
  
  
  ggplot(allAUC_concat, aes(Model, Tissue, fill = AUROC)) +
    geom_tile() +
    theme_minimal()+ 
    scale_fill_viridis_c(limits =c(0.45, 0.652)) +
    labs(tag = base)
  
  ggsave(filename = paste0("heatmap_AUROC_allModelsONallTissues_basewise", base, ".png"), path = "/cellnet/MutationModel/fig", 
         width = 12, height = 12, units = "cm")
  
})


countB <- lapply(tissues, function(tissue){
  load(paste0("data/rdata/", tissue, "/completeData_withBinwise.RData")) # Load data to use the model on
  return(table(data$ref))
})
names(countB) <- tissues



setwd("/cellnet/MutationModel/")

### Packages ###
library(tidyr)
library(ROCR)
library(dplyr)
library(tibble)
library(ggplot2)
###############

tissues = c("luad",  "skin", "colon", "ovary", "kidney", "prostate", "breast")

set.seed(22)

##### EXTRACT MINIMUM NUMBER OF POSITIONS AMONG ALL TISSUES #####
countData <- sapply(tissues, function(tissue){
  print(tissue)
  
  path = paste0("data/rdata/", tissue, "/logRegression")
  
  load(paste0(path, "/data_subset_pCutOff0.01.RData"))
  
  num <- c(nrow(subset_data_features_padj), table(subset_data_features_padj$mutated)[1], table(subset_data_features_padj$mutated)[2])
  names(num) <- c("total", "non-mutated", "mutated")
  return(num)
}) %>% as.data.frame()  

downsampling <- countData[which.min(countData["total",])]

##### DOWNSAPMLING TISSUE DATA #####
lapply(tissues, function(tissue){
  print(tissue)
  
  path = paste0("data/rdata/", tissue, "/logRegression")
  
  load(paste0(path, "/data_subset_pCutOff0.01.RData")) # Load full data set 
  
  chroms <- unique(subset_data_features_padj$chr)
  
  
  ### Data downsampling ###
  print("Creating data...")
  
  data_sample0 <- subset(subset_data_features_padj, mutated == "0") %>% rownames_to_column("rowname") %>% 
    sample_n(size = downsampling["non-mutated",]) %>% column_to_rownames("rowname") # Sample non-mutated (0) data

  data_sample1 <- subset(subset_data_features_padj, mutated == "1") %>% rownames_to_column("rowname") %>% 
    sample_n(size = downsampling["mutated",]) %>% column_to_rownames("rowname") # Sample mutated (1) data
  
  data_downsampling <- rbind(data_sample0, data_sample1)
  
  
  save(data_downsampling, file = paste0(path, "/data_downsampling_rep1.RData"))
  


  ### Create glm with leave-one-chromosome-out CV ###
  print("Model training ...")
  
  # Train model on downsampled data
  lapply(chroms, function(chrom){
    
    # Data for training, n-1 chromosomes 
    data_train <- subset(data_downsampling, chr != chrom, select = -c(chr))
    
    # glm
    logR <- glm(formula = mutated ~ ., data = data_train, family = binomial(link = "logit"))
    save(logR, file = paste0(path, "/", chrom,"/logR_downsampling_rep1_", chrom,".RData"))
  })

  

  ### Test model on remaining chromosome ###
  print("Model testing ...")
  
  lapply(chroms, function(chrom){
    load(paste0(path, "/", chrom,"/logR_downsampling_rep1_", chrom,".RData")) # load model trained on downsampled data
    
    # Model trained on downsampled data tested on downsampled data 
    data_test_dd <- subset(data_downsampling, chr == chrom, select = -c(chr, mutated))
    yhat <- predict(logR, newdata = data_test_dd, type = "response")
    save(yhat, file = paste0(path, "/", chrom,"/yhat_downsampling_rep1_",chrom,".RData"))

    # Model trained on downsampled data tested on full data 
    data_test_df <- subset(subset_data_features_padj, chr == chrom, select = -c(chr, mutated))
    yhat<- predict(logR, newdata = data_test_df, type = "response")
    save(yhat, file = paste0(path, "/", chrom,"/yhat_modelDown_testFull_",chrom,".RData"))
    
    # Test full model on downsampled data 
    load(paste0(path, "/", chrom,"/logR_featuresp0.01_", chrom,".RData")) # load model trained on full data
    
    data_test_fd <- subset(data_downsampling, chr == chrom, select = -c(chr, mutated))
    yhat <- predict(logR, newdata = data_test_fd, type = "response")
    save(yhat, file = paste0(path, "/", chrom,"/yhat_modelFull_testDown_",chrom,".RData"))
  })

  

  ### MODEL PERFORMANCE ###
  print("Model performance...")
  
  ## ROC performances per leave-one-chromosome-out CV 
  # Model: Downsampled, Test: Downsampled data
  ROC_CV <- lapply(chroms, function(chrom){ 
    
    load(paste0(path, "/", chrom,"/yhat_downsampling_rep1_",chrom,".RData"))
    
    data_test <- subset(data_downsampling,  chr == chrom)
    
    predict <- prediction(yhat, as.numeric(as.character(data_test$mutated)))
    
    return(performance(predict, "tpr", "fpr"))
  })
  names(ROC_CV) <- chroms
  save(ROC_CV, file = paste0(path, "/ROC_downsampling_rep1_CV.RData"))
  
  
  # Model: Downsampled, Test: Full data
  ROC_CV <- lapply(chroms, function(chrom){ 
    
    load(paste0(path, "/", chrom,"/yhat_modelDown_testFull_",chrom,".RData"))
    
    data_test <- subset(subset_data_features_padj,  chr == chrom)
    
    predict <- prediction(yhat, as.numeric(as.character(data_test$mutated)))
    
    return(performance(predict, "tpr", "fpr"))
  })
  names(ROC_CV) <- chroms
  save(ROC_CV, file = paste0(path, "/ROC_downsampling_modelDown_testFull_CV.RData"))
  
  
  # Model: Full, Test: Downsampled data
  ROC_CV <- lapply(chroms, function(chrom){ 
    
    load(paste0(path, "/", chrom,"/yhat_modelFull_testDown_",chrom,".RData"))
    
    data_test <- subset(data_downsampling,  chr == chrom)
    
    predict <- prediction(yhat, as.numeric(as.character(data_test$mutated)))
    
    return(performance(predict, "tpr", "fpr"))
  })
  names(ROC_CV) <- chroms
  save(ROC_CV, file = paste0(path, "/ROC_downsampling_modelFull_testDown_CV.RData"))
  
  
  ## Precision/Recall performances per leave-one-chromosome-out CV
  # Model: Downsampled, Test: Downsampled data
  PR_CV <- lapply(chroms, function(chrom){ 
    
    load(paste0(path, "/", chrom,"/yhat_downsampling_rep1_",chrom,".RData"))
    
    data_test <- subset(data_downsampling,  chr == chrom)
    
    predict <- prediction(yhat, as.numeric(as.character(data_test$mutated)))
    
    return(performance(predict, "prec", "rec"))
  })
  names(PR_CV) <- chroms
  save(PR_CV, file = paste0(path, "/PR_downsampling_rep1_CV.RData"))
  
  
  # Model: Downsampled, Test: Full data
  PR_CV <- lapply(chroms, function(chrom){ 
    
    load(paste0(path, "/", chrom,"/yhat_modelDown_testFull_",chrom,".RData"))
    
    data_test <- subset(subset_data_features_padj,  chr == chrom)
    
    predict <- prediction(yhat, as.numeric(as.character(data_test$mutated)))
    
    return(performance(predict, "prec", "rec"))
  })
  names(ROC_CV) <- chroms
  save(ROC_CV, file = paste0(path, "/ROC_downsampling_modelDown_testFull_CV.RData"))
  
  
  # Model: Full, Test: Downsampled data
  PR_CV <- lapply(chroms, function(chrom){ 
    
    load(paste0(path, "/", chrom,"/yhat_modelFull_testDown_",chrom,".RData"))
    
    data_test <- subset(data_downsampling,  chr == chrom)
    
    predict <- prediction(yhat, as.numeric(as.character(data_test$mutated)))
    
    return(performance(predict, "prec", "rec"))
  })
  names(ROC_CV) <- chroms
  save(ROC_CV, file = paste0(path, "/ROC_downsampling_modelFull_testDown_CV.RData"))


  ### CONCATENATED RESULTS ###
  print("Concatenate results...")
  
  # Model: Downsampled, Test: Downsampled data
    yhat_concat <- sapply(chroms, function(chrom){
    load(paste0(path, "/", chrom,"/yhat_downsampling_rep1_",chrom,".RData"))
    return(as.vector(yhat))
  }) %>% unlist() %>% as.data.frame() # All predictions combined
  
  save(yhat_concat, file = paste0(path, "/yhat_downsampling_rep1_concat.RData"))
  
  
  # Model: Downsampled, Test: Full data
  yhat_concat <- sapply(chroms, function(chrom){
    load(paste0(path, "/", chrom,"/yhat_modelDown_testFull_",chrom,".RData"))
    return(as.vector(yhat))
  }) %>% unlist() %>% as.data.frame() # All predictions combined
  
  save(yhat_concat, file = paste0(path, "/yhat_downsampling_modelDown_testFull_concat.RData"))
  
  
  # Model: Full, Test: Downsampled data
  yhat_concat <- sapply(chroms, function(chrom){
    load(paste0(path, "/", chrom,"/yhat_modelFull_testDown_",chrom,".RData"))
    return(as.vector(yhat))
  }) %>% unlist() %>% as.data.frame() # All predictions combined
  
  save(yhat_concat, file = paste0(path, "/yhat_downsampling_modelFull_testDown_concat.RData"))
  


  ### CONCATENATED EVALUATION ###

  # Model: Downsampled Test: Downsampeld data #
  load(file = paste0(path, "/yhat_downsampling_rep1_concat.RData"))
  predict_concat_dd <- prediction(yhat_concat, as.numeric(as.character(data_downsampling$mutated))) # Prediction data for ROC, PR and AUC
  
  # ROC
  ROC_concat <- performance(predict_concat_dd, "tpr", "fpr")
  save(ROC_concat, file = paste0(path, "/ROC_downsampling_rep1_concat.RData"))
  # AUC
  AUC_concat <- performance(predict_concat_dd, measure = "auc")
  save(AUC_concat, file = paste0(path, "/AUC_downsampling_rep1_concat.RData"))
  # PR
  PR_concat <- performance(predict_concat_dd, "prec", "rec")
  save(PR_concat, file = paste0(path, "/PR_downsampling_rep1_concat.RData"))
  
  
  # Model: Downsampled Test: Full data #
  load(file = paste0(path, "/yhat_downsampling_modelDown_testFull_concat.RData"))
  predict_concat_df <- prediction(yhat_concat, as.numeric(as.character(subset_data_features_padj$mutated)))
  
  # ROC
  ROC_concat <- performance(predict_concat_df, "tpr", "fpr")
  save(ROC_concat, file = paste0(path, "/ROC_downsampling_modelDown_testFull_concat.RData"))
  # AUC
  AUC_concat <- performance(predict_concat_df, measure = "auc")
  save(AUC_concat, file = paste0(path, "/AUC_downsampling_modelDown_testFull_concat.RData"))
  # PR
  PR_concat <- performance(predict_concat_df, "prec", "rec")
  save(PR_concat, file = paste0(path, "/PR_downsampling_modelDown_testFull_concat.RData"))
  
  
  # Model: Full Test: Downsampeld data #
  load(file = paste0(path, "/yhat_downsampling_modelFull_testDown_concat.RData"))
  predict_concat_fd <- prediction(yhat_concat, as.numeric(as.character(data_downsampling$mutated))) # Prediction data for ROC, PR and AUC
  
  # ROC
  ROC_concat <- performance(predict_concat_fd, "tpr", "fpr")
  save(ROC_concat, file = paste0(path, "/ROC_downsampling_modelFull_testDown_concat.RData"))
  # AUC
  AUC_concat <- performance(predict_concat_fd, measure = "auc")
  save(AUC_concat, file = paste0(path, "/AUC_downsampling_modelFull_testDown_concat.RData"))
  
  # PR
  PR_concat <- performance(predict_concat_fd, "prec", "rec")
  save(PR_concat, file = paste0(path, "/PR_downsampling_modelFull_testDown_concat.RData"))
  
})



##### CORSS-TISSUE APPLICATION #####
for(tissue_model in tissues){
  print(tissue_model)
  
  lapply(tissues[-grep(tissue_model, tissues)], function(tissue_apply){
    print(tissue_apply)

    ##### Load data and  paths for processed data and figures #####
    path_model = paste0("data/rdata/", tissue_model, "/logRegression")
    path_apply = paste0("data/rdata/", tissue_apply, "/logRegression")
    
    # Get significant features used for the model
    load(file = paste0(path_model, "/data_downsampling_rep1.RData")) 
    model_features <- colnames(data_downsampling)
    
    # Get the sample IDs of the previous downsampled data (which only contained significant features of the same-tissue model)
    load(file = paste0(path_apply, "/data_downsampling_rep1.RData")) 
    samples <- rownames(data_downsampling)
    
    # Get data set with all features followed by downscaling so that the exact same samples are selected as before -> downsampled data set with all features instead of only the sig. ones
    load(paste0("data/rdata/", tissue_apply, "/completeData_withBinwise.RData")) # Load data to use the model on
    data <- cbind(removed$chr, data) # Data with chromosomes
    colnames(data)[1] <- "chr"
    
    data_downsampling_apply <- data[match(samples, rownames(data)), ] # Downsample data to the same IDs as used to learn the model
    
    chroms <- unique(data_downsampling_apply$chr) # List of all chromosomes
    #################################################################################
  
    print("Creating data...")
    # Create data with only signifcant features to apply the model on 
    # Downsampled data with matching features
    data_applymodel_down <- sapply(model_features, function(feature){ # Select all features used in the model
      if(feature %in% colnames(data)){
        return(data_downsampling_apply[feature])
      } else {
        return(rep(x = 0, nrow(data_downsampling_apply))) # If feature is not in the data, create a dummy variable with only zeros
      }
    }) %>% as.data.frame()
    
    colnames(data_applymodel_down) <- model_features
    
    # Full data with matching features 
    data_applymodel_full <- sapply(model_features, function(feature){ # Select all features used in the model 
      if(feature %in% colnames(data)){
        return(data[feature])
      } else {
        return(rep(x = 0, nrow(data))) # If feature is not in the data, create a dummy variable with only zeros
      }
    }) %>% as.data.frame()
    
    colnames(data_applymodel_full) <- model_features
    
    
    print("Model testing...")
    ##### MODEL TESTING #####
    
    lapply(chroms, function(chrom){
      load(file = paste0(path_model, "/", chrom,"/logR_downsampling_rep1_", chrom,".RData")) # Load model trained on downsampled data  
      
      # Model trained on downsampled data tested on downsampled data
      data_test_dd <- subset(data_applymodel_down, chr == chrom, select = -c(chr, mutated))
      yhat <- predict(logR, newdata = data_test_dd, type = "response")
      save(yhat, file = paste0(path_model, "/", chrom,"/yhat_downsampling_rep1_modelON",tissue_apply,"_",chrom,".RData"))
      
      # Model trained on downsampled data tested on full data 
      data_test_df <- subset(data_applymodel_full, chr == chrom, select = -c(chr, mutated))
      yhat<- predict(logR, newdata = data_test_df, type = "response")
      save(yhat, file = paste0(path_model, "/", chrom,"/yhat_modelDown_testFull_modelON",tissue_apply,"_",chrom,".RData"))
      
      
      # Test full model on downsampled data 
      load(paste0(path_model, "/", chrom,"/logR_featuresp0.01_", chrom,".RData")) # Load model trained on full data
      
      data_test_fd <- subset(data_applymodel_down, chr == chrom, select = -c(chr, mutated))
      yhat <- predict(logR, newdata = data_test_fd, type = "response")      
      save(yhat, file = paste0(path_model, "/", chrom,"/yhat_modelFull_testDown_modelON",tissue_apply,"_",chrom,".RData"))
      
    })
    

    
    print("Concatenating results and calculating Area under the ROC Curve...")
    ##### CONCATENATED RESULTS #####
    
    # Model: Downsampled, Test: Downsampled data
    yhat_concat <- sapply(chroms, function(chrom){
      load(file = paste0(path_model, "/", chrom,"/yhat_downsampling_rep1_modelON",tissue_apply,"_",chrom,".RData"))
      return(as.vector(yhat))
    }) %>% unlist() %>% as.data.frame() # All predictions combined
    
    save(yhat_concat, file = paste0(path_model, "/yhat_downsampling_rep1_concat_modelON",tissue_apply,".RData"))
    
    ## AUROC 
    predict_concat_dd <- prediction(yhat_concat, as.numeric(as.character(data_applymodel_down$mutated))) 
    AUC_concat <- performance(predict_concat_dd, measure = "auc")
    save(AUC_concat, file = paste0(path_model, "/AUC_downsampling_rep1_concat_modelON",tissue_apply,".RData"))
    
    
    
    
    # Model: Downsampled, Test: Full data
    yhat_concat <- sapply(chroms, function(chrom){
      load(paste0(path_model, "/", chrom,"/yhat_modelDown_testFull_modelON",tissue_apply,"_",chrom,".RData"))
      return(as.vector(yhat))
    }) %>% unlist() %>% as.data.frame() # All predictions combined
    
    save(yhat_concat, file = paste0(path_model, "/yhat_downsampling_modelDown_testFull_concat_modelON",tissue_apply,".RData"))
    
    ## AUROC 
    predict_concat_df <- prediction(yhat_concat, as.numeric(as.character(data_applymodel_full$mutated))) 
    AUC_concat <- performance(predict_concat_df, measure = "auc")
    save(AUC_concat, file = paste0(path_model, "/AUC_downsampling_modelDown_testFull_concat_modelON",tissue_apply,".RData"))
    
    
    
    
    # Model: Full, Test: Downsampled data
    yhat_concat <- sapply(chroms, function(chrom){
      load(paste0(path_model, "/", chrom,"/yhat_modelFull_testDown_modelON",tissue_apply,"_",chrom,".RData"))
      return(as.vector(yhat))
    }) %>% unlist() %>% as.data.frame() # All predictions combined
    
    save(yhat_concat, file = paste0(path_model, "/yhat_downsampling_modelFull_testDown_concat_modelON",tissue_apply,".RData"))
    
    ## AUROC 
    predict_concat_fd <- prediction(yhat_concat, as.numeric(as.character(data_applymodel_down$mutated))) 
    AUC_concat <- performance(predict_concat_fd, measure = "auc")
    save(AUC_concat, file = paste0(path_model, "/AUC_downsampling_modelFull_testDown_concat_modelON",tissue_apply,".RData"))
    
    })
}


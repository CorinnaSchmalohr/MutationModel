setwd("/cellnet/MutationModel/")

### Packages ###
library(tidyr)
library(ROCR)
library(ggplot2)
library(ModelMetrics)
library(gridBase)
library(grid)
###############

##### Select tissue #####
tissues = c("luad",  "skin", "colon", "ovary", "kidney", "prostate", "breast")
tissue = tissues[1]
#########################

pCutOff <- 0.01

##### CHROMOSOME-WISE EVALUATION #####
lapply(tissues, function(tissue){
  print(tissue)
  
  ##### Load preprocessed data with bins and create paths for data and figures #####
  path = paste0("data/rdata/", tissue, "/logRegression")
  
  load(paste0(path, "/data_subset_pCutOff",pCutOff,".RData"))
  
  chroms = paste0("chr", seq(1,22)) # List of all chromosomes

  
  ##### MODEL PERFORMANCE #####
  
  ## LogLoss error on testing chromosome prediction VS true labels ##
  logloss_CV <- sapply(chroms, function(chrom){
    
    load(paste0(path, "/", chrom,"/yhat_featuresp",pCutOff,"_",chrom, ".RData"))
    
    yhat <- as.vector(yhat)
    
    data_test <- subset(subset_data_features_padj,  chr == chrom) 
    
    return(logLoss(actual = as.numeric(as.character(data_test$mutated)), predicted = yhat, distribution = "binomial"))
  })
  save(logloss_CV, file = paste0(path, "/logloss_CV_pCutOff",pCutOff,".RData"))
  
  
  
  ## ROC performances per leave-one-chromosome-out CV ##
  ROC_CV <- lapply(chroms, function(chrom){ 
    
    load(paste0(path, "/", chrom,"/yhat_featuresp",pCutOff,"_",chrom, ".RData"))
    
    data_test <- subset(subset_data_features_padj,  chr == chrom)
    
    predict <- prediction(yhat, as.numeric(as.character(data_test$mutated)))
    
    return(performance(predict, "tpr", "fpr"))
  })
  names(ROC_CV) <- chroms
  save(ROC_CV, file = paste0(path, "/ROC_CV_pCutOff",pCutOff,".RData"))
  
  
  
  ## Precision/Recall (PR) for leave-one-chromosome-out CV ##
  PR_CV <- lapply(chroms, function(chrom){ 
    
    load(paste0(path, "/", chrom,"/yhat_featuresp",pCutOff,"_",chrom, ".RData"))
    
    data_test <- subset(subset_data_features_padj,  chr == chrom)
    
    predict <- prediction(yhat, as.numeric(as.character(data_test$mutated)))
    
    return(performance(predict, "prec", "rec"))
  })
  names(PR_CV) <- chroms
  save(PR_CV, file = paste0(path, "/PR_CV_pCutOff",pCutOff,".RData"))
  
  
  
  ## Area under the ROC (AUROC) curve per leave-one-chromosme-out CV ##
  AUROC_CV <- sapply(chroms, function(chrom){
    
    load(paste0(path, "/", chrom,"/yhat_featuresp",pCutOff,"_",chrom, ".RData"))
    
    data_test <- subset(subset_data_features_padj,  chr == chrom)
    
    predict <- prediction(yhat, as.numeric(as.character(data_test$mutated)))
    
    performAUC <- performance(predict, measure = "auc")
    
    return(performAUC@y.values[[1]])
  })
  names(AUROC_CV) <- chroms
  save(AUROC_CV, file = paste0(path, "/AUROC_CV_pCutOff",pCutOff,".RData"))
  
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
  
  pCutOffs <- c("0.01", "0.05", "not filtered")
  
  load(paste0(path, "/ROC_CV_pCutOff",pCutOff,".RData"))       # Load ROC
  load(paste0(path, "/PR_CV_pCutOff",pCutOff,".RData"))        # Load PR

  # LogLoss for each filtering and w/o filtering
  all_LogLoss_CV <- sapply(pCutOffs, function(cutoff){
    print(cutoff)
    if (cutoff == "0.01"){
      load(paste0(path, "/logloss_CV_pCutOff",cutoff,".RData"))
      return(logloss_CV)
    } else if (cutoff == "0.05"){
      paste0(path, "/logloss_CV_pCutOff",cutoff,".RData")
      return(logloss_CV)
    } else if (cutoff == "not filtered"){
      load(paste0(path, "/logLoss_CV.RData"))
      return(logloss_CV)
    }
  })
  
  
  ##### Plots #####
  # ROC and PR 
  png(filename = paste0(path_fig, "/ROC_PR_pCutOff",pCutOff,"_CV_",tissue,".png"), width = 720*2, height = 720, pointsize = 22)
  par(mfrow = c(1,2), cex.lab = 1.4)
  
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
  
  # logLoss for each method
  png(paste0(path_fig, "/boxplot_logLoss_reducedVSfull_",tissue,".png"), width = 720*2, height = 720, pointsize = 22)
  boxplot(all_LogLoss_CV, ylab = "LogLoss error", xlab = "p-Value cut-off")
  dev.off()
})  



##### CONCATENATED EVALUATION #####
lapply(tissues, function(tissue){
  print(tissue)
  
  ##### Load preprocessed data with bins and create paths for data and figures #####
  path = paste0("data/rdata/", tissue, "/logRegression")
  
  load(paste0("data/rdata/", tissue, "/completeData_withBinwise.RData"))
  
  load(file = paste0(path, "/yhat_concat_pCutOff",pCutOff,".RData"))
  
  
  ##### Concatinated model performance #####
  predict_concat <- prediction(yhat_concat, as.numeric(as.character(data$mutated))) # Prediction data for ROC, PR and AUC
  
  # ROC
  ROC_concat <- performance(predict_concat, "tpr", "fpr")
  save(ROC_concat, file = paste0(path, "/ROC_concat_pCutOff",pCutOff,".RData"))
  
  # AUC
  AUC_concat <- performance(predict_concat, measure = "auc")
  save(AUC_concat, file = paste0(path, "/AUC_concat_pCutOff",pCutOff,".RData"))
  
  # PR
  PR_concat <- performance(predict_concat, "prec", "rec")
  save(PR_concat, file = paste0(path, "/PR_concat_pCutOff",pCutOff,".RData"))
  
})



##### PLOTS CONCATENATED #####
lapply(tissues, function(tissue){  
  print(tissue)
  
  ##### Paths & Data #####
  path = paste0("data/rdata/", tissue, "/logRegression")
  
  path_fig = paste0("fig/", tissue, "/logRegression")
  
  
  # Data
  load(paste0("data/rdata/", tissue, "/completeData_withBinwise.RData"))
  chroms = unique(removed$chr)

  colors_chr <- as.vector(rainbow(22))

  load(file = paste0(path, "/ROC_concat_pCutOff",pCutOff,".RData"))
  load(file = paste0(path, "/PR_concat_pCutOff",pCutOff,".RData"))
  
  pCutOffs <- c("0.01", "0.05", "not filtered")

  
  ##### Plots #####

  # Grouped violin plot 
  yhat_compare <- sapply(pCutOffs, function(cutoff){
    print(cutoff)
    if (cutoff == "0.01"){
      load(paste0(path, "/yhat_concat_pCutOff",cutoff,".RData"))
      return(yhat_concat)
      
    }  else if (cutoff == "0.05"){
      load(paste0(path, "/yhat_concat_pCutOff",cutoff,".RData"))
      return(yhat_concat)
      
    } else if (cutoff == "not filtered"){
      load(paste0(path, "/yhat_concat.RData"))
      return(yhat_concat)
      
    }
  }) %>% as.data.frame()
  names(yhat_compare) <- pCutOffs
  yhat_compare <- cbind(position = rownames(data), label = data$mutated, yhat_compare)
  yhat_compare <- gather(data = yhat_compare, key = pcutoff, value = prediction, "0.01":"not filtered")
  
  
  ggplot(yhat_compare, aes(x=label, y=prediction, fill = pcutoff), ylim = c(0,1)) + 
    geom_violin(position = "dodge") + 
    labs(x = "True label", y = "Predicted outcome") +
    scale_y_continuous(limits = c(0, 1)) + 
    scale_fill_viridis_d()+ 
    theme_classic() +
    theme(text = element_text(size=15), panel.grid.major = element_line(size = 0.3, linetype = 'solid', colour = "lightgrey"), 
          panel.grid.minor = element_line(size = 0.3, linetype = 'solid', colour = "lightgrey"))+ 
    theme(legend.position="bottom")
  
  
  ggsave(filename = paste0("violin_predVStrue_reducedVSfull_concat_",tissue,".png"), path = path_fig,
         width = 30, height = 20, dpi = 300, units = "cm")  
  
  
  # ROC and PR plot
  png(paste0(path_fig,"/ROC_PR_concat_pCutOff",pCutOff,"_",tissue,".png"), width = 720*2, height = 720, pointsize = 22)
  par(mfrow = c(1,2), cex.lab = 1.4)
  
  # ROC
  plot(ROC_concat, lwd = 2.5)
  abline(coef = c(0,1), lty = 2)
  
  # PR
  plot(PR_concat, xlim = c(0,1), ylim = c(0,1), lwd = 2.5)
  dev.off()
  
  
  ##### All three in one plot #####
  png(filename = paste0(path_fig, "/ROC_PR_violin_concat_pCutOff",pCutOff,"_",tissue,".png"), width = 720*2, height = 720*2, pointsize = 22)
  par(mfrow = c(2,2), cex.lab = 1.4, cex = 1.01)
  
  # ROC
  plot(ROC_concat, lwd = 2.5)
  abline(coef = c(0,1), lty = 2)
  
  # PR
  plot(PR_concat, xlim = c(0,1), ylim = c(0,1), lwd = 2.5)
  
  # Mutation ratio
  plot.new()              
  vps <- baseViewports()
  pushViewport(vps$figure)
  vp1 <-plotViewport(c(1.8,1,0,1))
  
  p <- ggplot(yhat_compare, aes(x=label, y=prediction, fill = pcutoff), ylim = c(0,1)) + 
    geom_violin(position = "dodge") + 
    labs(x = "True label", y = "Predicted outcome") +
    scale_y_continuous(limits = c(0, 1)) + 
    scale_fill_viridis_d()+ 
    theme_classic() +
    theme(text = element_text(size=30), panel.grid.major = element_line(size = 0.3, linetype = 'solid', colour = "lightgrey"), 
          panel.grid.minor = element_line(size = 0.3, linetype = 'solid', colour = "lightgrey"))+ 
    theme(legend.position="bottom")
  
  print(p, vp = vp1)
  dev.off()
})  
setwd("/cellnet/MutationModel/")

### Packages ###
library(ROCR)
library(ggplot2)
library(gridBase)
library(grid)
###############

##### Select tissue #####
tissues = c("luad",  "skin", "colon", "ovary", "kidney", "prostate", "breast")
tissue = tissues[1]
#########################



##### CHROMOSOME-WISE EVALUATION #####
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
  
  data_wChr <- cbind(removed$chr, data) # Dataset with corresponding chromosome assignmet
  colnames(data_wChr)[1] <- "chr"
  
  chroms = paste0("chr", seq(1,22)) # List of all chromosomes
  

  
  ## Mutation ratio per chromosome ##
  mutation_ratio <- sapply(chroms, function(chrom){
    subset_chrom <- subset(data_wChr,  chr == chrom) # Subset data chromosome-wise
    mut <- sum(as.numeric(as.character(subset_chrom$mutated))) # Number of all mutations
    non_mut <- length(as.numeric(as.character(subset_chrom$mutated)))- mut # Number of all not-mutated loci
    ratio <- mut/non_mut
    return(c(mut, non_mut, ratio)) 
  })
  rownames(mutation_ratio) <- c("mutated", "non-mutated", "ratio")
  save(mutation_ratio, file = paste0(path, "/mutation_ratio.RData"))
  

  
  ##### MODEL PERFORMANCE #####
  
  ## LogLoss error on testing chromosome prediction VS true labels ##
  logloss_CV <- sapply(chroms, function(chrom){
    
    load(paste0(path, "/", chrom,"/yhat_",chrom,".RData"))
    
    yhat <- as.vector(yhat)
    
    data_test <- subset(data_wChr,  chr == chrom) 
    
    return(logLoss(actual = as.numeric(as.character(data_test$mutated)), predicted = yhat, distribution = "binomial"))
  })
  save(logloss_CV, file = paste0(path, "/logloss_CV.RData"))
  
  
  
  ## ROC performances per leave-one-chromosome-out CV ##
  ROC_CV <- lapply(chroms, function(chrom){ 

    load(paste0(path, "/", chrom,"/yhat_",chrom,".RData"))
    
    data_test <- subset(data_wChr,  chr == chrom)
    
    predict <- prediction(yhat, as.numeric(as.character(data_test$mutated)))
    
    return(performance(predict, "tpr", "fpr"))
  })
  names(ROC_CV) <- chroms
  save(ROC_CV, file = paste0(path, "/ROC_CV.RData"))
  
  
  
  ## Precision/Recall (PR) for leave-one-chromosome-out CV ##
  PR_CV <- lapply(chroms, function(chrom){ 
    
    load(paste0(path, "/", chrom,"/yhat_",chrom, ".RData"))
    
    data_test <- subset(data_wChr,  chr == chrom)
    
    predict <- prediction(yhat, as.numeric(as.character(data_test$mutated)))
    
    return(performance(predict, "prec", "rec"))
  })
  names(PR_CV) <- chroms
  save(PR_CV, file = paste0(path, "/PR_CV.RData"))



  ## Area under the ROC (AUROC) curve per leave-one-chromosme-out CV ##
  AUROC_CV <- sapply(chroms, function(chrom){
    
    load(paste0(path, "/", chrom,"/yhat_",chrom,".RData"))
    
    data_test <- subset(data_wChr,  chr == chrom)
    
    predict <- prediction(yhat, as.numeric(as.character(data_test$mutated)))
    
    performAUC <- performance(predict, measure = "auc")
    
    return(performAUC@y.values[[1]])
  })
  names(AUROC_CV) <- chroms
  save(AUROC_CV, file = paste0(path, "/AUROC_CV.RData"))
  
})  


lm_eqn <- function(m){
  eq <- substitute(~~italic(r)^2~"="~r2, 
                   list(r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}


##### PLOT ##### 
lapply(tissues, function(tissue){  
  print(tissue)
  
  ##### Paths & Data #####
  path = paste0("data/rdata/", tissue, "/logRegression")
  
  dir.create(paste0("fig/", tissue, "/logRegression"), showWarnings=F)
  path_fig = paste0("fig/", tissue, "/logRegression")
  
  
  # Data
  chroms = paste0("chr", seq(1,22)) # List of all chromosomes
  colors_chr <- as.vector(rainbow(22))

  load(paste0(path, "/ROC_CV.RData"))       # Load ROC
  load(paste0(path, "/PR_CV.RData"))        # Load PR
  load(paste0(path, "/AUROC_CV.RData"))     # Load AUROC
  load(paste0(path, "/mutation_ratio.RData"))   # Mutation ratio

  
  ##### Plots #####
  # ROC and PR 
  png(filename = paste0(path_fig, "/ROC_PR_CV_",tissue,".png"), width = 720*2, height = 720, pointsize = 22)
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
  
  
  # Mutation ratio VS AUROC
  data_plot_AUC_ratio <- data.frame(mutation_ratio[3,], colSums(mutation_ratio[1:2,]), AUROC_CV)
  names(data_plot_AUC_ratio) <- c("ratio", "size", "AUC")
  
  ggplot(data = data_plot_AUC_ratio, aes(ratio, AUC)) + 
    labs(x = "TP/TN mutation ratio", y = "AUROC") + 
    geom_smooth(method = "lm", alpha = 0.3, se = FALSE) + 
    geom_point() +
    theme_classic() +
    theme(text = element_text(size=15), panel.grid.major = element_line(size = 0.3, linetype = 'solid', colour = "lightgrey"), 
          panel.grid.minor = element_line(size = 0.3, linetype = 'solid', colour = "lightgrey"))
  
  ggsave(filename = paste0("mutRatio_per_AUROC_CV_",tissue,".png"), path = path_fig,
         width = 19, height = 19, dpi = 400, units = "cm")
  

  
  ##### All three in one plot #####
  png(filename = paste0(path_fig, "/ROC_PR_size_mutRatio_CV_",tissue,".png"), width = 720*2, height = 720*2, pointsize = 22)
  layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE))
  par(cex.lab = 1.4, cex = 1.01)
  
  
  # ROC
  plot(ROC_CV$chr1, col = colors_chr[1])
  col_num <- 2
  for (i in ROC_CV[2:length(ROC_CV)]){
    plot(i, add = TRUE, col = colors_chr[col_num])
    col_num <- col_num + 1
  }
  abline(coef = c(0,1), lty = 2)
  legend("bottomright", legend = chroms,bg="transparent", col = colors_chr, ncol = 2, lwd = 5, lty=1, box.lty=0, seg.len = 0.8, 
         y.intersp = 1, x.intersp = 0.8, text.width = c(0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1))
  
  # PR
  plot(PR_CV$chr1, col = colors_chr[1], xlim = c(0,1), ylim = c(0,1))
  col_num <- 2
  for (i in PR_CV[2:length(PR_CV)]){
    plot(i, add = TRUE, col = colors_chr[col_num], xlim = c(0,1), ylim = c(0,1))
    col_num <- col_num + 1
  }
  
  
  plot.new()              
  vps <- baseViewports()
  pushViewport(vps$figure) 
  vp1 <-plotViewport(c(1.8,1,0,1))
  # Chromosome size
  size_lm=lm(size~AUC,data=data_plot_AUC_ratio)

  pp <- ggplot(data = data_plot_AUC_ratio, aes(size, AUC)) + 
    labs(x = "Data size", y = "AUROC") + 
    geom_smooth(method = "lm", alpha = 0.3, se = FALSE) + 
    geom_point(size = 3) +
    theme_classic() +
    theme(text = element_text(size=30))+ 
    geom_text(x=8500, y=0.68, label = lm_eqn(size_lm), parse = TRUE, size = 9)
  
  
  # Mutation ratio
  tp_lm=lm(ratio~AUC,data=data_plot_AUC_ratio)

  p <- ggplot(data = data_plot_AUC_ratio, aes(ratio, AUC)) + 
    labs(x = "TP/TN mutation ratio", y = "AUROC") + 
    geom_smooth(method = "lm", alpha = 0.3, se = FALSE) + 
    geom_point(size = 3) +
    theme_classic() +
    theme(text = element_text(size=30))+ 
    geom_text(x=0.75, y=0.68, label = lm_eqn(tp_lm), parse = TRUE, size = 9)

  
  print(pp+p, vp = vp1)
  dev.off()
  
})  



##### CONCATENATED EVALUATION #####
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
  
  data_wChr <- cbind(removed$chr, data) # Dataset with corresponding chromosome assignmet
  colnames(data_wChr)[1] <- "chr"
  
  chroms = unique(removed$chr)
  
  
  
  # Concatinate yhat
  yhat_concat <- sapply(chroms, function(chrom){
    print(chrom)
    load(paste0(path, "/", chrom,"/yhat_",chrom,".RData"))
    return(as.vector(yhat))
  }) %>% unlist() %>% as.data.frame() # All predictions combined
  
  save(yhat_concat, file = paste0(path, "/yhat_concat.RData"))
  
  
  
  ##### Concatinated model performance #####
  predict_concat <- prediction(yhat_concat, as.numeric(as.character(data$mutated))) # Prediction data for ROC, PR and AUC
  
  # ROC
  ROC_concat <- performance(predict_concat, "tpr", "fpr")
  save(ROC_concat, file = paste0(path, "/ROC_concat.RData"))
  
  # AUC
  AUC_concat <- performance(predict_concat, measure = "auc")
  save(AUC_concat, file = paste0(path, "/AUC_concat.RData"))
  
  # PR
  PR_concat <- performance(predict_concat, "prec", "rec")
  save(PR_concat, file = paste0(path, "/PR_concat.RData"))

})



##### PLOTS CONCATENATED #####
lapply(tissues, function(tissue){  
  print(tissue)
  
  ##### Paths & Data #####
  path = paste0("data/rdata/", tissue, "/logRegression")
  
  dir.create(paste0("fig/", tissue, "/logRegression"), showWarnings=F)
  path_fig = paste0("fig/", tissue, "/logRegression")
  
  
  # Data
  load(paste0("data/rdata/", tissue, "/completeData_withBinwise.RData"))
  
  chroms = paste0("chr", seq(1,22)) # List of all chromosomes
  colors_chr <- as.vector(rainbow(22))
  colors_t <- as.vector(rainbow(7))
  names(colors_t) <- tissues
  
  load(paste0(path, "/yhat_concat.RData"))
  load(file = paste0(path, "/ROC_concat.RData"))
  load(file = paste0(path, "/PR_concat.RData"))
  
  
  
  ##### Plots #####
  # Violin Plot: predicted VS true values
  
  data_all_pred_label <- cbind(yhat_concat, data$mutated) # Prediction VS true label 
  colnames(data_all_pred_label) <- c("prediction", "label")
  
  ggplot(data_all_pred_label, aes(label, prediction), ylim = c(0,1)) + 
    geom_violin(fill=colors_t[tissue], alpha = 0.6) + 
    labs(x = "True label", y = "Predicted outcome") +
    scale_y_continuous(limits = c(0, 1)) +
    theme_classic() +
    theme(text = element_text(size=15), panel.grid.major = element_line(size = 0.3, linetype = 'solid', colour = "lightgrey"), 
          panel.grid.minor = element_line(size = 0.3, linetype = 'solid', colour = "lightgrey"))
  
  ggsave(filename = paste0("violin_yhatVStrue_concat_",tissue,".png"), path = paste0(path_fig),
         width = 19, height = 19, dpi = 400, units = "cm")
  
  
  
  # ROC and PR plot
  png(paste0(path_fig,"/ROC_PR_concat_",tissue,".png"), width = 720*2, height = 720, pointsize = 22)
  par(mfrow = c(1,2), cex.lab = 1.4)
  
  # ROC
  plot(ROC_concat, lwd = 2.5)
  abline(coef = c(0,1), lty = 2)
  
  # PR
  plot(PR_concat, xlim = c(0,1), ylim = c(0,1), lwd = 2.5)
  dev.off()
  
  
  ##### All three in one plot #####
  png(filename = paste0(path_fig, "/ROC_PR_violin_concat_",tissue,".png"), width = 720*2, height = 720*2, pointsize = 22)
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
  
  p <- ggplot(data_all_pred_label, aes(label, prediction), ylim = c(0,1)) + 
    geom_violin(fill=colors_t[tissue], alpha = 0.6) + 
    labs(x = "True label", y = "Predicted outcome") +
    scale_y_continuous(limits = c(0, 1)) +
    theme_classic() +
    theme(text = element_text(size=30), panel.grid.major = element_line(size = 0.3, linetype = 'solid', colour = "lightgrey"), 
          panel.grid.minor = element_line(size = 0.3, linetype = 'solid', colour = "lightgrey"))
  
  print(p, vp = vp1)
  dev.off()
})  
  





##### BASE-WISE EVALUATION #####
lapply(tissues, function(tissue){
  print(tissue)
  
  ##### Load preprocessed data with bins and create paths for data and figures #####
  path = paste0("data/rdata/", tissue, "/logRegression")
  
  dir.create(paste0("fig/", tissue, "/logRegression"), showWarnings=F)
  path_fig = paste0("fig/", tissue, "/logRegression")
  
  load(paste0("data/rdata/", tissue, "/completeData_withBinwise.RData"))
  
  bases <- c("A", "C", "T", "G")
  colors_bases <- as.vector(rainbow(4))
  
  load(paste0(path, "/yhat_concat.RData"))
  
  
  ##### Concatinated base-wise model performance #####
  
  # ROC
  ROC_bases <- lapply(bases, function(base){ 
    
    data_yhat <- as.data.frame(cbind(yhat_concat, data))
    data_yhat_subset <- subset(data_yhat, ref == base)
    
    predict <- prediction(data_yhat_subset$., as.numeric(as.character(data_yhat_subset$mutated)))
    
    return(performance(predict, "tpr", "fpr"))
  })
  names(ROC_bases) <- bases
  
  # PR
  PR_bases <- lapply(bases, function(base){ 
    
    data_yhat <- as.data.frame(cbind(yhat_concat, data))
    data_yhat_subset <- subset(data_yhat, ref == base)
    
    predict <- prediction(data_yhat_subset$., as.numeric(as.character(data_yhat_subset$mutated)))
    
    return(performance(predict, "prec", "rec"))
  })
  names(PR_bases) <- bases
  
  # AUROC
  AUROC_bases <- lapply(bases, function(base){ 
    
    data_yhat <- as.data.frame(cbind(yhat_concat, data))
    data_yhat_subset <- subset(data_yhat, ref == base)
    
    predict <- prediction(data_yhat_subset$., as.numeric(as.character(data_yhat_subset$mutated)))
    
    return(performance(predict, "auc"))
  })
  names(AUROC_bases) <- bases
  save(AUROC_bases, file = paste0(path, "/AUROC_concat_basewise.RData"))
  
  
  
  ##### Plot #####
  png(filename = paste0(path_fig, "/ROC_PR_basewise_",tissue,".png"), width = 720*2, height = 720, pointsize = 22)
  par(mfrow = c(1,2), cex.lab = 1.4)
  
  plot(ROC_bases$A, col = colors_bases[1], lwd = 2.5)
  sapply(seq(2,length(ROC_bases)), function(i){
    plot(ROC_bases[[i]], add = TRUE, col = colors_bases[i], lwd = 2.5)
  })
  abline(coef = c(0,1), lty = 2)
  legend("bottomright", legend = bases, col = colors_bases, lty=1, cex=1.3, lwd = 3,bg="transparent")

  plot(PR_bases$A, col = colors_bases[1], xlim = c(0,1), ylim = c(0,1), lwd = 2.5)
  sapply(seq(2,length(PR_bases)), function(i){
    plot(PR_bases[[i]], add = TRUE, col = colors_bases[i], xlim = c(0,1), ylim = c(0,1), lwd = 2.5)
  })
  
  dev.off()

})

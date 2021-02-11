setwd("/cellnet/MutationModel/")
### Packages ###
library(ggplot2)
library(tidyr)
###############

##### Select tissue #####
tissues = c("luad",  "skin", "colon", "ovary", "kidney", "prostate", "breast")
tissue = tissues[1]
########################################

##### Data with bins #####
path_lasso = paste0("data/nina/rdata_withBinwise/", tissue)

path = paste0("data/rdata/", tissue, "/logRegression")

path_fig = paste0("fig/", tissue, "/logRegression")

load(paste0("data/rdata/", tissue, "/completeData_withBinwise.RData"))

data$context = NULL
data$pentamer = NULL
data$trimer = NULL
data$septamer = NULL
data$inexon = NULL

###########################

chroms = unique(removed$chr)

mod_methods = c("lasso", "glm")

##### LogLoss for all regularization methods ##### 
all_LogLoss_cv_error <- sapply(mod_methods, function(method){
  print(method)
  if (method == "lasso"){
    load(paste0(path_lasso, "/logLoss_",method,".RData"))
    return(loco_logLoss_cv_error)
  } else {
    paste0(path, "/logloss_CV_pCutOff0.01.RData")    
    return(logloss_CV)
  }
})


# Plot all logLoss for each method
png(paste0(path_fig, "/boxplot_logLoss_logRVSlasso.png"), width = 720*2, height = 720, pointsize = 22)
boxplot(all_LogLoss_cv_error, ylab = "LogLoss error")

lapply(c(1:2), function(i){
  text(y = boxplot.stats(all_LogLoss_cv_error[,i])$stats[c(1,3,5)]-0.0007, 
       labels = round(boxplot.stats(all_LogLoss_cv_error[,i])$stats[c(1,3,5)], digits=4), x = i)
})

lapply(c(1:2), function(i){
  text(y = boxplot.stats(all_LogLoss_cv_error[,i])$out, 
       labels = names(boxplot.stats(all_LogLoss_cv_error[,i])$out), x = i+0.08)
})

dev.off()



# ROC for glm and lasso
col2 = c("orangered", "blue")

png(filename = paste0(path_fig, "/ROC_concat_glmVSlasso.png"), width = 720, height = 720, pointsize = 22)

load(file = paste0(path_lasso, "/performROC_concat_lasso.RData"))
plot(performanceROC_concat, col = col2[1], alpha = 0.2, lwd = 5)
load(paste0(path, "/ROC_CV_pCutOff0.01.RData"))
plot(ROC_concat, add = TRUE, col = col2[2], alpha = 0.2, lwd = 5, lty = 2)

abline(coef = c(0,1))
legend("bottomright", legend = mod_methods, col = col2, cex=1.3, lwd = 5)

dev.off()



all_AUC_concat <- sapply(mod_methods, function(method){
  print(method)
  if (method == "lasso"){
    load(paste0(path_lasso, "/performAUC_concat_",method,".RData"))
    return(as.numeric(performanceAUC_concat@y.values))
    
  } else {
    load(paste0(path, "/AUROC_CV_pCutOff0.01.RData"))
    return(as.numeric(AUC_concat@y.values))
  }
})



# AUC for lasso and glm
all_AUC_concat <- as.data.frame(all_AUC_concat)
colnames(all_AUC_concat) <- "AUC" 
all_AUC_concat$regMethod <- rownames(all_AUC_concat)
rownames(all_AUC_concat) <- NULL

ggplot(as.data.frame(all_AUC_concat), aes(regMethod, AUC)) + 
  geom_point() + 
  geom_text(label = round(all_AUC_concat$AUC, digits = 4), 
            nudge_y = -0.0025) + 
  ylim(c(0.6,0.68))

ggsave("AUC_CV_logRVSlasso.png", path = paste0(path_fig),
       width = 20, height = 20, dpi = 300, units = "cm")


##### Grouped violin plot #####
yhat_compare <- sapply(mod_methods, function(method){
  print(method)
  if (method == "lasso"){
    load(paste0(path_lasso, "/yhat_concat_lasso.RData"))
    return(yhat_concat)
    
  } else {
    load(paste0(path, "/yhat_concat_pCutOff0.01.RData"))
    return(yhat_concat)
  }
})
yhat_compare <- cbind(position = rownames(data), label = data$mutated, as.data.frame(yhat_compare))
yhat_compare <- gather(data = yhat_compare, key = type, value = prediction, lasso:glm)


ggplot(yhat_compare, aes(x=label, y=prediction, fill = type), ylim = c(0,1)) + 
  geom_violin(position = "dodge", alpha=0.5) + 
  labs(x = "True value", y = "Predicted value") +
  scale_y_continuous(limits = c(0, 1)) + 
  scale_fill_viridis_d()+ 
  theme_light()

ggsave(filename = paste0("violin_predVStrue_logRVSlasso.png"), path = paste0(path_fig),
       width = 30, height = 20, dpi = 300, units = "cm")




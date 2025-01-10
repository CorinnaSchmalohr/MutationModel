source("scripts/05_analysis/00_NamesAndColors.R")
tissues = c("brain","breast", "colon","kidney","liver", "lung","ovary", #"esophagus",
            "prostate", "skin")
t2T = t2T[tissues]
library(ROCR)
library(ggplot2)
cr= "chr1"

# only train data subsampled #####
subTrain = sapply(tissues, function(tissue2predict){
  cat("data: ", tissue2predict, "\n")
  load(paste0("data/procData/traindata/traindata_processed_",
              tissue2predict, ".RData"))
  tempData = dat[datchroms == cr,]
  truePreds = tempData$mutated
  cat("model: ")
  perfs = sapply(tissues, function(predictingTissue){
    cat(predictingTissue, " ")
    load(paste0("data/rdata/GLMmodel/", predictingTissue, "_subsampledNsamples_", cr, "_sig.RData"))
    load(paste0("data/procData/traindata/traindata_processed_",
                predictingTissue, "_subsampledNsamples.RData"))
    
    # Just transfer on the sig. predictors
    load(file = paste0("data/rdata/GLMmodel/", predictingTissue,"_pvals_subsampledNsamples.RData"))
    predictors = names(rowMeans(pvals)[rowMeans(pvals) < 0.05])
    testData = tempData[,colnames(tempData) %in% predictors]
    toAdd = predictors[!predictors %in% colnames(testData)]
    for(x in toAdd){
      testData[x] = rep(mean(dat[,x]), nrow(testData))
    }
    yHat = predict(logR, newdata = testData, type = "response")
    temp = prediction(yHat, truePreds)
    roc = performance(temp, "tpr", "fpr")
    pr = performance(temp,"prec", "rec")
    auc = performance(temp,"auc")
    perf = list(roc = roc, pr = pr, auc = auc)
    return(perf)
  }, simplify=F)
  cat("\n")
  return(perfs)
}, simplify = F)
save(subTrain, file = "data/rdata/GLMCrossTissue_subTrainNsamples.RData")
#####

# only test data subsampled #####
subTest = sapply(tissues, function(tissue2predict){
  cat("data: ", tissue2predict, "\n")
  load(paste0("data/procData/traindata/traindata_processed_",
              tissue2predict, "_subsampledNsamples.RData"))
  tempData = dat[datchroms == cr,]
  truePreds = tempData$mutated
  cat("model: ")
  perfs = sapply(tissues, function(predictingTissue){
    cat(predictingTissue, " ")
    load(paste0("data/procData/traindata/traindata_processed_",
                predictingTissue, ".RData"))
    load(paste0("data/rdata/GLMmodel/", predictingTissue,
                "_", cr, "_sig.RData"))
    
    # Just transfer on the sig. predictors
    load(file = paste0("data/rdata/GLMmodel/", predictingTissue,"_pvals.RData"))
    predictors = names(rowMeans(pvals)[rowMeans(pvals) < 0.05])
    testData = tempData[,colnames(tempData) %in% predictors]
    toAdd = predictors[!predictors %in% colnames(testData)]
    for(x in toAdd){
      testData[x] = rep(mean(dat[,x]), nrow(testData))
    }
    yHat = predict(logR, newdata = testData, type = "response")
    temp = prediction(yHat, truePreds)
    roc = performance(temp, "tpr", "fpr")
    pr = performance(temp,"prec", "rec")
    auc = performance(temp,"auc")
    perf = list(roc = roc, pr = pr, auc = auc)
    return(perf)
  }, simplify=F)
  cat("\n")
  return(perfs)
}, simplify = F)
save(subTest, file = "data/rdata/GLMCrossTissue_subTestNsamples.RData")
#####

# both subsampled #####
subBoth = sapply(tissues, function(tissue2predict){
  cat("data: ", tissue2predict, "\n")
  load(paste0("data/procData/traindata/traindata_processed_",
              tissue2predict, "_subsampledNsamples.RData"))
  tempData = dat[datchroms == cr,]
  truePreds = tempData$mutated
  cat("model: ")
  perfs = sapply(tissues, function(predictingTissue){
    cat(predictingTissue, " ")
    load(paste0("data/rdata/GLMmodel/", predictingTissue, "_subsampledNsamples_", cr, "_sig.RData"))
    load(paste0("data/procData/traindata/traindata_processed_",
                predictingTissue, "_subsampledNsamples.RData"))
    
    # Just transfer on the sig. predictors
    load(file = paste0("data/rdata/GLMmodel/", predictingTissue,"_pvals_subsampledNsamples.RData"))
    predictors = names(rowMeans(pvals)[rowMeans(pvals) < 0.05])
    testData = tempData[,colnames(tempData) %in% predictors]
    toAdd = predictors[!predictors %in% colnames(testData)]
    for(x in toAdd){
      testData[x] = rep(mean(dat[,x]), nrow(testData))
    }
    yHat = predict(logR, newdata = testData, type = "response")
    temp = prediction(yHat, truePreds)
    roc = performance(temp, "tpr", "fpr")
    pr = performance(temp,"prec", "rec")
    auc = performance(temp,"auc")
    perf = list(roc = roc, pr = pr, auc = auc)
    return(perf)
  }, simplify=F)
  cat("\n")
  return(perfs)
}, simplify = F)
save(subBoth, file = "data/rdata/GLMCrossTissue_subBothNsamples.RData")
#####

# collect ROCs #####
load("data/rdata/GLMCrossTissue_subTrainNsamples.RData")
load("data/rdata/GLMCrossTissue_subTestNsamples.RData")
load("data/rdata/GLMCrossTissue_subBothNsamples.RData")

ROCs_subTrain = data.frame(do.call(rbind,lapply(tissues, function(tissue2predict){
  t(sapply(tissues, function(predictingTissue){
    temp = subTrain[[tissue2predict]][[predictingTissue]]$auc@y.values[[1]]
    return(c("trainedOn" = predictingTissue, 
             "testedOn" = tissue2predict, "AUROC" = temp))
  }, USE.NAMES= F))
})), stringsAsFactors=F)
ROCs_subTrain$AUROC = as.numeric(ROCs_subTrain$AUROC)
ROCs_subTrain$rank = NA

ROCs_subTest = data.frame(do.call(rbind,lapply(tissues, function(tissue2predict){
  t(sapply(tissues, function(predictingTissue){
    temp = subTest[[tissue2predict]][[predictingTissue]]$auc@y.values[[1]]
    return(c("trainedOn" = predictingTissue, 
             "testedOn" = tissue2predict, "AUROC" = temp))
  }, USE.NAMES= F))
})), stringsAsFactors=F)
ROCs_subTest$AUROC = as.numeric(ROCs_subTest$AUROC)
ROCs_subTest$rank = NA

ROCS_subBoth = data.frame(do.call(rbind,lapply(tissues, function(tissue2predict){
  t(sapply(tissues, function(predictingTissue){
    temp = subBoth[[tissue2predict]][[predictingTissue]]$auc@y.values[[1]]
    return(c("trainedOn" = predictingTissue, 
             "testedOn" = tissue2predict, "AUROC" = temp))
  }, USE.NAMES= F))
})), stringsAsFactors=F)
ROCS_subBoth$AUROC = as.numeric(ROCS_subBoth$AUROC)
ROCS_subBoth$rank = NA

# Add rowrank
for(tissue in tissues){
  ROCs_subTrain[ROCs_subTrain$testedOn == tissue,]$rank[order(-ROCs_subTrain[ROCs_subTrain$testedOn == tissue,]$AUROC)] = 1:length(tissues)
  ROCs_subTest[ROCs_subTest$testedOn == tissue,]$rank[order(-ROCs_subTest[ROCs_subTest$testedOn == tissue,]$AUROC)] = 1:length(tissues)
  ROCS_subBoth[ROCS_subBoth$testedOn == tissue,]$rank[order(-ROCS_subBoth[ROCS_subBoth$testedOn == tissue,]$AUROC)] = 1:length(tissues)
}

#####

# visualize #####
###
upperL = 0.69
lowerL = 0.48
###

ggplot(ROCs_subTrain, aes(trainedOn, testedOn, fill= AUROC)) +
  geom_tile() +
  scale_fill_viridis_c(limits=c(lowerL,upperL)) +
  geom_text(aes(label = rank), color = "white", size = 2.5)+
  labs(y = "Mutations and predictors tested on", x = "Model trained on") +  
  scale_x_discrete(labels=t2T, guide = guide_axis(n.dodge = 1, angle = 40)) + 
  scale_y_discrete(labels = t2T)+
  theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0), 
        panel.background = element_rect(fill = NA),
        axis.ticks.length= unit(0, "cm"),
        axis.line= element_line(linewidth = NA))
#ggsave(filename="fig/CrossTissue_subTrainNsamples.png", 
#       width = 12, height = 12, units = "cm")
ggsave(filename="fig/CrossTissue_subTrainNsamples.pdf", 
       width = 13, height = 10, units = "cm", bg = "white", dpi = 200)

ggplot(ROCs_subTest, aes(trainedOn, testedOn, fill= AUROC)) +
  geom_tile() +
  scale_fill_viridis_c(limits=c(lowerL,upperL)) +
  geom_text(aes(label = rank), color = "white", size = 2.5)+
  labs(y = "Mutations and predictors tested on", x = "Model trained on") +  
  scale_x_discrete(labels=t2T, guide = guide_axis(n.dodge = 1, angle = 40)) + 
  scale_y_discrete(labels = t2T)+
  theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0), 
        panel.background = element_rect(fill = NA),
        axis.ticks.length= unit(0, "cm"),
        axis.line= element_line(linewidth = NA))
#ggsave(filename="fig/CrossTissue_subTestNsamples.png", 
#       width = 12, height = 12, units = "cm")
ggsave(filename="fig/CrossTissue_subTestNsamples.pdf", 
       width = 13, height = 10, units = "cm", bg = "white", dpi = 200)

ggplot(ROCS_subBoth, aes(trainedOn, testedOn, fill= AUROC)) +
  geom_tile() +
  scale_fill_viridis_c(limits=c(lowerL,upperL)) +
  geom_text(aes(label = rank), color = "white", size = 2.5)+
  labs(y = "Mutations and predictors tested on", x = "Model trained on") +  
  scale_x_discrete(labels=t2T, guide = guide_axis(n.dodge = 1, angle = 40)) + 
  scale_y_discrete(labels = t2T)+
  theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0), 
        panel.background = element_rect(fill = NA),
        axis.ticks.length= unit(0, "cm"),
        axis.line= element_line(linewidth = NA))
#ggsave(filename="fig/CrossTissue_subBothNsamples.png", 
#       width = 12, height = 12, units = "cm")
ggsave(filename="fig/CrossTissue_subBothNsamples.pdf", 
       width = 13, height = 10, units = "cm", bg = "white", dpi = 200)
#####


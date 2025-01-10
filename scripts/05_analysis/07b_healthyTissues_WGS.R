# preparation #####
.libPaths(new = "/data/public/cschmalo/R-4.1.2/")
library(ROCR)
library(ggplot2)
library(sinaplot)
library(ranger)
library(plotrix)
source("./scripts/05_analysis/00_NamesAndColors.R")
modelTissues = c("brain","breast","esophagus",
                 "kidney", "liver",# "ovary",
                 "prostate", "skin")
testTissues =c("adipocytes", "adrenal_gland", "bladder", 
               "bonemarrow", "brain", "breast", "colon", 
               "esophagus", "fibroblast", "heart",  "kidney",
               "liver", "lung", "pancreas", "placenta", "prostate", "rectum", 
               "skeletal_muscle", "skin", "small_intestine", "spleen", 
               "stomach", "testicle", "thyroid", "tonsil", "ureter") 
nThread = 14
capitalize <- function(x) {
  s <- strsplit(x, "_")[[1]]
  paste(toupper(substring(s, 1, 1)), substring(s, 2),
        sep = "", collapse = " ")
}
lightCol = function(x,alpha){
  x = col2rgb(x)
  rgb(red=x[1,], green=x[2,], blue=x[3,], alpha=alpha, maxColorValue=255)
}
dir.create("fig/healthyTissuesWGS", showWarnings = F)
plotEnding = "_20241209"
#####

# tissue-specific models #####
print("tissue-specific models")
perfsTissueSpecific = sapply(rev(modelTissues), function(tissue){
  print(tissue)
  # load healthy tissue trainingData
  load(paste0("data/MutTables/healthyTissuesWGS/", tissue, "_Muts_mapped_processed.RData")) #dat,datchroms
  # load model
  load(paste0("data/Modeling/WholeGenomeData/RF/", tissue, "_finalModel.RData"))
  # predict
  yhat = predict(rf, data = dat, num.threads = nThread, type = "response")
  rm(rf, datchroms); gc()
  # compute performance
  temp = prediction(pred = yhat$predictions[,2],labels = dat$mutated)
  rm(yhat, dat);
  roc = performance(temp,  "tpr", "fpr")
  auc = performance(temp, "auc")
  pr = performance(temp, "prec", "rec")
  res = list(roc = roc, pr = pr, auc = auc)
  return(res)
}, simplify = F)
save(perfsTissueSpecific, file="data/Modeling/healthyTissues/perfsTissueSpecific_WGS.RData")
#####

# allTissue model #####
print("general model")
load("data/Modeling/WholeGenomeData/RF/TissueCombination_finalModel.RData")

perfsAllTissueModel = sapply(testTissues, function(tissue){
  print(tissue)
  # load healthy tissue trainingData
  load(paste0("data/MutTables/healthyTissuesWGS/", tissue, 
              "_Muts_mapped_processed.RData")) #dat,datchroms
  
  # predict
  yhat = predict(rf, data = dat, num.threads = nThread, type = "response")
  # compute performance
  temp = prediction(pred = yhat$predictions[,2],labels = dat$mutated)
  roc = performance(temp,  "tpr", "fpr")
  auc = performance(temp, "auc")
  pr = performance(temp, "prec", "rec")
  return(list(roc = roc, pr = pr, auc = auc))
}, simplify = F)
save(perfsAllTissueModel, file="data/Modeling/healthyTissues/perfsAllTissueModel_WGS.RData")
#####


# load/compute training performance of general and tissue-specific models #####
print("compute performance")
load("data/Modeling/WholeGenomeData/RF/TissueCombination_predictions.RData") # predictions
GeneralPredConcat = do.call(rbind,predictions)
temp = prediction(GeneralPredConcat$pred, GeneralPredConcat$label)
GeneralSelfPerf = list(roc = performance(temp, "tpr", "fpr"),
                       pr = performance(temp,"prec", "rec"),
                       auc = performance(temp,"auc")@y.values[[1]])

load("data/Modeling/WholeGenomeData/RF/ROC_PR_RF_concat.RData")

#####

# plotting #####
print("plotting")
load("data/Modeling/healthyTissues/perfsTissueSpecific_WGS.RData")
load("data/Modeling/healthyTissues/perfsAllTissueModel_WGS.RData")
# ROC, and PR comparing for each tissue: 
# general model, training perf of general model, 
# tissue-specific model (when available), perf on tissue-specific training
# ROC
png(paste0("fig/healthyTissuesWGS/ROC", plotEnding, ".png"),
    width = 1200, height = 800, pointsize = 15)
par(mfrow = c(4,7), mar = c(1,1,2,1), oma = c(3,3,0,0))
dumpVar = sapply(1:length(testTissues), function(i){
  tissue = testTissues[i]
  plot(NA, xlim = c(0,1), ylim = c(0,1),
       xlab = "", ylab = "", las = 1, xaxt = "n", yaxt = "n", main = capitalize(tissue))
  abline(0,1, col = "dark grey", lty = 2, lwd = 1.5)
  if((i-1)%%7 == 0){
    axis(2, las = 1)
    title(ylab = "TPR", xpd = NA)
  }
  if(i>(3*7)){
    axis(1, las = 1)
    title(xlab = "FPR", xpd = NA)
  }
  # general model
  lines(perfsAllTissueModel[[tissue]]$roc@x.values[[1]],
        perfsAllTissueModel[[tissue]]$roc@y.values[[1]])
  # general model training performance
  lines(GeneralSelfPerf$roc@x.values[[1]],
        GeneralSelfPerf$roc@y.values[[1]], lty = 2)
  # tissue-specific model
  if(tissue %in% modelTissues){
    lines(perfsTissueSpecific[[tissue]]$roc@x.values[[1]],
          perfsTissueSpecific[[tissue]]$roc@y.values[[1]], col = tissueCols[tissue])
    lines(ROC_PR_RF_concat[[tissue]]$roc@x.values[[1]],
          ROC_PR_RF_concat[[tissue]]$roc@y.values[[1]], col = tissueCols[tissue], lty = 2)
  }
})
plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
legend("topleft", lty = c(1,2,1,2), col = c("black", "black", "red", "red"),xpd = NA,
       legend = c("General model", "Training performance general model", 
                  "Tissue-specific model", "Training performance tissue-specific model"))
dev.off()
# PR
png(paste0("fig/healthyTissuesWGS/PR", plotEnding, ".png"),
    width = 1200, height = 800, pointsize = 15)
par(mfrow = c(4,7), mar = c(1,1,2,1), oma = c(3,3,0,0))
dumpVar = sapply(1:length(testTissues), function(i){
  tissue = testTissues[i]
  plot(NA, xlim = c(0,1), ylim = c(0,1),
       xlab = "", ylab = "", las = 1, xaxt = "n", yaxt = "n", main = capitalize(tissue))
  if((i-1)%%7 == 0){
    axis(2, las = 1)
    title(ylab = "Precision", xpd = NA)
  }
  if(i>(3*7)){
    axis(1, las = 1)
    title(xlab = "Recall", xpd = NA)
  }
  # general model
  lines(perfsAllTissueModel[[tissue]]$pr@x.values[[1]],
        perfsAllTissueModel[[tissue]]$pr@y.values[[1]])
  # general model training performance
  lines(GeneralSelfPerf$pr@x.values[[1]],
        GeneralSelfPerf$pr@y.values[[1]], lty = 2)
  # tissue-specific model
  if(tissue %in% modelTissues){
    lines(perfsTissueSpecific[[tissue]]$pr@x.values[[1]],
          perfsTissueSpecific[[tissue]]$pr@y.values[[1]], col = tissueCols[tissue])
    lines(ROC_PR_RF_concat[[tissue]]$pr@x.values[[1]],
          ROC_PR_RF_concat[[tissue]]$pr@y.values[[1]], col = tissueCols[tissue], lty = 2)
  }
})
plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
legend("topleft", lty = c(1,2,1,2), col = c("black", "black", "red", "red"),xpd = NA,
       legend = c("General model", "Training performance general model", 
                  "Tissue-specific model", "Training performance tissue-specific model"))
dev.off()

# AUC general model
png(paste0("fig/healthyTissuesWGS/AUC", plotEnding, ".png"),
    width = 1200, height = 800, pointsize = 15)
AUCs = rbind("General" = sapply(perfsAllTissueModel, function(x)x$auc@y.values[[1]]),
             "Tissue-specific" = sapply(perfsTissueSpecific, function(x)x$auc@y.values[[1]])[testTissues])
par(mar = c(3,7,0.5,0.5))
AUCs = AUCs[,ncol(AUCs):1]
axisbuffer = 0.02
offset = 0.4
barplot(AUCs-offset, 
        beside = T, horiz = T, las = 1, legend.text = T, 
        col = rbind("grey20",  tissueCols[colnames(AUCs)]),
        names.arg = sapply((colnames(AUCs)), capitalize), xaxt = "n")
# axis(1)
temp =axTicks(1)
axis(1,at = temp, labels = c(0,temp[-1]+offset))
axis.break(axis = 1, breakpos = 0.01, style = "gap", bgcol = "grey")
dev.off()

# AUCs vs nMuts
nMutsHealthy = sapply(testTissues, function(tissue){
  load(paste0("data/MutTables/healthyTissuesWGS/", tissue, "_Muts_mapped_processed.RData")) #dat,datchroms
  return(nrow(dat))
})
png(paste0("fig/healthyTissuesWGS/AUCvsNmuts", plotEnding, ".png"),
    width = 600, height = 600, pointsize = 15)
plot(nMutsHealthy, AUCs["General", names(nMutsHealthy)], ylab = "AUC", xlab = "Test data size")
points(nMutsHealthy, AUCs[2, names(nMutsHealthy)], 
       col = tissueCols[names(nMutsHealthy)], pch = 19)
legend("bottomright", c("General model", "Tissue-specific model"), col = c("black", "grey30"), pch = c(1,19))
dev.off()
####

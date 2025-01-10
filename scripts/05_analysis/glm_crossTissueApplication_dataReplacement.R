source("scripts/05_analysis/00_NamesAndColors.R")
library(ROCR)
library(ggplot2)
cr= "chr1"
load("data/rdata/dataInfos.RData")
mindatsize = min(sapply(dataInfos, function(x)x$nMuts))

# variant 1: e.g. test lung model on lung mutations with kidney predictors #####
print("variant 1")
CrossTissuePreds = sapply(tissues, function(muttissue){
  cat("muts and model from ", muttissue, "\n")
  # load glm
  load(paste0("data/rdata/GLMmodel/", muttissue,
              "_", cr, "_sig.RData"))
  # Get sig. predictors
  load(file = paste0("data/rdata/GLMmodel/", muttissue,"_pvals.RData"))
  predictors = names(rowMeans(pvals)[rowMeans(pvals) < 0.05])
  
  cat("predictors from ")
  perfs = sapply(tissues, function(predtissue){
    cat(predtissue, " ")
    if(predtissue == muttissue){
      load(paste0("data/procData/traindata/traindata_processed_",
                  muttissue, ".RData"))
    }else{
      load(paste0("data/procData/traindata_crosstissue/traindata_crosstissue_",
                  muttissue, "With",predtissue,"Preds_processed", ".RData"))
    }
    
    
    testData = dat[datchroms == cr,]
    truePreds = testData$mutated
    # add dummy variables for predictors present only in predtissue
    #colnames(testData)[colnames(testData) == "aPhased_repeats"] = "aPhased_repeates"
    testData = testData[,colnames(testData) %in% predictors]
    toAdd = predictors[!predictors %in% colnames(testData)]
    load(paste0("data/procData/traindata/traindata_processed_",
                muttissue, ".RData"))
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
save(CrossTissuePreds, file = "data/rdata/GLMCrossTissuePreds.RData")
#####


# variant 2: e.g. test lung model on kidney mutations with lung preds #####
print("variant 2")
CrossTissueMuts = sapply(tissues, function(muttissue){
  cat("preds and model from ", muttissue, "\n")
  # load glm
  load(paste0("data/rdata/GLMmodel/", muttissue,
              "_", cr, "_sig.RData"))
  # Get sig. predictors
  load(file = paste0("data/rdata/GLMmodel/", muttissue,"_pvals.RData"))
  predictors = names(rowMeans(pvals)[rowMeans(pvals) < 0.05])
  
  cat("muts from ")
  perfs = sapply(tissues, function(predtissue){
    cat(predtissue, " ")
    if(predtissue == muttissue){
      load(paste0("data/procData/traindata/traindata_processed_",
                  predtissue, ".RData"))
    }else{
      load(paste0("data/procData/traindata_crosstissue/traindata_crosstissue_",
                  predtissue, "With",muttissue,"Preds_processed", ".RData"))
    }
    testData = dat[datchroms == cr,]
    truePreds = testData$mutated
    # add dummy variables for predictors present only in predtissue
    #colnames(testData)[colnames(testData) == "aPhased_repeats"] = "aPhased_repeates"
    testData = testData[,colnames(testData) %in% predictors]
    toAdd = predictors[!predictors %in% colnames(testData)]
    load(paste0("data/procData/traindata/traindata_processed_",
                muttissue, ".RData"))
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
save(CrossTissueMuts, file = "data/rdata/GLMCrossTissueMuts.RData")
#####


# visualize results #####
###
upperL = 0.69
lowerL = 0.48
###

# variant 1
#load("data/rdata/GLMCrossTissuePreds.RData")
RocsPreds = data.frame(do.call(rbind,lapply(tissues, function(muttissue){
  t(sapply(tissues, function(predtissue){
    temp = CrossTissuePreds[[muttissue]][[predtissue]]$auc@y.values[[1]]
    return(c("predictorSource" = predtissue, 
             "MutationAndModelSource" = muttissue, "AUROC" = temp))
  }, USE.NAMES= F))
})), stringsAsFactors=F)
RocsPreds$AUROC = as.numeric(RocsPreds$AUROC)
RocsPreds$rank <- NA
for(tissue in tissues){
  RocsPreds[RocsPreds$predictorSource == tissue,]$rank[order(-RocsPreds[RocsPreds$predictorSource == tissue,]$AUROC)] = 1:length(tissues)
}

ggplot(RocsPreds, aes(x = predictorSource, y = MutationAndModelSource, fill= AUROC)) +
  geom_tile() +
  scale_fill_viridis_c(limits=c(lowerL,upperL)) +
  labs(x = "Predictors tested on", y = "Model & mutations trained and tested on") + 
  scale_x_discrete(labels=t2T, guide = guide_axis(n.dodge = 1, angle = 40)) + 
  scale_y_discrete(labels = t2T)+
  geom_text(aes(label = rank), color = "white", size = 2.5)+
  theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0), 
        panel.background = element_rect(fill = NA),
        axis.ticks.length= unit(0, "cm"),
        axis.line= element_line(linewidth = NA))
#ggsave(filename="fig/CrossTissuePreds_glm_191022.png", 
#       width = 12, height = 12, units = "cm")
ggsave(filename="fig/CrossTissuePreds.pdf", 
       width = 13, height = 10, units = "cm", bg = "white", dpi = 200)

# it looks like we have slightly higher performance than for the
# original cross-tissue application. Is that the case?
#load("data/rdata/GLMperformanceCrossTissue.RData")
#ROCs = data.frame(do.call(rbind,lapply(tissues, function(tissue2predict){
#  t(sapply(tissues, function(predictingTissue){
#    temp = GLMperformanceCrossTissue[[tissue2predict]][[predictingTissue]]$auc@y.values[[1]]
#    return(c("trainedOn" = predictingTissue, 
#             "testedOn" = tissue2predict, "AUROC" = temp))
#  }, USE.NAMES= F))
#})), stringsAsFactors=F)
#ROCs$AUROC = as.numeric(ROCs$AUROC)
#all.equal(ROCs$trainedOn, RocsPreds$predictorSource)
#all.equal(ROCs$testedOn, RocsPreds$MutationAndModelSource)
#
#table(ROCs$AUROC <= RocsPreds$AUROC)
#png("fig/crossTissue_origVsRocsPreds_191022.png")
#plot(ROCs$AUROC, RocsPreds$AUROC, xlim = c(0.47,0.7), ylim = c(0.47,0.7), 
#     xlab = "original AUROCs", ylab = "AUROCs from exchanged predictors",
#     col = (RocsPreds$AUROC > ROCs$AUROC)+3)
#abline(0,1)
#text(x = 0.55, y = 0.68, paste0("n = ",sum(RocsPreds$AUROC > ROCs$AUROC)), col = 4, font = 2)
#text(x = 0.68, y = 0.55, paste0("n = ",sum(RocsPreds$AUROC <= ROCs$AUROC)), col = 3, font = 2)
#
#dev.off()


# variant 2
#load("data/rdata/GLMCrossTissueMuts.RData")
RocsMuts = data.frame(do.call(rbind,lapply(tissues, function(muttissue){
  t(sapply(tissues, function(predtissue){
    temp = CrossTissueMuts[[muttissue]][[predtissue]]$auc@y.values[[1]]
    return(c("MutationSource" = predtissue, 
             "PredictorAndModelSource" = muttissue, "AUROC" = temp))
  }, USE.NAMES= F))
})), stringsAsFactors=F)
RocsMuts$AUROC = as.numeric(RocsMuts$AUROC)
RocsMuts$rank <- NA
for(tissue in tissues){
  RocsMuts[RocsMuts$MutationSource == tissue,]$rank[order(-RocsMuts[RocsMuts$MutationSource == tissue,]$AUROC)] = 1:length(tissues)
}
ggplot(RocsMuts, aes(y = MutationSource, x = PredictorAndModelSource, fill= AUROC)) +
  geom_tile() +
  scale_fill_viridis_c(limits=c(lowerL,upperL)) +
  labs(y = "Mutations tested on", x = "Model & predictors trained and tested on") + 
  scale_x_discrete(labels=t2T, guide = guide_axis(n.dodge = 1, angle = 40)) + 
  scale_y_discrete(labels = t2T)+
  geom_text(aes(label = rank), color = "white", size = 2.5)+
  theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0), 
        panel.background = element_rect(fill = NA),
        axis.ticks.length= unit(0, "cm"),
        axis.line= element_line(linewidth = NA))
#ggsave(filename="fig/CrossTissueMuts_191022.png", 
#       width = 12, height = 12, units = "cm")
ggsave(filename="fig/CrossTissueMuts.pdf", 
       width = 13, height = 10, units = "cm", bg = "white", dpi = 200)
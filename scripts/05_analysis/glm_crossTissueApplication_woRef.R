library(ROCR)
library(ggplot2)
library("gridExtra")   
source("scripts/05_analysis/00_NamesAndColors.R")

# Test remove ref base and look if cross-tissue performance changes 
tissues = c("lung", "kidney")

for(tissue in tissues){
  print(tissue)
  # prepare data for this tissue######
  load(paste0(paste0("data/procData/traindata/traindata_processed_", tissue, ".RData")))
  dat = dat[colnames(dat) != "ref_CG"]
  chroms = unique(datchroms)
  #####
  # glm #####
  print("training models")
  temp = lapply(chroms, function(cr){
    cat(cr, ' ')
    trainData = dat[datchroms != cr,]
    logR = glm(formula = mutated ~ ., data = trainData, 
               family = binomial(link = "logit"))
    save(logR, file = paste0("data/rdata/GLMmodel/", tissue, "_", cr, "_woRef.RData"))
    return(NULL)
  })
  cat('\n')
  #####
  
  # calculate variable  p-values ####
  print("pvals")
  pvals <- sapply(chroms, function(cr){
    cat(cr, ' ')
    load(paste0("data/rdata/GLMmodel/", tissue, "_", cr, "_woRef.RData")) 
    temp = drop1(object = logR, test = "LRT")
    drop1_features <- setNames(temp$`Pr(>Chi)`, rownames(temp))
    return(drop1_features)
  })
  pvals <- pvals[-1,] # remove first empty row
  pvals <- as.data.frame(pvals)
  save(pvals, file = paste0("data/rdata/GLMmodel/", tissue,
                            "_pvals_woRef.RData"))
  #####
  
  # Remove non-significant features and retrain glm #####
  #load(file = paste0("data/rdata/GLMmodel/", tissue,"_pvals.RData"))
  meanPvals = rowMeans(pvals)
  sigFreatures = names(meanPvals[meanPvals < 0.05])
  dat = cbind(dat[sigFreatures], mutated = dat$mutated)
  
  # glm ##
  print("training significant models")
  temp = lapply(chroms, function(cr){
    cat(cr, ' ')
    trainData = dat[datchroms != cr,]
    logR = glm(formula = mutated ~ ., data = trainData, 
               family = binomial(link = "logit"))
    save(logR, file = paste0("data/rdata/GLMmodel/", tissue, "_", cr, "_sig_woRef.RData"))
    return(NULL)
  })
  cat('\n')
  
  # predict on held-out chromosomes #
  print("significant predictions")
  predictions = lapply(chroms, function(cr){
    cat(cr, ' ')
    load(paste0("data/rdata/GLMmodel/", tissue, "_", cr, "_sig_woRef.RData"))
    testData = dat[datchroms == cr,]
    yhat = predict(logR, newdata = testData, type = "response")
    temp = data.frame(pred  = yhat,label = testData$mutated)
    return(temp)
  })
  names(predictions) = chroms
  save(predictions, file = paste0("data/rdata/GLMmodel/", tissue,
                                  "_predictions_sig_woRef.RData"))
  cat('\n')
  #####
  
  print("done with this tissue")
}


# apply models to data from other tissues and get performance #####
cr = "chr1"
GLMperformanceCrossTissue_woRef = sapply(tissues, function(tissue2predict){
  cat("data: ", tissue2predict, "\n")
  load(paste0("data/procData/traindata/traindata_processed_",
              tissue2predict, ".RData"))
  dat = dat[colnames(dat) != "ref_CG"]
  tempData = dat[datchroms == cr,]
  truePreds = tempData$mutated
  cat("model: ")
  perfs = sapply(tissues, function(predictingTissue){
    cat(predictingTissue, " ")
    load(paste0("data/rdata/GLMmodel/", predictingTissue, "_", 
                cr, "_sig_woRef.RData"))
    load(paste0("data/procData/traindata/traindata_processed_",
                predictingTissue, ".RData"))
    dat = dat[colnames(dat) != "ref_CG"]
    # Just transfer on the sig. predictors
    load(file = paste0("data/rdata/GLMmodel/", predictingTissue,"_pvals_woRef.RData"))
    predictors = names(rowMeans(pvals)[rowMeans(pvals) < 0.05])
    testData = tempData[,colnames(tempData) %in% predictors]
    toAdd = predictors[!predictors %in% colnames(testData)]
    for(x in toAdd){ # Add mean value of predictors missing for the tissue to predict on (from predicting tissue)
      testData[x] = rep(mean(dat[,x]), nrow(testData))
    }
    # Predict on new transfered data 
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
PerfsSelf_woRef = sapply(tissues, function(tissue){
  GLMperformanceCrossTissue_woRef[[tissue]][[tissue]]
}, simplify = F)
load(file = "data/rdata/GLMperformanceCrossTissue.RData")
#####

# image of ROCs
ROCs = data.frame(do.call(rbind,lapply(tissues, function(tissue2predict){
  t(sapply(tissues, function(predictingTissue){
    temp = GLMperformanceCrossTissue[[tissue2predict]][[predictingTissue]]$auc@y.values[[1]]
    temp_woRef = GLMperformanceCrossTissue_woRef[[tissue2predict]][[predictingTissue]]$auc@y.values[[1]]
    return(c("trainedOn" = predictingTissue, 
             "testedOn" = tissue2predict, "AUROC" = temp, "AUROC_woRef" = temp_woRef))
  }, USE.NAMES= F))
})), stringsAsFactors=F)
ROCs$AUROC = as.numeric(ROCs$AUROC)
ROCs$AUROC_woRef = as.numeric(ROCs$AUROC_woRef)

ROCs$rank <- NA
for(tissue in tissues){
  ROCs[ROCs$testedOn == tissue,]$rank[order(-ROCs[ROCs$testedOn == tissue,]$AUROC)] = 1:length(tissues)
}
ROCs$rank_woRef <- NA
for(tissue in tissues){
  ROCs[ROCs$testedOn == tissue,]$rank_woRef[order(-ROCs[ROCs$testedOn == tissue,]$AUROC_woRef)] = 1:length(tissues)
}

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

ggp1 = ggplot(ROCs, aes(trainedOn, testedOn, fill= AUROC_woRef)) +
  geom_tile(show.legend = FALSE) +
  scale_fill_viridis_c(limits=c(0.48,0.69)) +
  theme_minimal() +
  theme(axis.text.x=element_text(angle = 46, vjust = 1, hjust = 1), 
        plot.background = element_rect(fill = 'white', colour = 'white'))+
  labs(y = "Tested on ", x = "Trained on") +  
  scale_x_discrete(labels=t2T) + 
  scale_y_discrete(labels = t2T)+
  geom_text(aes(label = round(AUROC_woRef, 4)), color = "white", size = 5)+ 
  ggtitle("Data without reference base")
ggp2 = ggplot(ROCs, aes(trainedOn, testedOn, fill= AUROC)) +
  geom_tile() +
  scale_fill_viridis_c(limits=c(0.48,0.69)) +
  theme_minimal() +
  theme(axis.text.x=element_text(angle = 46, vjust = 1, hjust = 1), 
        plot.background = element_rect(fill = 'white', colour = 'white'))+
  labs(y = "Tested on ", x = "Trained on") +  
  scale_x_discrete(labels=t2T) + 
  scale_y_discrete(labels = t2T)+
  geom_text(aes(label = round(AUROC, 4)), color = "white", size = 5)+ 
  ggtitle("Data with reference base")
legend = get_legend(ggp2)
ggp2 = ggp2 + theme(legend.position="none")

grid.arrange(ggp1, ggp2, legend, ncol = 3, widths=c(2.3, 2.3, 1))  

g <- arrangeGrob(ggp1, ggp2, legend, ncol = 3, widths=c(2.3, 2.3, 1))
ggsave(filename="fig/CrossTissue_ComparisonwoRef_191022.png", g,  
       width = 24, height = 12, units = "cm")
  

library(tidyr)
library(readxl)
library(ggplot2)
library(gridExtra)
library(tidyverse)
source("./scripts/05_analysis/00_NamesAndColors.R")
dir.create("/cellnet/MutationModel/fig/RangeAnalysis/", showWarnings = F)
ranges = c("1Mb" = "500kb", "100kb" = "50kb","10kb" = "5kb", "1kb" =  "500", 
           "100bp" = "50bp", "10bp" = "5bp", "none" = "0bp") # Range to the left and right around the pos
rarePredictors = c("ZBTB33", "YY1", "TAF1", "SP1",
                   "RXRA", "REST", "RAD21",
                   "NR2F2", "MAX", "JUND", "HNF4G",
                   "HNF4A", "GABPA", "FOXA2", 
                   "FOXA1", "EGR1", "ATF3")



# create function to plot the rank comparisons #####
plotComparison = function(method = c("Pearson", "Spearman", "GLM", "IndGLM"), 
                          plotName = NA, 
                          tissues, ranges){
  if(is.na(plotName)){plotName = method}
  plotDat = lapply(tissues, function(t){
    print(t)
    valueTable = lapply(names(ranges), function(r){
      load(file=paste0("data/procData/traindata/traindata_processed_",
                       t, "_BinAnalysis_",r,".RData"))
      dat_bin = dat_bin[datchroms_bin == "chr1",]
      features_range = paste(features, r, sep = "_")  
      if(method == "Pearson"){
        value = suppressWarnings(abs(cor(dat_bin[features_range], 
                                        as.numeric(paste(dat_bin$mutated))))) 
      } else if(method == "Spearman"){
        value = suppressWarnings(abs(cor(dat_bin[features_range], 
                                 as.numeric(paste(dat_bin$mutated)),
                                 method = "spearman"))) 
      } else if(method == "GLM"){
        logR = glm(formula = mutated ~ ., data = dat_bin, 
                   family = binomial(link = "logit"))
        p = summary(logR)$coefficients[,"Pr(>|z|)"]
        value = p[features_range]
      } else if (method == "IndGLM"){
        value = sapply(features_range, function(feat){
          logR = glm(formula = dat_bin$mutated ~ dat_bin[,feat], 
                     family = binomial(link = "logit"))
          if(is.na(logR$coefficients[2])){return(NA)}
          summary(logR)$coefficients["dat_bin[, feat]","Pr(>|z|)"]
        })
      }
      res = data.frame(features = factor(features, levels = features),
                       range = r, value = unname(value))
      return(res)
    })  
    valueTable = do.call(rbind, valueTable)
    valueTable$range = factor(valueTable$range, levels = names(ranges),labels = ranges)
    valueTable$tissue = t2T[t]
    return(valueTable)
  })
  plotDat = do.call(rbind, plotDat)
  plotDat = plotDat[!plotDat$features %in% rarePredictors,]
  # get the rank within each tissue for each predictor 
  plotDat$ranks = NA
  for(feat in unique(plotDat$features)){
    for(tis in unique(plotDat$tissue)){
      if(method %in% c("Pearson", "Spearman")){
        plotDat[plotDat$features == feat & 
                  plotDat$tissue == tis,"ranks"] = 
          rank(-plotDat[plotDat$features == feat & 
                          plotDat$tissue == tis,"value"],
               na.last = "keep")
      } else{
        plotDat[plotDat$features == feat & 
                  plotDat$tissue == tis,"ranks"] = 
          rank(plotDat[plotDat$features == feat & 
                          plotDat$tissue == tis,"value"],
               na.last = "keep")
      }
      
    }
  }
  # and plot it
  {# first the raw values
  ggplot(plotDat, aes(range, tissue)) +
    geom_tile(aes(fill = value)) + 
    ylab("Tissue") + 
    xlab("Range") +
    scale_fill_gradient(low="white", high="red")+
    theme_minimal() + 
    facet_wrap(~features,  strip.position = "top") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=8),
          axis.text.y = element_text(size=8),
          strip.clip = "off",
          strip.text = element_text(hjust = 1, size=10),
          legend.key.size = unit(10,"pt"),
          legend.text = element_text(size=7),
          legend.title = element_text(size=7),
          legend.box.spacing = unit(0,"pt"))
  ggsave(filename = paste0("fig/RangeAnalysis/",plotName,"Values.png"), 
         width = 25, height = 25, units = "cm")
  # then with correlations summarized (median) over tissues
  plotDat %>% group_by(features, range) %>%
    summarize(medianValue = median(value, na.rm = T)) %>%
    ungroup()   %>%
    ggplot(aes(range, features)) +
    geom_tile(aes(fill = medianValue))+
    ylab("Feature") + 
    xlab("Range") +
    scale_fill_gradient(low="white", high="red")+
    theme(legend.key.size = unit(10,"pt"),
          legend.text = element_text(size=7),
          legend.title = element_text(size=7),
          legend.box.spacing = unit(0,"pt"))
  ggsave(filename = paste0("fig/RangeAnalysis/",plotName,"ValueSummary.png"), 
         width = 15, height = 25, units = "cm")
  # then the ranks
  ggplot(plotDat, aes(range, tissue)) +
    geom_tile(aes(fill = ranks)) + 
    ylab("Tissue") + 
    xlab("Range") +
    scale_fill_gradient(low="white", high="red")+
    theme_minimal() + 
    facet_wrap(~features,  strip.position = "top") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=8),
          axis.text.y = element_text(size=8),
          strip.clip = "off",
          strip.text = element_text(hjust = 1, size=10),
          legend.key.size = unit(10,"pt"),
          legend.text = element_text(size=7),
          legend.title = element_text(size=7),
          legend.box.spacing = unit(0,"pt"))
  ggsave(filename = paste0("fig/RangeAnalysis/",plotName,"Ranks.png"), 
         width = 25, height = 25, units = "cm")
  # and the ranks summarized (median)
  plotDat %>% group_by(features, range) %>% 
    summarize(medianRank = median(ranks, na.rm = T)) %>%
    ungroup()   %>%
    ggplot(aes(range, features)) +
    geom_tile(aes(fill = medianRank))  +
    ylab("Features") + 
    xlab("Range") +
    scale_fill_gradient(low="white", high="red")+
    theme(legend.key.size = unit(10,"pt"),
          legend.text = element_text(size=7),
          legend.title = element_text(size=7),
          legend.box.spacing = unit(0,"pt"))
  ggsave(filename = paste0("fig/RangeAnalysis/",plotName,"RankSummary.png"), 
         width = 15, height = 25, units = "cm")
  }
}
#####


# Create different plots #####
plotComparison(method = "Spearman", tissues = tissues, ranges = ranges)
plotComparison(method = "IndGLM", tissues = tissues, ranges = ranges)
plotComparison(method = "Pearson", tissues = tissues, ranges = ranges)
plotComparison(method = "GLM", tissues = tissues, ranges = ranges)
#####

#####
# tissueCorrelations = lapply(tissues, function(t){
#   print(t)
#   correlations = lapply(names(ranges), function(r){
#     load(file=paste0("data/procData/traindata/traindata_processed_",
#                      t, "_BinAnalysis_",r,".RData"))
#     dat_bin = dat_bin[datchroms_bin == "chr1",]
#     features_range = paste(features, r, sep = "_")  
#     cors = suppressWarnings(abs(cor(dat_bin[features_range], 
#                    as.numeric(paste(dat_bin$mutated))))) # Calc pearson correlation coefficient
#     res = data.frame(features = factor(features), correlation = cors, range = r)
#     return(res)
#   })  
#   correlations = do.call(rbind, correlations)
#   correlations$range = factor(correlations$range, levels = names(ranges),labels = ranges)
#   correlations$tissue = t
#   return(correlations)
# })
# tissueCorrelations = do.call(rbind, tissueCorrelations)
# tissueCorrelations$tissue = t2T[tissueCorrelations$tissue]
# tissueCorrelations = tissueCorrelations[!tissueCorrelations$features %in% rarePredictors,]
# # get the rank of correlation within each tissue for each predictor 
# tissueCorrelations$ranks = NA
# for(feat in unique(tissueCorrelations$features)){
#   for(tis in unique(tissueCorrelations$tissue)){
#     tissueCorrelations[tissueCorrelations$features == feat & 
#                          tissueCorrelations$tissue == tis,"ranks"] = 
#       rank(-tissueCorrelations[tissueCorrelations$features == feat & 
#                                   tissueCorrelations$tissue == tis,"correlation"],
#             na.last = "keep")
#   }
# }
# and plot it
# first the raw correlation values
# ggplot(tissueCorrelations, aes(range, tissue)) +
#   geom_tile(aes(fill = correlation)) + 
#   ylab("Features") + 
#   xlab("Range") +
#   scale_fill_gradient(low="white", high="red")+
#   theme_minimal() + 
#   facet_wrap(~features,  strip.position = "top") +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=8),
#         axis.text.y = element_text(size=8),
#         strip.clip = "off",
#         strip.text = element_text(hjust = 1, size=10),
#         legend.key.size = unit(10,"pt"),
#         legend.text = element_text(size=7),
#         legend.title = element_text(size=7),
#         legend.box.spacing = unit(0,"pt"))
# ggsave(filename = paste0("compareCorrelationPerWindowsize_allTissues.png"), 
#        path = "/cellnet/MutationModel/fig/RangeAnalysis/", 
#        width = 25, height = 25, units = "cm")
# # then with correlations summarized (median) over tissues
# tissueCorrelations %>% group_by(features, range) %>%
#   summarize(medianCorrelation = median(correlation, na.rm = T)) %>%
#   ungroup()   %>%
#   ggplot(aes(range, features)) +
#   geom_tile(aes(fill = medianCorrelation))+
#   theme(legend.key.size = unit(10,"pt"),
#         legend.text = element_text(size=7),
#         legend.title = element_text(size=7),
#         legend.box.spacing = unit(0,"pt"))
# ggsave(filename = paste0("compareMedianCorrelationPerWindowsize_allTissues.png"), 
#        path = "/cellnet/MutationModel/fig/RangeAnalysis/", 
#        width = 25, height = 25, units = "cm")
# # then the ranks
# ggplot(tissueCorrelations, aes(range, tissue)) +
#   geom_tile(aes(fill = ranks)) + 
#   ylab("Features") + 
#   xlab("Range") +
#   scale_fill_gradient(low="white", high="red")+
#   theme_minimal() + 
#   facet_wrap(~features,  strip.position = "top") +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=8),
#         axis.text.y = element_text(size=8),
#         strip.clip = "off",
#         strip.text = element_text(hjust = 1, size=10),
#         legend.key.size = unit(10,"pt"),
#         legend.text = element_text(size=7),
#         legend.title = element_text(size=7),
#         legend.box.spacing = unit(0,"pt"))
# ggsave(filename = paste0("compareCorrelationRanksPerWindowsize_allTissues.png"), 
#        path = "/cellnet/MutationModel/fig/RangeAnalysis/", 
#        width = 25, height = 25, units = "cm")
# # and the ranks summarized (median)
# tissueCorrelations %>% group_by(features, range) %>% 
#   summarize(medianRank = median(ranks, na.rm = T)) %>%
#   ungroup()   %>%
#   ggplot(aes(range, features)) +
#   geom_tile(aes(fill = medianRank))  +
#   theme(legend.key.size = unit(10,"pt"),
#         legend.text = element_text(size=7),
#         legend.title = element_text(size=7),
#         legend.box.spacing = unit(0,"pt"))
# ggsave(filename = paste0("compareMedianCorrelationRanksPerWindowsize_allTissues.png"), 
#        path = "/cellnet/MutationModel/fig/RangeAnalysis/", 
#        width = 25, height = 25, units = "cm")
######

# # Compare spearman correlation ####
# tissueSpearman = lapply(tissues, function(t){
#   print(t)
#   
#   correlations = sapply(names(ranges), function(r){
#     load(file=paste0("data/procData/traindata/traindata_processed_",
#                      t, "_BinAnalysis_",r,".RData"))
#     features_range = paste(features, r, sep = "_")  
#     dat_bin = dat_bin[datchroms_bin == "chr1",]
#     
#     cors = suppressWarnings(abs(cor(dat_bin[features_range], 
#                                     as.numeric(paste(dat_bin$mutated)),
#                                     method = "spearman"))) 
#     return(setNames(cors, features))
#   }, simplify = T) 
#   temp = data.frame(features = row.names(correlations), correlations, check.names = F)
#   longCorrelations = pivot_longer(temp, cols = !features, 
#                                   names_to = "range", values_to = "correlation")
#   longCorrelations$features = factor(longCorrelations$features, 
#                                      levels = rownames(correlations))
#   longCorrelations$range = factor(longCorrelations$range, levels = names(ranges))
#   longCorrelations$tissue = t
#   return(longCorrelations)
# })
# 
# # create a figure with all tissues
# tissueSpearman = do.call(rbind, tissueSpearman)
# tissueSpearman$tissue = t2T[tissueSpearman$tissue]
# tissueSpearman = tissueSpearman[!tissueSpearman$features %in% rarePredictors,]
# tissueSpearman$ranks = NA
# for(feat in unique(tissueSpearman$features)){
#   for(tis in unique(tissueSpearman$tissue)){
#     tissueSpearman[tissueSpearman$features == feat & 
#                      tissueSpearman$tissue == tis,"ranks"] = 
#       rank(- tissueSpearman[tissueSpearman$features == feat & 
#                               tissueSpearman$tissue == tis,"correlation"],
#            na.last = "keep")
#   }
# }
# tissueSpearman$range = factor(tissueSpearman$range, 
#                                   levels = names(ranges), 
#                                   labels = ranges)
# ggplot(tissueSpearman) +
#   geom_tile(aes(range, tissue, fill = ranks)) + 
#   ylab("Features") + 
#   xlab("Range") +
#   facet_wrap(~features, scales = "free", strip.position = "top") +
#   scale_fill_gradient(low="white", high="red")+
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=8),
#         axis.text.y = element_text(size=8),
#         strip.clip = "off",
#         strip.text = element_text(hjust = 1, size=10))
# ggsave(filename = paste0("compareSpearmanPerWindowsize_allTissues.png"), 
#        path = "/cellnet/MutationModel/fig/RangeAnalysis/", 
#        width = 25, height = 25, units = "cm")
# 
# #####

# ## Compare GLM p-values per range #####
# # the warning message "fitted probabilities numerically 0 or 1 occurred" can
# # be ignored - it is due to few available levels for some of the predictors
# # those get thrown out when we filter for p.values
# tissueImportances = lapply(tissues, function(t){
#   print(t)
#   importances = lapply(names(ranges), function(r){
#     load(file=paste0("data/procData/traindata/traindata_processed_",
#                      t, "_BinAnalysis_",names(ranges[r]),".RData"))
#     features_range = setNames(features, paste(features, r, sep = "_")) 
#     dat_bin = dat_bin[datchroms_bin == "chr1",]
#     logR = glm(formula = mutated ~ ., data = dat_bin, 
#                family = binomial(link = "logit"))
#     p = summary(logR)$coefficients[,"Pr(>|z|)"]
#     p = p[names(p) %in% names(features_range)]
#     res = data.frame(features = features_range[names(p)],
#                      range = r,
#                      pvalue = p)
#     res$features = factor(res$features, levels = features)
#     return(res)
#   }) 
#   importances=as.data.frame(do.call(rbind, importances))
#   importances$tissue = t
#   return(importances)
# })
# tissueImportances = do.call(rbind, tissueImportances)
# tissueImportances$range = factor(tissueImportances$range, 
#                                  levels = names(ranges), 
#                                  labels = ranges)
# tissueImportances$tissue = t2T[tissueImportances$tissue]
# tissueImportances = tissueImportances[!tissueImportances$features %in% rarePredictors,]
# tissueImportances$ranks = NA
# for(feat in unique(tissueImportances$features)){
#   for(tis in unique(tissueImportances$tissue)){
#     tissueImportances[tissueImportances$features == feat & 
#                         tissueImportances$tissue == tis,"ranks"] = 
#       rank( tissueImportances[tissueImportances$features == feat & 
#                                 tissueImportances$tissue == tis,"pvalue"],
#             na.last = "keep")
#   }
# }
# ggplot(tissueImportances) +
#   geom_tile(aes(range, tissue, fill = ranks)) + 
#   ylab("Features") + 
#   xlab("Range") +
#   facet_wrap(~features, scales = "free", strip.position = "top") +
#   scale_fill_gradient(low="white", high="red")+
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=8),
#         axis.text.y = element_text(size=8),
#         strip.clip = "off",
#         strip.text = element_text(hjust = 1, size=10))
# 
# ggsave(filename = paste0("comparePvaluePerWindowsize_allTissues.png"), 
#        path = "/cellnet/MutationModel/fig/RangeAnalysis/", 
#        width = 25, height = 25, units = "cm")
#####

# ## Compare individual GLM p-values per range #####
# #GCcontent   1Mb       0  Colon   3.0 individual GLM p-value
# tissueIndPvals = lapply(tissues, function(t){
#   print(t)
#   # load once to get predictors
#   load(file=paste0("data/procData/traindata/traindata_processed_",
#                    t, "_BinAnalysis_100bp.RData"))
#   importances = lapply(names(ranges), function(r){
#     load(file=paste0("data/procData/traindata/traindata_processed_",
#                      t, "_BinAnalysis_",r,".RData"))
#     dat_bin = dat_bin[datchroms_bin == "chr1",]
#     features_range = setNames(features, paste(features, r, sep = "_")) 
#     indPvals = sapply(names(features_range), function(feat){
#       logR = glm(formula = dat_bin$mutated ~ dat_bin[,feat], 
#                  family = binomial(link = "logit"))
#       if(is.na(logR$coefficients[2])){return(NA)}
#       summary(logR)$coefficients["dat_bin[, feat]","Pr(>|z|)"]
#     })
#     data.frame(features = factor(features_range[names(indPvals)], 
#                                        levels = features_range),
#                      range = r,
#                      pvalue = indPvals)
#   }) 
#   importances=as.data.frame(do.call(rbind, importances))
#   importances$tissue = t
#   return(importances)
# })
# tissueIndPvals = do.call(rbind, tissueIndPvals)
# tissueIndPvals$pvalue[tissueIndPvals$pvalue == 0] = 2e-16
# tissueIndPvals$range = factor(tissueIndPvals$range, 
#                                  levels = names(ranges), 
#                                  labels = ranges)
# tissueIndPvals$tissue = t2T[tissueIndPvals$tissue]
# tissueIndPvals = tissueIndPvals[!tissueIndPvals$features %in% rarePredictors,]
# tissueIndPvals$ranks = NA
# for(feat in unique(tissueIndPvals$features)){
#   for(tis in unique(tissueIndPvals$tissue)){
#     tissueIndPvals[tissueIndPvals$features == feat & 
#                      tissueIndPvals$tissue == tis,"ranks"] = 
#       rank( tissueIndPvals[tissueIndPvals$features == feat & 
#                              tissueIndPvals$tissue == tis,"pvalue"],
#             na.last = "keep")
#   }
# }
# ggplot(tissueIndPvals) +
#   geom_tile(aes(range, tissue, fill = ranks)) + 
#   ylab("Features") + 
#   xlab("Range") +
#   facet_wrap(~features, scales = "free", strip.position = "top") +
#   scale_fill_gradient(low="white", high="red")+
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=8),
#         axis.text.y = element_text(size=8),
#         strip.clip = "off",
#         strip.text = element_text(hjust = 1, size=10))
# ggsave(filename = paste0("compareIndivPvaluePerWindowsize_allTissues.png"), 
#        path = "fig/RangeAnalysis/", 
#        width = 25, height = 25, units = "cm")
#####

# # for each predictor, plot the four metrics together #####
# cor = cbind(tissueCorrelations,
#             method = "Pearson correlation")
# colnames(cor)[colnames(cor) == "correlation"] = "measure"
# spear = cbind(tissueSpearman,
#               method = "Spearman correlation")
# colnames(spear)[colnames(spear) == "correlation"] = "measure"
# pval = cbind(tissueImportances,
#              method = "multiple GLM p-value")
# colnames(pval)[colnames(pval) == "pvalue"] = "measure"
# indPval =  cbind(tissueIndPvals,
#                  method = "individual GLM p-value")
# colnames(indPval)[colnames(indPval) == "pvalue"] = "measure"
# dat = rbind(cor, spear, pval, indPval)
# save(dat, file = "temp/temp.RData")
# feats = unique(dat$features)
# dumpVar = sapply(feats, function(feat){
#   dumpVar2 = sapply(unique(dat$method), function(meth){
#     subDat = dat[dat$features == feat & dat$method == meth,]
#     if(meth %in% c("multiple GLM p-value" ,"individual GLM p-value")){
#       subDat$measure = -log10(subDat$measure)
#     }
#     dumpVar3 = lapply(c("measure", "ranks"), function(toPlot){
#       # get the range with the best average value
#       if(toPlot == "measure"){
#         best = names(which.max(sapply(as.character(unique(subDat$range)), function(ran){
#           median(subDat[subDat$range == ran,toPlot], na.rm = T)
#         })))
#       } else{
#         best = names(which.min(sapply(as.character(unique(subDat$range)), function(ran){
#           median(subDat[subDat$range == ran,toPlot], na.rm = T)
#         })))
#       }
#       
#       p1 = ggplot(subDat, aes(tissue, range)) +
#         geom_tile(aes(fill = .data[[toPlot]])) + 
#         ylab("Range") +
#         xlab("Tissue") +
#         scale_fill_gradient(low="white", high="red")+
#         theme_minimal() +
#         theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=7),
#               axis.text.y = element_text(size=7),
#               strip.clip = "off",
#               strip.text = element_text(hjust = 0, size=10),
#               strip.text.y = element_text(angle = 0),
#               legend.key.size = unit(10,"pt"),
#               legend.text = element_text(size=7),
#               legend.title = element_text(size=7),
#               legend.box.spacing = unit(0,"pt"))
#       temp = p1 +  geom_tile(data = subDat %>% filter(range == best),
#                   fill = NA, colour = "black", linewidth = 2)
#       markingBox = layer_data(temp, 2) %>%
#         summarise(xmin = min(xmin), xmax = max(xmax),
#                   ymin = min(ymin), ymax = max(ymax))
#       p1 = p1 + geom_rect(data = markingBox,
#                 aes(xmin = xmin, xmax = xmax,
#                     ymin = ymin, ymax = ymax),
#                 inherit.aes = FALSE, fill = NA, 
#                 colour = "black", linewidth = 1)
#       return(p1)
#     })
#     grid.arrange(grobs = dumpVar3, nrow = 1, top = as.character(meth))
#   }) 
#   grid.arrange(grobs = dumpVar2, nrow = 2, top = as.character(feat))
# })
# multiPage = marrangeGrob(grobs = dumpVar, nrow = 1, ncol = 1, top = "")
# ggsave("fig/RangeAnalysis/compareMeasuresPerWindowsize_allTissues.pdf",
#        multiPage)
# # ggsave("fig/RangeAnalysis/compareMeasuresPerWindowsize_allTissues.png", 
#        # width = 25, height = 100, units = "cm")
#####


# ## Compare RF importances ######
# for(t in tissues){
#   print(t)
#   
#   # Load mapped data once to get the features of interest
#   load(file=paste0("data/procData/traindata/traindata_processed_",
#                      t, "_BinAnalysis_100bp.RData"))
#   
#   importances = sapply(seq(length(ranges),1), function(r){
#     features_range = paste(features, names(ranges[r]), sep = "_")  
#       
#     load(file = paste0("data/rdata/RFmodel/BinAnalysis/", t, "_",names(ranges[r]),
#                       "_importances.RData"))
#     
#     imp_ranges = imp[rownames(imp) %in% features_range,]
#     return(rowMeans(imp_ranges))
#   }) 
#   colnames(importances)= names(ranges)
#   importances = cbind(features, importances) %>% as.data.frame()
#   
#   importances_gatherd <- gather(data = importances, key = range, value = importance, 1:length(ranges)+1)
#   importances_gatherd$range = factor(importances_gatherd$range, levels = rev(names(ranges)))
#   importances_gatherd$importance = as.numeric(importances_gatherd$importance)
#   
#   ggplot(importances_gatherd, aes(range, features, fill = importance)) +
#     geom_tile() +
#     scale_fill_gradient(low="white", high="navy", name = "Gini importance")+
#     theme_minimal()+
#     ylab("Features") + 
#     xlab("Range")
#   
#   ggsave(filename = paste0("compareRFimpPerWindowsize_",t,".png"), path = "/cellnet/MutationModel/fig/RangeAnalysis/", 
#          width = 15, height = 15, units = "cm")
# }
#####

# ## Compare RF importances of all features #####
# for(t in tissues){
#   print(t)
#   
#   # Load mapped data once to get the features of interest
#   importances = sapply(seq(length(ranges),1), function(r){
#     load(file = paste0("data/rdata/RFmodel/BinAnalysis/", t, "_",names(ranges[r]),
#                        "_importances.RData"))
#     
#     #imp_ranges = imp[rownames(imp) %in% features_range,]
#     return(rowMeans(imp))
#   }) 
#   colnames(importances)= names(ranges)
#   importances = cbind(features = rownames(importances), importances) %>% as.data.frame()
#   importances$features = as.character(strsplit(importances$features, "_100bp"))
#   
#   importances_gatherd <- gather(data = importances, key = range, value = importance, 1:length(ranges)+1)
#   importances_gatherd$range = factor(importances_gatherd$range, levels = rev(names(ranges)))
#   importances_gatherd$importance = as.numeric(importances_gatherd$importance)
#   
#   ggplot(importances_gatherd, aes(range, features, fill = importance)) +
#     geom_tile() +
#     scale_fill_gradient(low="white", high="navy", name = "Gini importance")+
#     theme_minimal()+
#     ylab("Features") + 
#     xlab("Range")
#   
#   ggsave(filename = paste0("compareRFimpPerWindowsize_",t,"_allFeatures.png"), path = "/cellnet/MutationModel/fig/RangeAnalysis/", 
#          width = 15, height = 15, units = "cm")
#   
# }
# #####



# ### Feature correlation between tissue spec. features #####
# tissue_comb = expand.grid(tissues, tissues, stringsAsFactors = F)
# tissue_comb = tissue_comb[!duplicated(t(apply(tissue_comb, 1, sort))),]
# 
# 
# # Getting list of all tissue-specific features
# features = sapply(tissues, function(t){
#   tab = read_xlsx("data/rawdata/dataMapping_BinAnalysis.xlsx", 
#                   sheet=t, col_names=T)
#   return(tab$abbreviation[tab$sourceTissue == t & tab$range == "x"])
# }) %>% unlist() %>% unname() %>% unique()
# features = features[!features %in% rarePredictors]
# 
# 
# # Get correlation between features for each tissue combination
# corTab = sapply(tissues, function(t1){
#   load(paste0("data/procData/traindata/traindata_processed_",
#               t1, ".RData"))
#   dat1 = dat
#   tab1 = read_xlsx("data/rawdata/dataMapping_BinAnalysis.xlsx", 
#                   sheet=t1, col_names=T)
#   subT = sapply(tissues, function(t2){
#     if(t1 == t2){
#       return(NULL)
#     } else{
#       load(paste0("data/procData/traindata_crosstissue/traindata_crosstissue_",
#                   t1, "With",t2,"Preds_processed", ".RData"))
#       tab2 = read_xlsx("data/rawdata/dataMapping_BinAnalysis.xlsx", 
#                        sheet=t2, col_names=T)
#       features = tab1$abbreviation[tab1$abbreviation %in% tab2$abbreviation &
#                                      !is.na(tab1$sourceTissue)]
#       present = features[features %in% colnames(dat1) & features %in% colnames(dat)]
#     }
#   })
# })
# cor_list = lapply(features, function(f){
#   print(f)
#   
#   corF = data.frame()
#   
#   for(t1 in tissues){
#     # Data of tissue 1
#     load(paste0("data/procData/traindata/traindata_processed_",
#                 t1, ".RData"))
#     dat1 = dat
#     for(t2 in tissues){
#       if(t1 == t2){
#         c = 1
#       } else {   # Data of tissue 2 mapped on tissue 1 positions
#         load(paste0("data/procData/traindata_crosstissue/traindata_crosstissue_",
#                     t1, "With",t2,"Preds_processed", ".RData"))
#         # Correlation for given data
#         if(is.null(dat1[[f]]) | is.null(dat[[f]])){
#           c = NA
#         } else {
#           c = cor(dat1[[f]], dat[[f]])
#         }
#       }
#       corF = rbind(corF, c(f, t1, t2, c))
#       colnames(corF) = c("feature", "tissue1", "tissue2", "correlation")
#     }
#   }
#   return(corF)
# }) 
# ######
# 
# # For each duplicated correlation (t1 - t2, t2 -t1) take the mean correlation ####
# cor = lapply(seq(length(features)), function(f){
#   print(features[f])
#   data = cor_list[[f]]
#   cor_temp = data.frame()
#   
#   
#   for(i in seq(nrow(tissue_comb))){
#     t1 = as.character(tissue_comb[i,][1])
#     t2 = as.character(tissue_comb[i,][2])
#     
#     if(t1 != t2){
#       temp1 = data[data$tissue1 == t1 & data$tissue2 == t2,]
#       temp2 = data[data$tissue2 == t1 & data$tissue1 == t2,]
#       
#       mean_cor = mean(as.numeric(temp1$correlation), as.numeric(temp2$correlation), na.rm = T)
#       temp1$correlation = mean_cor
#       
#       cor_temp = rbind(cor_temp, temp1)
#       
#     } else {
#       temp = data[data$tissue1 == t1 & data$tissue2 == t2,]
#       cor_temp = rbind(cor_temp, temp)
#     }
#   }
#   
#   return(cor_temp)
# })
# #####
# 
# # Plot heatmap for each feature #####
# for(p in seq(length(features))){
#   ggplot(cor[[p]], aes(tissue1, tissue2, fill = as.numeric(correlation))) +
#     geom_tile() + ggtitle(p2P[cor[[p]][1,1]])+
#     scale_x_discrete(labels=t2T)+
#     scale_y_discrete(labels=t2T)+
#     theme_minimal()+
#     theme(axis.title.x = element_blank(),
#           axis.title.y = element_blank(),
#           axis.ticks = element_blank(),
#           legend.position = c(0.15, 0.8))+
#     scale_fill_gradient(low="white", high="red", name = "Pearson\nCorrelation")
#   
#   ggsave(filename = paste0("allTissue_correlation_", cor[[p]][1,1], ".png"), 
#          path = "/cellnet/MutationModel/fig/featureExchange/", 
#          width = 15, height = 15, units = "cm")
# }
#####

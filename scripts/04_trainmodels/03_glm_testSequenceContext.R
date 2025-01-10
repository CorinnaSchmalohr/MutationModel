# get example data #####
predictorOrder = read.table("scripts/04_trainmodels/predictorOrder_glm.txt")[,1]
tissue = "esophagus"
cr = "chr1"
load(paste0(paste0("data/procData/traindata/traindata_processed_", tissue, ".RData")))
trainData = dat[datchroms != cr,]
testData = dat[datchroms == cr,]
#####


# test if the data is balanced in terms of reference base: #####
table(dat$ref_CG[dat$mutated == 0])
table(dat$ref_CG[dat$mutated == 1])

table(trainData$ref_CG[trainData$mutated == 0])
table(trainData$ref_CG[trainData$mutated == 1])
# yes #####


# method of testing and predictor order #####
# https://stats.stackexchange.com/questions/312137/p-value-from-a-binomial-model-glm-for-a-binomial-predictor
model1 = glm(formula = mutated ~ ., data = trainData, 
             family = binomial(link = "logit"))
res1 = anova(model1, test = "LRT") 
res1 = setNames(res1$`Pr(>Chi)`, rownames(res1))[-1]
# res2 = anova(model1, test = "Chisq") # same as LRT
# res2 =  setNames(res2$`Pr(>Chi)`, rownames(res2))
res3 = summary(model1)$coefficients[,4]
res3 = setNames(res3, gsub("CG1", "CG",names(res3)))
res4 = drop1(object = model1, test = "LRT")
res4 = setNames(res4$`Pr(>Chi)`, rownames(res4))
# second model with reordered predictors
model2 = glm(formula = mutated ~ ., data = trainData[,ncol(trainData):1], 
             family = binomial(link = "logit"))
res5 = anova(model2, test = "LRT") # this is wrong, depends on predictor order
res5 = setNames(res5$`Pr(>Chi)`, rownames(res5))
# res6 = anova(model2, test = "Chisq")
# res6 = setNames(res6$`Pr(>Chi)`, rownames(res6))
res7 = summary(model2)$coefficients[,4]
res7 = setNames(res7, gsub("CG1", "CG",names(res7)))
res8 = drop1(object = model2, test = "LRT")
res8 = setNames(res8$`Pr(>Chi)`, rownames(res8))
comparison = data.frame(anova = res1, 
                        summary = res3[names(res1)],
                        drop1 = res4[names(res1)], 
                        reordered_anova = res5[names(res1)], 
                        reordered_summary = res7[names(res1)],
                        reordered_drop1 = res8[names(res1)])
cols = rep(1,length(res1))
cols[names(res1) %in% c("precedingBase_CG", "ref_CG","followingBase_CG")] = 2:4
png("fig/SequenceContextTestPvalue.png", height = 800, width = 800, pointsize = 20)
plot(comparison, col = cols, main = "p-value")
dev.off()
rm(model1, model2, res1, res2, res3, res4, res5, res6, res7, res8, cols, comparison)
# the method of p-value calculation does not matter. Using anova() depends on predictor order #####


# does the term order in the formula influence the model? ####
# From ?glm
#The terms in the formula will be re-ordered so that main effects come first, 
#followed by the interactions, all second-order, all third-order and so on: 
#to avoid this pass a terms object as the formula.

f1 = terms(mutated ~ ., data = trainData,
          keep.order = TRUE)
model1 = glm(formula = f1, data = trainData, family = binomial("logit"))
f2 = terms(mutated ~ ., data = trainData[,ncol(trainData):1],
           keep.order = TRUE)
model2 = glm(formula = f1, data = trainData, family = binomial("logit"))
all.equal(coef(summary(model1))[,4],coef(summary(model2))  [,4] ) # ändert nichts
rm(f1, f2, model1, model2)
# no #####

# mult. testing correction #####
model = glm(formula = mutated ~ ., data = trainData, 
                family = binomial(link = "logit"))
p.adjust(coef(summary(model))[,4], method = "fdr")
rm(model)
# --> still significant #####



# Is it still significant if I remove all non-significant predictors? #####
model = glm(formula = mutated ~ ., data = trainData, 
            family = binomial(link = "logit"))
pvals = coef(summary(model))[,4]
signPreds = names(which(pvals <= 0.05))[-1]
signPreds = gsub("CG1", "CG", signPreds)
reducedData = trainData[,c("mutated", signPreds)]
redModel = glm(formula = mutated ~ ., data = reducedData, 
               family = binomial(link = "logit"))
summary(redModel)
# yes! #####


# test if removing context influences model performance #####
library(ROCR)
tissues = c("brain","breast", "colon","esophagus","kidney","liver", "lung","ovary",
            "prostate", "skin")
minusContext = sapply(tissues, function(tissue){
  print(tissue)
  load(paste0(paste0("data/procData/traindata/traindata_processed_", tissue, ".RData")))
  chroms = unique(datchroms)
  temp = t(sapply(chroms, function(cr){
    trainData = dat[datchroms != cr,]
    testData = dat[datchroms == cr,]
    
    reducedData = trainData[,!colnames(trainData) %in% c("precedingBase_CG", "ref_CG", "followingBase_CG")]
    
    # train models
    realModel = glm(formula = mutated ~ ., data = trainData, 
                    family = binomial(link = "logit"))
    reducedModel = glm(formula = mutated ~ ., data = reducedData, 
                       family = binomial(link = "logit"))
    
    # get predictions
    yhatReal = predict(realModel, newdata = testData, type = "response")
    yhatReduced = predict(reducedModel, newdata = testData, type = "response")
    predReal =  prediction(yhatReal, testData$mutated)
    # ROCreal = performance(predReal, "tpr", "fpr")
    predReduced = prediction(yhatReduced, testData$mutated)
    # ROCreduced = performance(predReduced, "tpr", "fpr")
    # plot(ROCreal, main = tissue)
    # plot(ROCreduced, add = T, col = "red")
    # abline(0,1, col = "grey", lty = 2)
    aucReal = performance(predReal,"auc")@y.values[[1]]
    aucReduced = performance(predReduced,"auc")@y.values[[1]]
    return(c(fullModel = aucReal, reducedModel = aucReduced))
  }))
  return(temp)
}, simplify = F)

png("fig/SequenceContextRemoved.png", width = 800, height = 700, pointsize = 20)
par(mfrow = c(3,4), mar = c(1,3,3,0), oma = c(5,3,0,0))
sapply(tissues, function(tissue){
  plot(minusContext[[tissue]], main = tissue, las = 1, xlab = "", ylab = "")
  abline(0,1)
})
mtext("Full model", side = 1, outer = T, line = 2)
mtext("reduced model", side = 2, outer = T, line = 1)
dev.off()
# a little bit, for most tissues #####


# TODO Why difference between tissues?#####
# compare colon and liver: which predictors differ?
#####

# is it overfitting? compare test and train error #####
library(ROCR)
tissues = c("brain","breast", "colon","esophagus","kidney","liver", "lung","ovary",
            "prostate", "skin")
testTrain = sapply(tissues, function(tissue){
  print(tissue)
  load(paste0(paste0("data/procData/traindata/traindata_processed_", tissue, ".RData")))
  chroms = unique(datchroms)
  par(mfrow = c(3,3))
  temp = t(sapply(chroms, function(cr){
    trainData = dat[datchroms != cr,]
    testData = dat[datchroms == cr,]
    
    realModel = glm(formula = mutated ~ ., data = trainData, 
                    family = binomial(link = "logit"))
    
    # get predictions
    yhatTrain = realModel$fitted.values
    yhatTest = predict(realModel, newdata = testData, type = "response")
    predTrain =  prediction(yhatTrain, trainData$mutated)
    predTest =  prediction(yhatTest, testData$mutated)
    # ROCtrain = performance(predTrain, "tpr", "fpr")
    # ROCtest = performance(predTest, "tpr", "fpr")
    # plot(ROCtrain)
    # plot(ROCtest, add = T, col = "red")
    # abline(0,1, col = "grey", lty = 2)
    # legend("bottomright", legend = c("train data", "test data"), col = 1:2, lty = 1)
    aucTrain = performance(predTrain,"auc")@y.values[[1]]
    aucTest = performance(predTest,"auc")@y.values[[1]]
    return(c(train = aucTrain, test = aucTest))
  }))
  return(temp)
}, simplify = F)
png("fig/TestvsTrainPerformance.png")
par(mfrow = c(3,4), mar = c(2,3,2,0), oma = c(2,2,0,0))
sapply(tissues, function(tissue){
  plot(testTrain[[tissue]], main = tissue, las = 1)
  abline(0,1)
})
mtext(text = "training data AUC",side = 1, line = 0.5, outer = T)
mtext(text = "test data AUC",side = 2, line = 0.5, outer = T)
dev.off()
rm(testTrain)
#####


# is it because of sequence context as a factor/binary predictor? #####
tempData = trainData
tempData$ref_CG = as.numeric(tempData$ref_CG)  
model1 = glm(formula = mutated ~ ., data = tempData, 
            family = binomial(link = "logit"))
summary(model1)
model2 = glm(formula = mutated ~ ., data = trainData, 
            family = binomial(link = "logit"))
summary(model2)

# add balanced test predictor to trainData
tempData$test = NA
tempData$test[tempData$mutated == 0] = c(0,1) 
tempData$test[tempData$mutated == 1] = c(0,1)
table(tempData$test[tempData$mutated == 0])
table(tempData$test[tempData$mutated == 1])
logRT = glm(formula = mutated ~ ., data = tempData, 
            family = binomial("logit"))
summary(logRT)
# Test variable is not significant but ref still is


# permute sequence context 
tempData = trainData
tempData$ref_CG = sample(tempData$ref_CG)  
model1 = glm(formula = mutated ~ ., data = tempData, 
             family = binomial(link = "logit"))
summary(model1)
# no longer significant


# center and scale
tempData = trainData
tempData[,c("precedingBase_CG", "ref_CG", "followingBase_CG")] = 
  apply(tempData[,c("precedingBase_CG", "ref_CG", "followingBase_CG")], 2,function(x){
  x = as.numeric(x)
  x = (x-min(x, na.rm = T))/(max(x, na.rm = T)-min(x, na.rm = T)) # scale from 0 to 1
  x = log(x +1) # log
  x[is.na(x)] = mean(x,na.rm = T) # replace missing
  return(x)
})
model1 = glm(formula = mutated ~ ., data = tempData, 
             family = binomial(link = "logit"))
summary(model1)
# still significant

rm(tempData, model1, model2, model, temp, tempData, logRT)
#####


# demonstrate the dependence on number of other predictors #####
# add more and more predictors to the model
nPreds = sapply(colnames(dat)[53:1], function(pred){
  i = which(colnames(dat) == pred)
  tempDat= dat[datchroms != cr,i:ncol(dat)]
  model = glm(formula = mutated ~ ., data = tempDat, 
              family = binomial(link = "logit"))
  pvals = summary(model)$coefficients[,4]
  return(c(pvals[c("precedingBase_CG1", "ref_CG1", "followingBase_CG1")]))
})
# calculate the difference in refBase p-value that adding each predictor caused
steps = diff(-log(nPreds[2,]))
# and get the p-value of the predictors in the full model
fullModel = glm(formula = mutated ~ ., data = trainData, 
                family = binomial(link = "logit"))
pvals = coef(summary(fullModel))[names(steps),4]

# removing predictors one by one
res = sapply(colnames(dat)[1:52], function(pred){
  i = which(colnames(dat) == pred)
  tempDat= dat[datchroms != cr,-i]
  model = glm(formula = mutated ~ ., data = tempDat, 
              family = binomial(link = "logit"))
  pvals = summary(model)$coefficients[,4]
  return(pvals[c("precedingBase_CG1", "ref_CG1", "followingBase_CG1")])
})

# calculate correlation
tempData= as.data.frame(apply(dat,2,as.numeric))
corr = cor(tempData$ref_CG, tempData[,1:52])[1,]

# visualize
png("fig/SequenceContextOtherPredictors.png", 
    width = 1200, height = 900, pointsize = 25)
layout(cbind(c(1,1),c(2,3), c(4,4)))
par(mar = c(4,8,1,1))
temp = barplot(-log(nPreds[2,]), horiz = T, las = 1,cex.names = 0.7, 
               xlab = "-log(p-value)")
abline(h = temp[which(steps>1)]+0.5, col = 2:5, lty =1, lwd = 2)
cols = setNames(rep(1,length(steps)), names(steps))
cols[steps > 1] = 2:5
par(mar = c(4,4,1,1))
plot(-log(pvals[names(steps)]), steps, xlab = "predictor -log(p-value) in full model", 
     ylab = "increase in -log(p-value)\n of referene Base", las = 1, 
     col = cols, mgp = c(2,1,0))
plot(corr^2, steps, xlab = "R² with reference Base", 
     ylab = "increase in -log(p-value)\n of reference Base", 
     col = cols, mgp = c(2,1,0))
par(mar = c(4,8,1,1))
barplot(res[2,names(pvals)],  las = 1, xlab = "refBase p-value",
     col = cols, mgp = c(3,1,0), horiz = T, cex.names = 0.7)
dev.off()

rm(cols, corr, nPreds, fullModel, pvals, res, steps, temp, tempData)
#####


# check correlation between predictors #####
tempData= as.data.frame(apply(dat,2,as.numeric))
corr = cor(tempData)
toRemove = which(corr["ref_CG",]^2 > 0.01)
toRemove = toRemove[names(toRemove) != "ref_CG"]
reducedData = trainData[,!(colnames(trainData) %in% toRemove)]
reducedModel = glm(formula = mutated ~ ., data = reducedData, 
                   family = binomial(link = "logit"))
summary(reducedModel)

rm(cols, corr, fullModel, nPreds, pvals, reducedData, reducedModel, 
   res, steps, temp, tempData, toRemove)
######




# TODO does removing preceding and following have an  impact on sig. of ref? #####
tempData = trainData[-c(43,45)]
tempModel = glm(formula = mutated ~ ., data = tempData, 
                   family = binomial(link = "logit"))
fullModel = glm(formula = mutated ~ ., data = trainData, 
                family = binomial(link = "logit"))
coef(summary(tempModel))["ref_CG1",]
coef(summary(fullModel))["ref_CG1",]
# no (as expected) ######

# allow for interaction #####
normalModel = glm(formula = mutated ~ ., data = trainData, 
                family = binomial(link = "logit"))
intModel = glm(formula = mutated ~ . + ref_CG*., data = trainData, 
               family = binomial(link = "logit"))
normalPvals = coef(summary(normalModel))[,4]
temp = coef(summary(intModel))[,4]
singlePvals = temp[-grep(":", names(temp))]
intPvals = temp[grep(":", names(temp))]
names(intPvals) = sapply(names(intPvals), function(x){
  x = strsplit(x, split = ":")[[1]]
  return(x[x!="ref_CG1"])
})
# TODO, this plot has to be finalized
png("fig/ContextInteraction.png", height = 1200, width = 600, pointsize = 25)
library(gplots)
compar = rbind(normalPvals = (normalPvals<0.05)+0,
               newPvals = (singlePvals[names(normalPvals)]<0.05)+0,
               intPvals = (intPvals[names(normalPvals)]<0.05)+0)
htmap = heatmap.2(t((compar)), dendrogram = "none", trace = "none", key = F, breaks = 3,
                  density.info = "none", 
                  lhei = c(1,10), lwid = c(0.5,10), 
                  Rowv = F, Colv = F,key.par=list(mar = c(1,1,1,1)),
                  sepwidth=c(0.001,0.01),sepcolor="black",
                  colsep=1:nrow(compar),rowsep=1:ncol(compar),
                  margin=c(5, 10), cexCol = 1, srtCol = 45, na.color = "grey")
legend("top", fill = c("yellow", "red"), legend = c("yes", "no"), 
       title = "significant", xpd = NA, inset = -0.12, ncol = 2, y.intersp = 0.8)
dev.off()

# model without the features that interact
toRemove = names(which(intPvals < 0.1))
reducedData = trainData[,!colnames(trainData) %in% toRemove]
reducedModel = glm(formula = mutated ~ ., data = reducedData, 
                   family = binomial(link = "logit"))
summary(reducedModel)
#####




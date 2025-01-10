args = as.numeric(commandArgs(trailingOnly=T))
tissues = c("brain","breast", "colon","esophagus","kidney","liver", "luad","ovary",
            "prostate", "skin")
tissue = tissues[args]
print(tissue)
dir.create("data/rdata/GLMmodel/",showWarnings=F)




# prepare data for this tissue######
load(paste0(paste0("data/procData/traindata_processed_", tissue, ".RData")))
chroms = unique(datchroms)
#####


# glm #####
print("training models")
temp = lapply(chroms, function(cr){
   cat(cr, ' ')
   trainData = dat[datchroms != cr,]
   logR = glm(formula = mutated ~ ., data = trainData, 
              family = binomial(link = "logit"))
   save(logR, file = paste0("data/rdata/GLMmodel/", tissue, "_", cr, ".RData"))
   return(NULL)
})
cat('\n')
#####



# get variable importances ####
print("var Importances")
imp = sapply(chroms, function(cr){
   cat(cr, ' ')
   load(paste0("data/rdata/GLMmodel/", tissue, "_", cr, ".RData"))
   return(logR$coefficients)
})
save(imp, file = paste0("data/rdata/GLMmodel/", tissue,
                        "_importances.RData"))
cat('\n')
#####


# predict on held-out chromosomes #####
print("predictions")
predictions = lapply(chroms, function(cr){
   cat(cr, ' ')
   load(paste0("data/rdata/GLMmodel/", tissue, "_", cr, ".RData"))
   testData = dat[datchroms == cr,]
   yhat = predict(logR, newdata = testData, type = "response")
   temp = data.frame(pred  = yhat,label = testData$mutated)
   return(temp)
})
names(predictions) = chroms
save(predictions, file = paste0("data/rdata/GLMmodel/", tissue,
                                "_predictions.RData"))
cat('\n')
#####


# calculate variable  p-values ####
pvals <- sapply(chroms, function(cr){
   load(paste0("data/rdata/GLMmodel/", tissue, "_", cr, ".RData")) 
   drop1_features <- drop1(object = logR, test = "LRT")$`Pr(>Chi)`
   return(drop1_features)
})
pvals <- pvals[-1,] # remove first empty row
pvals <- as.data.frame(pvals)
save(pvals, file = paste0("data/rdata/GLMmodel/", tissue,
                          "_pvals.RData"))
#####

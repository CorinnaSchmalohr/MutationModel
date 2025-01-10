setwd("/cellnet/MutationModel/")

### Packages ###
library(dplyr)
library(data.table)
library(tibble)
library(ROCR)
library(ggplot2)
library(gridBase)
library(grid)
library(tidyr)
###############

tissues = c("luad",  "skin", "colon", "ovary", "kidney", "prostate", "breast")

set.seed(1)

##### CHROMOSOME-WISE EVALUATION FOR AN EVEN AMOUNT OF TP&TN #####
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
  
  chroms <- unique(removed$chr)
  
  load(paste0(path, "/mutation_ratio.RData"))
  
  
  
  ##### Data with the same number of TN&TP per chromosome #####
  print("Creating data...")
  
  data_evenTPTN <- lapply(chroms, function(chrom){
    num <- min(mutation_ratio[c("mutated", "non-mutated"),]) # Smallest number of TP or TN per chromosme within the whole data
    
    data_oneChr <- subset(data_wChr, chr == chrom) # Subset data chromosme-wise
    
    data_oneChr_sample0 <- subset(data_oneChr, mutated == "0") %>% rownames_to_column("rowname") %>% 
      sample_n(size = num) %>% column_to_rownames("rowname") # Sample non-mutated (0) data
    
    data_oneChr_sample1 <- subset(data_oneChr, mutated == "1") %>% rownames_to_column("rowname") %>% 
      sample_n(size = num) %>% column_to_rownames("rowname") #Sample mutated (1) data
    
    bound_data <- rbind(data_oneChr_sample0, data_oneChr_sample1)
    
    return(bound_data)
  })
  
  data_evenTPTN <- as.data.frame(do.call(rbind, data_evenTPTN))
  save(data_evenTPTN, file = paste0(path, "/data_evenTPTN.RData"))
  
  
  
  ###### Create glm with leave-one-chromosome-out CV ######
  print("Model training and testing...")
  
  lapply(chroms, function(chrom){
    
    # Data for training, n-1 chromosomes 
    data_train <- subset(data_evenTPTN, chr != chrom, select = -c(chr))
    
    # glm
    logR <- glm(formula = mutated ~ ., data = data_train, 
                family = binomial(link = "logit"))
    save(logR, file = paste0(path, "/", chrom,"/logR_evenTPTN_", chrom,".RData"))
  })
  
  
  
  ##### Test model on remaining chromosome #####
  lapply(chroms, function(chrom){
    load(paste0(path, "/", chrom,"/logR_evenTPTN_", chrom,".RData"))
    
    # Data for testing
    data_test <- subset(data_evenTPTN, chr == chrom, select = -c(chr, mutated))
    
    # Prediction
    yhat <- predict(logR, newdata = data_test, type = "response")
    save(yhat, file = paste0(path, "/", chrom,"/yhat_evenTPTN_",chrom,".RData"))
  })  
  
  
  ##### MODEL PERFORMANCE #####
  print("Model performance...")
  
  # ROC performances per leave-one-chromosome-out CV 
  ROC_CV <- lapply(chroms, function(chrom){ 
    
    load(paste0(path, "/", chrom,"/yhat_evenTPTN_",chrom,".RData"))
    
    data_test <- subset(data_evenTPTN,  chr == chrom)
    
    predict <- prediction(yhat, as.numeric(as.character(data_test$mutated)))
    
    return(performance(predict, "tpr", "fpr"))
  })
  names(ROC_CV) <- chroms
  save(ROC_CV, file = paste0(path, "/ROC_evenTPTN_CV.RData"))
  
  
  # Precision/Recall performances per leave-one-chromosome-out CV
  PR_CV <- lapply(chroms, function(chrom){ 
    
    load(paste0(path, "/", chrom,"/yhat_evenTPTN_",chrom,".RData"))
    
    data_test <- subset(data_evenTPTN,  chr == chrom)
    
    predict <- prediction(yhat, as.numeric(as.character(data_test$mutated)))
    
    return(performance(predict, "prec", "rec"))
  })
  names(PR_CV) <- chroms
  save(PR_CV, file = paste0(path, "/PR_evenTPTN_CV.RData"))
  
  
  
  ##### CONCATENATED RESULTS #####
  print("Concatenate results...")
  
  yhat_concat <- sapply(chroms, function(chrom){
    load(paste0(path, "/", chrom,"/yhat_evenTPTN_",chrom,".RData"))
    return(as.vector(yhat))
  }) %>% unlist() %>% as.data.frame() # All predictions combined
  
  save(yhat_concat, file = paste0(path, "/yhat_evenTPTN_concat.RData"))
  
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
  
  load(paste0(path, "/ROC_evenTPTN_CV.RData"))       # Load ROC
  load(paste0(path, "/PR_evenTPTN_CV.RData"))        # Load PR

  
  ##### Plots #####
  # ROC and PR 
  png(filename = paste0(path_fig, "/ROC_PR_evenTPTN_CV_",tissue,".png"), width = 720*2, height = 720, pointsize = 22)
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
}) 



####### Select tissue #######
tissue = tissues[1]

path = paste0("data/rdata/", tissue, "/logRegression")
path_fig = paste0("fig/", tissue, "/logRegression")
#############################

load(paste0("data/rdata/", tissue, "/completeData_withBinwise.RData"))
  
data$context = NULL
data$pentamer = NULL
data$trimer = NULL
data$septamer = NULL
data$inexon = NULL
  
data_wChr <- cbind(removed$chr, data) # Dataset with corresponding chromosome assignmet
colnames(data_wChr)[1] <- "chr"  

chroms <- unique(removed$chr)

load(paste0(path, "/mutation_ratio.RData"))


###### CHECK FOR DATA SEELCTION RANDOMNESS #####
lapply(1:5, function(i){  
  
  print(i)
  set.seed(i)
  
  ##### Data with the same number of TN&TP per chromosome #####
  data_evenTPTN <- lapply(chroms, function(chrom){
    num <- min(mutation_ratio[c("mutated", "non-mutated"),]) # Smallest number of TP or TN per chromosme within the whole data
    
    data_oneChr <- subset(data_wChr, chr == chrom) # Subset data chromosme-wise
    
    data_oneChr_sample0 <- subset(data_oneChr, mutated == "0") %>% rownames_to_column("rowname") %>% 
      sample_n(size = num) %>% column_to_rownames("rowname") # Sample non-mutated (0) data
    
    data_oneChr_sample1 <- subset(data_oneChr, mutated == "1") %>% rownames_to_column("rowname") %>% 
      sample_n(size = num) %>% column_to_rownames("rowname") #Sample mutated (1) data
    
    bound_data <- rbind(data_oneChr_sample0, data_oneChr_sample1)
    
    return(bound_data)
  })
  
  data_evenTPTN <- as.data.frame(do.call(rbind, data_evenTPTN))
  
  
  
  # Create glm with leave-one-chromosome-out CV
  lapply(chroms, function(chrom){
    
    # Data for training, n-1 chromosomes 
    data_train <- subset(data_evenTPTN, chr != chrom, select = -c(chr))
    
    # glm
    logR <- glm(formula = mutated ~ ., data = data_train, 
                family = binomial(link = "logit"))
    save(logR, file = paste0(path, "/", chrom,"/logR_evenTPTN_testr_", chrom,".RData"))
  })
  
  
  # Test model on remaining chromosome 
  lapply(chroms, function(chrom){
    
    load(paste0(path, "/", chrom,"/logR_evenTPTN_testr_", chrom,".RData"))
    
    # Data for testing
    data_test <- subset(data_evenTPTN, chr == chrom, select = -c(chr, mutated))
    
    # Prediction
    yhat <- predict(logR, newdata = data_test, type = "response")
    save(yhat, file = paste0(path, "/", chrom,"/yhat_evenTPTN_testr_",chrom, ".RData"))
  })
  
  
  
  ##### Model performances #####
  # AUROC performances
  AUROC_CV <- lapply(chroms, function(chrom){ 

    load(paste0(path, "/", chrom,"/yhat_evenTPTN_testr_",chrom, ".RData"))
    
    data_test <- subset(data_evenTPTN,  chr == chrom)
    
    predict <- prediction(yhat, as.numeric(as.character(data_test$mutated)))
    
    return(performance(predict, "auc"))
  })
  names(AUROC_CV) <- chroms
  save(AUROC_CV, file = paste0(path, "/AUROC_evenTPTN_testr",i,"_CV.RData"))
  
})



# Get all calculated AUC for each run
chroms <- paste0("chr", 1:22)

allAUC <- sapply(1:5, function(i){
  load(paste0(path, "/AUROC_evenTPTN_testr",i,"_CV.RData"))
    
  sapply(chroms, function(chrom){
    print(chrom)
    return(AUROC_CV[[chrom]]@y.values[[1]])
 })
}) %>% as.data.frame()

  
### ANOVA and tukey test (Check, if the differences of the means for each chromosome are significant) ###
allAUC_gatherd <- gather(allAUC, key = "sampling", value = "AUROC", V1:V5)
allAUC_gatherd$chromosome <- rep(chroms, times = 5)

summary(aov(AUROC ~ chromosome, data=allAUC_gatherd)) # There are significant differences between the chromosomes 
AUC_tukey <- TukeyHSD(aov(AUROC ~ chromosome, data=allAUC_gatherd))

as.data.frame(subset(AUC_tukey$chromosome[,4], subset = AUC_tukey$chromosome[,4] <= 0.05))

write.csv(as.data.frame(subset(AUC_tukey$chromosome[,4], subset = AUC_tukey$chromosome[,4] <= 0.05))
          ,paste0(path, "/AUROC_tukeyp0.05_5samples.csv"), row.names = TRUE)  



### PLOTS ###
ggplot(data = allAUC_gatherd, aes(chromosome, AUROC)) +
  geom_hline(yintercept = mean(allAUC_gatherd$AUROC),colour='blue',size = 0.3) +
  geom_boxplot() +    
  labs(x = "Testing chromosome",y = "AUROC") +
  theme_classic() +
  theme(text = element_text(size=20))

ggsave(filename = paste0("/AUROC_evenTPTN_5samplings_",tissue,".png"), path = path_fig,
       width = 40, height = 20, dpi = 400, units = "cm")



##### All three in one plot #####
load(paste0(path, "/ROC_evenTPTN_CV.RData"))    # Load ROC
load(paste0(path, "/PR_evenTPTN_CV.RData"))     # Load PR

chroms = names(ROC_CV)
chrom_num = as.numeric(sort(as.character(c(1:22)))) # numbers for sorting

colors_chr <- as.vector(rainbow(22))

load(paste0(path, "/AUROC_evenTPTN_testr1_CV.RData"))     # Load AUROC
AUROC_full_evenTPTN = sapply(AUROC_CV, function(i){
  return(i@y.values[[1]])
})
load(paste0(path, "/AUROC_CV.RData"))     # Load AUROC
AUROC_CV = AUROC_CV[sort(names(AUROC_CV))]

AUROC_full_evenTPTN = as.data.frame(t(rbind(AUROC_full_evenTPTN, AUROC_CV)))
AUROC_full_evenTPTN = rownames_to_column(AUROC_full_evenTPTN)
AUROC_full_evenTPTN[,"Difference"] = AUROC_full_evenTPTN$AUROC_CV - AUROC_full_evenTPTN$AUROC_full_evenTPTN

colnames(AUROC_full_evenTPTN) = c("Chromosome", "AUROC_evenTPTN", "AUROC_fullData", "Difference")
AUROC_full_evenTPTN <- gather(data = AUROC_full_evenTPTN, key = Method, value = AUROC, "AUROC_evenTPTN": "Difference")
AUROC_full_evenTPTN = cbind(AUROC_full_evenTPTN, chrom_num)
AUROC_full_evenTPTN = AUROC_full_evenTPTN[order(AUROC_full_evenTPTN$chrom_num),]
AUROC_full_evenTPTN$Chromosome = as.factor(AUROC_full_evenTPTN$Chromosome)



png(filename = paste0(path_fig, "/ROC_PR_AUROCcomp_evenTPTN_CV_",tissue,".png"), width = 720*2, height = 720*2, pointsize = 22)
layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE))
par(cex.lab = 1.4, cex = 1.01)

# ROC
plot(ROC_CV$chr1, col = colors_chr[1])
sapply(seq(2,length(ROC_CV)), function(i){
  plot(ROC_CV[[i]], add = TRUE, col = colors_chr[i])
})
abline(coef = c(0,1), lty = 2)
legend("bottomright", legend = chroms,bg="transparent", col = colors_chr, ncol = 2, lwd = 5, lty=1, box.lty=0, seg.len = 0.8, 
       y.intersp = 1, x.intersp = 0.8, text.width = c(0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1))
#PR
plot(PR_CV$chr1, col = colors_chr[1], xlim = c(0,1), ylim = c(0,1))
sapply(seq(2,length(PR_CV)), function(i){
  plot(PR_CV[[i]], add = TRUE, col = colors_chr[i], xlim = c(0,1), ylim = c(0,1))
})

plot.new()              
vps <- baseViewports()
pushViewport(vps$figure) 
vp1 <-plotViewport(c(1.8,1,0,1))

# Compare AUROC before and after even TP and TN
g = ggplot(data = AUROC_full_evenTPTN, aes(x=reorder(Chromosome, chrom_num), y=AUROC, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge") + 
  theme_minimal() + 
  theme(legend.position = "bottom", legend.box = "vertical") + 
  theme(text = element_text(size=28))+
  scale_fill_manual(values=c("#56B4E9", "#E69F00", "#999999"),
                    name = "", labels = c("AUROC of unsampled data", "AUROC of data with even TP and TN", "Difference"))+
  labs(x = "Chromosome")
print(g, vp = vp1)
dev.off()



####### Select tissue #######
tissue = tissues[1]

path = paste0("data/rdata/", tissue, "/logRegression")
path_fig = paste0("fig/", tissue, "/logRegression")
#############################

chroms <- unique(data_evenTPTN$chr)

###### ROC ON PERMUTATED DATA #####
load(paste0(path, "/data_evenTPTN.RData"))

ROC_random <- lapply(chroms, function(chrom){ 
  
  load(paste0(path, "/", chrom,"/yhat_evenTPTN_",chrom,".RData"))
  
  data_test <- subset(data_evenTPTN,  chr == chrom)
  
  predict <- prediction(yhat, as.numeric(as.character(sample(x = data_test$mutated, size = length(yhat)))))
  
  return(performance(predict, "tpr", "fpr"))
})
names(ROC_random) <- chroms

allAUC_gatherd$chrom_num = rep(x = 1:22, 5)


# Plot random ROC
colors_chr <- as.vector(rainbow(22))

png(filename = paste0(path_fig, "/ROC_evenTPTN_5repBoxplot_random_CV_",tissue,".png"), width = 720*2, height = 720*2, pointsize = 22)
layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
par(cex.lab = 1.4, cex = 1.01)

# All AUROC
plot.new() 
vps <- baseViewports()
pushViewport(vps$figure) 
vp1 <-plotViewport(c(1.8,1,0,1))

p <- ggplot(data = allAUC_gatherd, aes(reorder(chromosome, chrom_num), AUROC)) +
  geom_hline(yintercept = mean(allAUC_gatherd$AUROC),colour='blue',size = 0.3) +
  geom_boxplot() +    
  labs(x = "Testing chromosome",y = "AUROC") +
  theme_classic() +
  theme(text = element_text(size=28))

print(p, vp = vp1)

plot(ROC_random$chr1, col = colors_chr[1])
sapply(seq(2,length(ROC_random)), function(i){
  plot(ROC_random[[i]], add = TRUE, col = colors_chr[i])
})
abline(coef = c(0,1), lty = 2)
legend("bottomright", legend = chroms,bg="transparent", col = colors_chr, ncol = 2, lwd = 5, lty=1, box.lty=0, seg.len = 0.8, 
       y.intersp = 0.9, x.intersp = 0.8, text.width = c(0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1))

  
dev.off()

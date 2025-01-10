setwd("/cellnet/MutationModel/")

### Packages ###
library(ggplot2)
library(gtools)
library(tidyr)
library(gridBase)
library(grid)
library(tibble)
###############


##### Data with bins #####
path = "data/rdata/"
path_fig = "/cellnet/MutationModel/fig/"
###########################

chroms = paste0("chr",1:22)

pCutOff <- 0.01


##### Select tissue #####
tissues = c("luad",  "skin", "colon", "ovary", "kidney", "prostate", "breast")

tissue_comp = c(tissues[1], tissues[2])
#########################

colors_t <- as.vector(rainbow(length(tissues)))


##### MODEL EVALUATION COMPARISON #####

##### LogLoss for all regularization methods ##### 
all_LogLoss_CV <- sapply(tissues, function(tissue){
  print(tissue)
  load(paste0(path, tissue,"/logRegression/logloss_CV_pCutOff",pCutOff,".RData"))
  return(logloss_CV)
}) %>% as.data.frame()

all_LogLoss_aov <- gather(all_LogLoss_CV, key = tissue, value = LogLoss, luad:breast)

summary(aov(all_LogLoss_aov$LogLoss ~ all_LogLoss_aov$tissue))
TukeyHSD(aov(all_LogLoss_aov$LogLoss ~ all_LogLoss_aov$tissue))

# Plot all logLoss for each method
png(paste0(path_fig, "boxplot_logLoss_pCutOff",pCutOff,"_alltissues.png"), width = 720*2, height = 720, pointsize = 22)
boxplot(all_LogLoss_CV, ylab = "LogLoss")
dev.off()



##### AUC VS absolute number of mutations per tissue 
all_AUC_concat <- sapply(tissues, function(tissue){
  print(tissue)
  load(paste0(path, tissue, "/logRegression/AUC_concat_pCutOff",pCutOff,".RData")) # AUC concat
  load(paste0(path,tissue,"/completeData_withBinwise.RData")) # Tissue dataset
  
  
  num_mut <- table(data$mutated)[["1"]] # Number of all mutations 
  num_all <-nrow(data)
  return(c(AUC_concat@y.values[[1]], num_mut, num_all))
}) %>% t() %>% as.data.frame()

colnames(all_AUC_concat) <- c("AUC", "Absolute number of mutations", "Absolute number of all observations")

# Plot
ggplot(data = all_AUC_concat, aes(`Absolute number of all observations`, AUC)) + 
  labs(x = "Absolute number of observations", y = "AUROC") + 
  geom_smooth(method = "lm", alpha = 0.3, se = FALSE) + 
  geom_point() +
  theme_classic() +
  theme(text = element_text(size=15), panel.grid.major = element_line(size = 0.3, linetype = 'solid', colour = "lightgrey"), 
        panel.grid.minor = element_line(size = 0.3, linetype = 'solid', colour = "lightgrey")) +
  geom_text(aes(label=tissues),hjust=-0.2, vjust=0)

ggsave(filename = "/datasizeVSAUROC_alltissues.png", path = paste0(path_fig),
       width = 30, height = 20, dpi = 300, units = "cm") 


###### ROC and PR #####
## ROC data
allROC_concat <- sapply(tissues, function(tissue){
  print(tissue)
  load(paste0(path, tissue,"/logRegression/ROC_concat_pCutOff",pCutOff,".RData"))
  return(ROC_concat)
})

## PR data
allPR_concat <- sapply(tissues, function(tissue){
  print(tissue)
  load(paste0(path, tissue,"/logRegression/PR_concat_pCutOff",pCutOff,".RData"))
  return(PR_concat)
})



# Plot
png(filename = paste0(path_fig, "ROC_PR_alltissues.png"), width = 720*2, height = 720, pointsize = 22)
par(mfrow = c(1,2), cex.lab = 1.4)

#ROC
plot(allROC_concat$luad, col = colors_t[1], lwd = 2)
sapply(seq(2,length(allROC_concat)), function(i){
  plot(allROC_concat[[i]], add = TRUE, col = colors_t[i], lwd = 2)
})
abline(coef = c(0,1), lty = 2)
legend(0.75, 0.45, legend = tissues, bg="transparent", col = colors_t, lwd = 5, lty=1, box.lty=0, seg.len = 0.8, 
       y.intersp = 1, x.intersp = 0.8, text.width = c(0.1,0.1,0.1,0.1,0.1,0.1,0.1))


#PR
plot(allPR_concat$luad, col = colors_t[1], xlim = c(0,1), ylim = c(0,1), lwd = 2.5)
sapply(seq(2,length(allPR_concat)), function(i){
  plot(allPR_concat[[i]], add = TRUE, col = colors_t[i], xlim = c(0,1), ylim = c(0,1), lwd = 2.5)
})
dev.off()



##### ROC, PR and AUC in one plot #####
png(filename = paste0(path_fig, "ROC_PR_AUC_alltissues.png"), width = 720*2, height = 720*2, pointsize = 22)
layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE))
par(cex.lab = 1.2, cex = 1.2)


#ROC
plot(allROC_concat$luad, col = colors_t[1], lwd = 2)
sapply(seq(2,length(allROC_concat)), function(i){
  plot(allROC_concat[[i]], add = TRUE, col = colors_t[i], lwd = 2)
})
abline(coef = c(0,1), lty = 2)
legend(0.75, 0.45, legend = tissues, bg="transparent", col = colors, lwd = 5, lty=1, box.lty=0, seg.len = 0.8, 
       y.intersp = 1, x.intersp = 0.8, text.width = c(0.1,0.1,0.1,0.1,0.1,0.1,0.1))


#PR
plot(allPR_concat$luad, col = colors_t[1], xlim = c(0,1), ylim = c(0,1), lwd = 2.5)
sapply(seq(2,length(allPR_concat)), function(i){
  plot(allPR_concat[[i]], add = TRUE, col = colors_t[i], xlim = c(0,1), ylim = c(0,1), lwd = 2.5)
})


# Mutation ratio per AUROC
plot.new()              
vps <- baseViewports()
pushViewport(vps$figure) 
vp1 <-plotViewport(c(1.8,1,0,1))

p <- ggplot(data = all_AUC_concat, aes(`Absolute number of all observations`, AUC)) + 
  labs(x = "Absolute number of observations", y = "AUROC") + 
  geom_smooth(method = "lm", alpha = 0.3, se = FALSE) + 
  geom_point() +
  theme_classic() +
  geom_text(aes(label=tissues),hjust=-0.2, vjust=0, size = 9)+
  theme(text = element_text(size=30))


print(p, vp = vp1)
dev.off()




##### COEFFICIENT AND FEATURE COMPARISSON #####
# List of all features with p <= 0.01 in the compared tissues 
all_tissues_p0.01 <- sapply(tissues, function(tissue){
  print(tissue)
  load(paste0(path, tissue,"/data_subset_pCutOff0.01.RData"))
  return(colnames(subset_data_features_padj)[-c(1,2)]) # All features w/o mutated and chr
})

df_all_tissues_p0.01 <- 
  for(l in 1:length(all_tissues_p0.01)){
    lapply(all_tissues_p0.01[l], function(x) write(x, paste0(path, tissues[l], "/feature_list_p0.01.txt"), append=TRUE, sep=";"))
  }



##### COMPARE COEFFICIENTS #####
## FILTERED Coefficients ##
all_Coefs_comp <- lapply(tissue_comp, function(tissue){ # Load all filtered coefficients for each tissue
  print(tissue)
  load(paste0(path, tissue,"/logRegression/coefficients_allChr_p",pCutOff,".RData"))
  return(coefficients_filteredPval)
})
names(all_Coefs_comp) <- tissue_comp


## Features ##
feature_comp <- sapply(rownames(all_Coefs_comp[[1]]), function(f){ # List of features that only appear in both tissues
  print(f)
  
  if(f %in% rownames(all_Coefs_comp[[2]])){
    return(f)
  }
})
feature_comp <- unique(unlist(feature_comp[!sapply(feature_comp,is.null)]))[-1]


all_Coefs_comp <- lapply(tissue_comp, function(i){ # Rearrange lists 
  all_Coefs_comp[[i]] <- cbind(predictors = as.factor(rownames(all_Coefs_comp[[i]])), all_Coefs_comp[[i]])
  rownames(all_Coefs_comp[[i]]) <- NULL
  all_Coefs_comp[[i]] <- gather(data = all_Coefs_comp[[i]], key = chromosome, value = coefficient, chr1:chr22)
})
names(all_Coefs_comp) <- tissue_comp



for(i in 1:length(tissue_comp)){ # Combine lists to one dataframe
  print(i)
  all_Coefs_comp[[i]]["tissue"] <- rep(tissue_comp[i], times = nrow(all_Coefs_comp[[i]]))
}
all_Coefs_comp <- rbind(all_Coefs_comp[[tissue_comp[1]]], all_Coefs_comp[[tissue_comp[2]]])


all_Coefs_comp_gg <- c()
for(feature in feature_luadskin){
  print(feature)
  for(row in 1:nrow(all_Coefs_comp))
    if(feature == all_Coefs_comp[row,]$predictors)
      all_Coefs_comp_gg <- rbind(all_Coefs_comp_gg, all_Coefs_comp[row,])
}
# Plot
# Jitter & boxplot
ggplot(all_Coefs_comp, aes(predictors, coefficient)) + 
  geom_boxplot() + 
  geom_jitter(shape = 20, aes(colour = tissue))+ 
  geom_hline(yintercept = 0) + 
  labs(x = "Feature",y = "Feature coefficient") +
  theme_classic()+
  theme(panel.grid.major = element_line(size = 0.3, linetype = 'solid', colour = "lightgrey"), 
        panel.grid.minor = element_line(size = 0.3, linetype = 'solid', colour = "lightgrey"))+
  theme(text = element_text(size=32), legend.position = c(0.15, 0.15), 
        legend.background = element_rect(fill="lightgrey",
                                         size=0.5, linetype="solid", 
                                         colour ="black"),
        axis.title=element_text(size=38), legend.text=element_text(size=34),
        legend.title=element_text(size=35))+
  coord_flip()+ 
  scale_color_manual(values=c("#E69F00","#56B4E9"), labels = c("lung", "skin"))+ 
  guides(colour = guide_legend(override.aes = list(size = 10)))



ggsave(filename = paste0("/all_coefs_",tissue_comp[1],"VS",tissue_comp[2],".pdf"), path = paste0(path_fig),
       width = 55, height = 68, dpi = 400, units = "cm")


# Scatterplot
all_Coefs_comp_gg <- spread(data = all_Coefs_comp_gg, key = tissue, value = coefficient)

ggplot(all_Coefs_comp_gg, aes(luad, skin)) + 
  geom_smooth(method = lm) +
  geom_point(shape = 20) +
  ggtitle(paste0("All significant variable coefficients of ", tissue_comp[1], " and ", tissue_comp[2]))+
  theme_classic()
ggsave(filename = paste0("/all_coefs_scatter_",tissue_comp[1],"VS",tissue_comp[2],".png"), path = paste0(path_fig),
       width = 40, height = 20, dpi = 300, units = "cm")


# Heatmap all filtered coefficients 
all_coeffs_mean <- lapply(tissues, function(tissue){
  print(tissue)
  load(paste0(path, tissue,"/logRegression/coefficients_allChr_p",pCutOff,".RData"))
  temp <- as.data.frame(rowMeans(coefficients_filteredPval))
  names(temp) <- tissue
  temp_names <- rownames_to_column(temp, var = "Features")
  
  return(temp_names)
})
names(all_coeffs_mean) = tissues


df_all_coeffs_mean <- all_coeffs_mean$luad
for(tissue in tissues[2:length(tissues)]){
  print(tissue)
  
  df_all_coeffs_mean <- merge(df_all_coeffs_mean, all_coeffs_mean[[tissue]], by='Features', all=TRUE)
}

df_all_coeffs_mean <- gather(df_all_coeffs_mean, key = "Tissue", value = "Coefficient", tissues)
df_all_coeffs_mean$Tissue <- factor(df_all_coeffs_mean$Tissue, levels = tissues)

# Plot
ggplot(df_all_coeffs_mean, aes(Tissue, Features, fill= Coefficient)) +
  geom_tile(color = "lightgrey") +
  scale_fill_gradient2(low = "red2", mid = "white", high = "blue", space = "Lab", 
                       name="coefficient") +
  theme_minimal()+
  theme(text = element_text(size=28))


ggsave(filename = paste0("/heatmap_FilteredCoeffs_alltissues.png"), path = paste0(path_fig),
       width = 45, height = 60, dpi = 400, units = "cm")



### ALL COEFFICIENTS (NOT FILTERED) ###
## Collect all coefficients
all_coeffs_mean <- lapply(tissues, function(tissue){
  print(tissue)
  load(paste0(path, tissue, "/logRegression/coefficients_allChr.RData"))
  temp <- as.data.frame(rowMeans(coefficients))
  names(temp) <- tissue
  temp_names <- rownames_to_column(temp, var = "Features")
  
  return(temp_names)
})
names(all_coeffs_mean) = tissues


df_all_coeffs_mean <- all_coeffs_mean$luad
for(tissue in tissues[2:length(tissues)]){
  print(tissue)
  
  df_all_coeffs_mean <- merge(df_all_coeffs_mean, all_coeffs_mean[[tissue]], by='Features', all=TRUE)
}

df_all_coeffs_mean <- gather(df_all_coeffs_mean, key = "Tissue", value = "Coefficient", tissues)
df_all_coeffs_mean$Tissue <- factor(df_all_coeffs_mean$Tissue, levels = tissues)


# Heatmap all coefficiens
ggplot(df_all_coeffs_mean, aes(Tissue, Features, fill= Coefficient)) +
  geom_tile(color = "lightgrey") + 
  theme_minimal()+
  scale_fill_gradient2(low = "red2", mid = "white", high = "blue", space = "Lab", 
                      name="coefficient") +
  theme(text = element_text(size=10)) 


ggsave(filename = paste0("/heatmap_allCoeffs_alltissues.png"), path = paste0(path_fig),
       width = 45, height = 60, dpi = 400, units = "cm")



## Collect all p-Values
all_pValFactor_mean <- lapply(tissues, function(tissue){
  print(tissue)
  load(paste0(path, tissue,"/logRegression/pValues_factors_allChr.RData"))
  temp <- as.data.frame(rowMeans(pValue_factors))
  names(temp) <- tissue
  temp_names <- rownames_to_column(temp, var = "Features")
  
  return(temp_names)
})
names(all_pValFactor_mean) = tissues


df_pValFactor_mean <- all_pValFactor_mean$luad
for(tissue in tissues[2:length(tissues)]){
  print(tissue)
  
  df_pValFactor_mean <- merge(df_pValFactor_mean, all_pValFactor_mean[[tissue]], by='Features', all=TRUE)
}


df_pValFactor_mean <- gather(df_pValFactor_mean, key = "Tissue", value = "p-Value", tissues)
df_pValFactor_mean$Tissue <- factor(df_pValFactor_mean$Tissue, levels = tissues)

df_all_coeffs_pVal_mean <- cbind(df_pValFactor_mean, coefficients = df_all_coeffs_mean$Coefficient)

# Vulcano plot
ggplot(df_all_coeffs_pVal_mean, aes(x=coefficients, y=-log10(`p-Value`), color = Tissue)) + 
  geom_point()+ 
  labs(y = "-log10(p-value)")+
  theme_minimal()+ 
  geom_text(aes(label = ifelse(coefficients >=0.2 | coefficients <= -0.2 | -log10(`p-Value`) >= 45, yes = Features, no = NA)),
            nudge_y = 4.5, size = 3)

ggsave(filename = paste0("/vulcano_alltissues.png"), path = paste0(path_fig),
       width = 20, height = 20, dpi = 400, units = "cm")



# Heatmap all coefficients with all significant marked 
ggplot(df_all_coeffs_pVal_mean, aes(Tissue, Features, fill= coefficients)) +
  geom_tile(color = "lightgrey") + 
  theme_minimal()+
  scale_fill_gradient2(low = "red2", mid = "white", high = "blue", space = "Lab", 
                       name="coefficient") +
  theme(text = element_text(size=28))+ 
  geom_text(aes(label = ifelse(`p-Value`<=0.01, yes = "*", no = NA)), color = "grey45", size = 17, nudge_y = -0.45)

ggsave(filename = paste0("/heatmap_allCoeffs_pMarked_alltissues.png"), path = paste0(path_fig),
       width = 45, height = 60, dpi = 400, units = "cm")



##### P-VALUE COMPARISON #####
# Heatmap 
all_pVals_mean <- lapply(tissues, function(t){
  print(t)
  load(file = paste0(path, t,"/logRegression/pValues_allChr_adjust.RData"))
  temp <- as.data.frame(rowMeans(pValues_adjust))
  names(temp) <- t
  temp_r <- rownames_to_column(temp, var = "Features")
  
  return(temp_r)
})
names(all_pVals_mean) = tissues


df_all_pVals_mean <- all_pVals_mean$luad
for(tissue in tissues[2:length(tissues)]){
  print(tissue)
  
  df_all_pVals_mean <- merge(df_all_pVals_mean, all_pVals_mean[[tissue]], by='Features', all=TRUE)
}


df_all_pVals_mean <- gather(df_all_pVals_mean, key = "Tissue", value = "p-Value", tissues)
df_all_pVals_mean$Tissue <- factor(df_all_pVals_mean$Tissue, levels = tissues)

# Plot
ggplot(df_all_pVals_mean, aes(Tissue, Features, fill= `p-Value`)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "red2", high = "white", limit = c(0,1), space = "Lab", 
                      name="p-Value") +
  theme_minimal()+
  theme(text = element_text(size=30))+ 
  geom_text(aes(label = ifelse(`p-Value`<=0.01, yes = "*", no = NA)), color = "grey", size = 17, nudge_y = -0.45)


ggsave(filename = paste0("/heatmap_pVals_pMarked_alltissues.png"), path = paste0(path_fig),
       width = 45, height = 60, dpi = 400, units = "cm")

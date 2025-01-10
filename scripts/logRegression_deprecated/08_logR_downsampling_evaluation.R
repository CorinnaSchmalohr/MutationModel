setwd("/cellnet/MutationModel/")

### Packages ###
library(ROCR)
library(ggplot2)
###############

tissues = c("luad",  "skin", "colon", "ovary", "kidney", "prostate", "breast")


##### EVALUATION OF THE MODEL ON DOWNSAMPLED DATA #####
### CV-wise performance of the  model trained on downsampled data ###
lapply(tissues, function(tissue){  
  print(tissue)
  
  ##### Paths & Data #####
  path = paste0("data/rdata/", tissue, "/logRegression")
  
  path_fig = paste0("fig/", tissue, "/logRegression")
  
  
  # Data
  chroms = paste0("chr", seq(1,22)) # List of all chromosomes
  colors_chr <- as.vector(rainbow(22))
  
  load(paste0(path, "/ROC_downsampling_rep1_CV.RData"))       # Load ROC
  load(paste0(path, "/PR_downsampling_rep1_CV.RData"))        # Load PR
  
  
  ##### Plots #####
  # ROC and PR 
  png(filename = paste0(path_fig, "/ROC_PR_downsampling_CV_",tissue,".png"), width = 720*2, height = 720, pointsize = 22)
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


### ALL TISSUE COMPARISON  ###
path = "data/rdata/"
path_fig = "/cellnet/MutationModel/fig/"

# ROC AND PR
allROC_concat <- sapply(tissues, function(tissue){
  print(tissue)
  load(paste0(path, tissue,"/logRegression/ROC_downsampling_rep1_concat.RData"))
  return(ROC_concat)
})


allPR_concat <- sapply(tissues, function(tissue){
  print(tissue)
  load(paste0(path, tissue,"/logRegression/PR_downsampling_rep1_concat.RData"))
  return(PR_concat)
})

# Plot
colors_t <- as.vector(rainbow(length(tissues)))

png(filename = paste0(path_fig, "ROC_PR_downsampling_alltissues.png"), width = 720*2, height = 720, pointsize = 22)
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



### COMPARE RESULTS BEFORE AND AFTER DOWNSAMPLING ###
# Compare ROC
allROC_downsampling_concat_allTissues<- sapply(tissues, function(t){
  path = paste0("data/rdata/", t, "/logRegression")
  load(file = paste0(path, "/ROC_downsampling_rep1_concat.RData"))
  return(ROC_concat) # Return different model on data performance
})

allROC_concat_allTissues <- sapply(tissues, function(t){
  path = paste0("data/rdata/", t, "/logRegression")
  load(file = paste0(path, "/ROC_concat_pCutOff0.01.RData"))
  return(ROC_concat) # Return ROC for model and data from the same tissue
})


colors_t <- as.vector(rainbow(length(tissues)))
names(colors_t) <- tissues


png(filename = paste0("/cellnet/MutationModel/fig/ROC_downsampling_allTissues.png"), width = 720, height = 720, pointsize = 22)

#ROC
plot(allROC_downsampling_concat_allTissues[[1]], col = colors_t[names(allROC_downsampling_concat_allTissues[1])], lwd = 2)
sapply(seq(2,length(allROC_downsampling_concat_allTissues)), function(i){
  plot(allROC_downsampling_concat_allTissues[[i]], add = TRUE, col = colors_t[names(allROC_downsampling_concat_allTissues[i])], lwd = 2)
})
sapply(seq(1,length(allROC_concat_allTissues)), function(i){
  plot(allROC_concat_allTissues[[i]], add = TRUE, col = colors_t[names(allROC_concat_allTissues[i])], lwd = 2, lty = 3)
})
abline(coef = c(0,1), lty = 2)
legend(0.7, 0.5, legend = tissues, bg="transparent", col = colors_t, lwd = 5, lty=1, box.lty=0, seg.len = 0.8, 
       y.intersp = 1, x.intersp = 0.8, text.width = c(0.1,0.1,0.1,0.1,0.1,0.1, 0.1))
legend(0.15, 0.18, legend = c("downsampled model", "full model"), bg="transparent", col = "darkgrey", lwd = 2, lty=c(1,3), box.lty=0, seg.len = 0.8, 
       y.intersp = 1, x.intersp = 0.8, text.width = c(0.1,0.1))

dev.off()  



# Compare AUROC
allAUC_concat_allTissues_comp <- sapply(tissues, function(t){
  path = paste0("data/rdata/", t, "/logRegression")
  load(file = paste0(path, "/AUC_downsampling_concat.RData"))
  AUROC_d = AUC_concat@y.values[[1]]
  
  load(file = paste0(path, "/AUC_downsampling_rep1_concat.RData"))
  AUROC_rep = AUC_concat@y.values[[1]]
  
  load(file = paste0(path, "/AUC_concat_pCutOff0.01.RData"))
  
  return(c(AUROC_downsampling = AUROC_d, AUROC_downsampling_rep = AUROC_rep,AUROC = AUC_concat@y.values[[1]])) # Return different model on data performance
}) %>% as.data.frame()

allAUC_concat_allTissues_comp <- rownames_to_column(allAUC_concat_allTissues_comp, var = "data")
allAUC_concat_allTissues_comp <- gather(data = allAUC_concat_allTissues_comp, key = Tissues, value = AUROC, "luad":"breast")

ggplot(allAUC_concat_allTissues_comp[allAUC_concat_allTissues_comp$data != "AUROC_downsampling_rep",], aes(x = Tissues, y = AUROC, fill=data)) + 
  geom_bar(position="dodge", stat="identity") + 
  theme_minimal()+
  theme(legend.position = "bottom", legend.box = "vertical")+
  scale_fill_manual(name = "", labels = c("AUROC full sized data", "AUROC reduced data"), 
                    values=c("#56B4E9", "#E69F00"))+
  theme(text = element_text(size=20))

ggsave(filename = paste0("/AUROC_downsampling_rep1_allTissues.png"), path = path_fig,
       width = 30, height = 20, dpi = 300, units = "cm") 




###### COMAPRE ALL DOWNSAMPLED MODEL CROSS-TISSUE-APPLICATIONS AT ONCE -> HEATMAPs #####
tissue_combs <- expand.grid(tissues, tissues)

all_AUROC = sapply(1:nrow(tissue_combs), function(t){
  tissue_model = as.character(tissue_combs[t,1])
  path_model = paste0("data/rdata/", tissue_model, "/logRegression")
  
  tissue_apply = as.character(tissue_combs[t,2])
  
  print(c(tissue_model, tissue_apply))
  
  if(tissue_model == tissue_apply){ # Load same-tissue AUROCs for all four test types
  print("Application on its own tissue ...")
      
    load(paste0(path_model, "/AUC_concat_pCutOff0.01.RData"))
    AUROC_modelFull_testFull = AUC_concat@y.values[[1]]
    
    load(paste0(path_model, "/AUC_downsampling_rep1_concat.RData"))
    AUROC_modelDown_testDown = AUC_concat@y.values[[1]]
    
    load(paste0(path_model, "/AUC_downsampling_modelDown_testFull_concat.RData"))
    AUROC_modelDown_testFull = AUC_concat@y.values[[1]]
    
    load(paste0(path_model, "/AUC_downsampling_modelFull_testDown_concat.RData"))
    AUROC_modelFull_testDown = AUC_concat@y.values[[1]]
    
    
  } else { #Load all cross-tissue AUROCs for all four test types 
    print("Cross-tissue appplication ...")
    
    load(paste0(path_model, "/AUC_concat_modelON",tissue_apply,".RData"))
    AUROC_modelFull_testFull = AUC_concat@y.values[[1]]
    
    load(paste0(path_model, "/AUC_downsampling_rep1_concat_modelON",tissue_apply,".RData"))
    AUROC_modelDown_testDown = AUC_concat@y.values[[1]]
    
    load(paste0(path_model, "/AUC_downsampling_modelDown_testFull_concat_modelON",tissue_apply,".RData"))
    AUROC_modelDown_testFull = AUC_concat@y.values[[1]]
    
    load(paste0(path_model, "/AUC_downsampling_modelFull_testDown_concat_modelON",tissue_apply,".RData"))
    AUROC_modelFull_testDown = AUC_concat@y.values[[1]]
    
  }
  
  return(c(tissue_model, tissue_apply, AUROC_modelFull_testFull, AUROC_modelDown_testDown, AUROC_modelDown_testFull, AUROC_modelFull_testDown))
}) %>% t() %>% as.data.frame()


colnames(all_AUROC) <- c("Model", "Tissue","AUROC_modelFull_testFull", "AUROC_modelDown_testDown", 
                         "AUROC_modelDown_testFull", "AUROC_modelFull_testDown")


all_AUROC = transform(all_AUROC, 
                      AUROC_modelFull_testFull = as.numeric(AUROC_modelFull_testFull), 
                      AUROC_modelDown_testDown = as.numeric(AUROC_modelDown_testDown), 
                      AUROC_modelDown_testFull = as.numeric(AUROC_modelDown_testFull), 
                      AUROC_modelFull_testDown = as.numeric(AUROC_modelFull_testDown))


# Manually plot all four heatmaps
ggplot(all_AUROC, aes(Model, Tissue, fill = AUROC_modelDown_testFull)) +
  geom_tile() +
  theme_minimal()+
  scale_fill_viridis_c(limits =c(0.40, 0.65), name = "AUROC") +
  ggtitle("Downsampled model tested on full data")

ggsave(filename = paste0("heatmap_AUROC_modelDown_testFull_allModelsONallTissues.png"), path = "/cellnet/MutationModel/fig", 
       width = 12, height = 12, units = "cm")


# Heatmaps next to each other 
require(gridExtra)

full_full = ggplot(all_AUROC, aes(Model, Tissue, fill = AUROC_modelFull_testFull)) +
  geom_tile() +
  theme_minimal()+
  scale_fill_viridis_c(limits =c(0.40, 0.65), name = "AUROC") +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(), 
        legend.position = "none", 
        text = element_text(size=18))


down_down = ggplot(all_AUROC, aes(Model, Tissue, fill = AUROC_modelDown_testDown)) +
  geom_tile() +
  theme_minimal()+
  scale_fill_viridis_c(limits =c(0.40, 0.65), name = "AUROC") +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(), 
        legend.position = "none", 
        text = element_text(size=18))

grid.arrange(full_full, down_down, ncol=2)
ggsave(filename = paste0("heatmap_full_down_allModelsONallTissues.pdf"),plot =  grid.arrange(full_full, down_down, ncol=2), 
       path = path_fig, width = 30, height = 20, dpi = 300, units = "cm")  

# library(ggplot2)
library(ranger)

tissues = c("luad",  "skin", "colon", "ovary",
            "kidney", "prostate", "breast")
for (tissue in tissues){
   print(tissue)
   filenames = c(load(paste0("data/rdata/", tissue, "/Muts.RData")),
                 load(paste0("data/rdata/", tissue, "/mers.RData")),
                 if(tissue != "kidney"){
                    load(paste0("data/rdata/", tissue, "/GTEx_eqtl.RData"))
                 },
                 load(paste0("data/rdata/", tissue, "/structures.RData")),
                 load(paste0("data/rdata/", tissue, "/structures_100bp.RData")),
                 load(paste0("data/rdata/", tissue, "/healthyExpr.RData")),
                 load(paste0("data/rdata/", tissue, "/cancerExpr.RData")),
                 load(paste0("data/rdata/", tissue, "/DNAaccessibility.RData"))[3],
                 load(paste0("data/rdata/", tissue, "/methylation.RData")),
                 if(tissue != "ovary"){
                    load(paste0("data/rdata/", tissue, "/methbank.RData"))
                 },
                 load(paste0("data/rdata/", tissue, "/GCcontent.RData")),
                 load(paste0("data/rdata/", tissue, "/DNAbinding.RData")),
                 load(paste0("data/rdata/", tissue, "/DNAbinding_100bp.RData")),
                 load(paste0("data/rdata/", tissue, "/UCSC_histones.RData")),
                 load(paste0("data/rdata/", tissue, "/UCSC_histones_100bp.RData")),
                 load(paste0("data/rdata/", tissue, "/Tfbs.RData")),
                 load(paste0("data/rdata/", tissue, "/Tfbs_100bp.RData")),
                 load(paste0("data/rdata/", tissue, "/conservation.RData")),
                 load(paste0("data/rdata/", tissue, "/conservation_100bp.RData")),
                 load(paste0("data/rdata/", tissue, "/replication.RData")),
                 load(paste0("data/rdata/", tissue, "/mappability_repeats.RData")),
                 load(paste0("data/rdata/", tissue, "/inexon.RData")),
                 load(paste0("data/rdata/", tissue, "/Nucleosome.RData")),
                 load(paste0("data/rdata/", tissue, "/Nucleosome_100bp.RData")),
                 if(tissue %in% c("luad", "ovary")){
                    load(paste0("data/rdata/", tissue, "/HiC.RData"))
                 } else if (tissue == "prostate"){
                    load(paste0("data/rdata/", tissue, "/HiC_interactions.RData"))
                 } else if (tissue %in% c("breast", "kidney", "colon", "skin")){
                    load(paste0("data/rdata/", tissue, "/HiC.RData"))
                 },
                 load(paste0("data/rdata/", tissue,
                             "/binwisePreds_allChrs.RData")))
   # if(tissue != "ovary"){meth = meth[,"methbank_100bp"]}
   # GCcontent = GCcontent[,"GCcontent_100bp"]
   data = data.frame(lapply(filenames, function(f){
      x = get(f)
      if(is.null(dim(x))){
         x = data.frame(x)
         colnames(x) = f
      }
      return(x)
   })) #, check.names=F
   
   data$context = NULL
   data$pentamer = NULL
   data$trimer = NULL
   data$septamer = NULL
   data$inexon = NULL
   data$replDirection = as.factor(data$replDirection)
   data$ref = factor(data$ref)
   data$mutated = as.factor(data$mutated)
   data$strand  = as.integer(data$strand == "+")
   miss = apply(data[,-(which(colnames(data) == "alt"))],
                1, function(x){sum(is.na(x))})
   data = data[miss == 0,]
   toRM = c("chr", "pos", "alt", "geneID", "baseGeneID")
   removed = data[,toRM]
   data = data[,(!colnames(data) %in% toRM)]
   
   toScale = sapply(data,is.numeric) & !sapply(data,is.integer)
   data[,toScale] = scale(data[,toScale])
   
   # random Forest #####
   
   dataSub = data[removed$chr == "chr1",]
   
   rf = ranger(mutated ~ ., data = dataSub, importance = 'permutation',
               write.forest = T, seed = 1234, num.threads =  10, num.trees=2000,
               respect.unordered.factors = 'partition', 
               scale.permutation.importance = T, probability = T)
   save(rf, file = paste0("data/rdata/", tissue,
                          "/RFmodel/testOnChr1_withBinWise.RData"))
   # rf$prediction.error
   imp = rf$variable.importance[order(names(rf$variable.importance))]
   resolutions = c("5bp", "10bp", "1kb","100bp","bins10kb",
                   "bins100kb", "bins1Mb")
   labels = t(sapply(names(imp), function(x){
      temp = strsplit(x,split = "_", fixed = T)[[1]]
      m = length(temp)
      if(temp[m] %in% resolutions){
         return(c(paste(temp[-m],collapse="_"), temp[m]))
      } else {
         return(c(x,"1bp"))
      } 
   }))
   dat = data.frame(names(imp), imp, labels, stringsAsFactors=F)
   dat$col  = dat$X1 %in% unique(dat$X1)[1:(length(unique(dat$X1))/2)]

   pdf(paste0("fig/", tissue, 
              "/RF_variable_importance_chr1Test.pdf"),
       height=20, width=8)
   # plot(ggplot(dat, aes(x=X1, y=imp, fill = X2))+
   ggplot(dat, aes(x=reorder(X1, imp), y=imp, fill = X2))+      
           # reorder(day, -perc), y = perc
      geom_bar(stat="identity",
               position = position_dodge2(preserve = "single"))+
      coord_flip() +
      theme_minimal() +
      theme(axis.text=element_text(hjust = 0.5, size=15))+
      ylab("RF permutation importance") +
      xlab("") +
      labs(fill = "res.") +
      facet_wrap(~col, scales="free_y")
   dev.off()
   
   subDat = dat[dat$X1 %in% dat$X1[duplicated(dat$X1)],]
   subDat$col  = subDat$X1 %in% unique(subDat$X1)[1:(length(unique(subDat$X1))/2)]
   
   pdf(paste0("fig/", tissue,
              "/RF_variable_importance_chr1Test_onlyMultiRes.pdf"), 
       height=20, width=8)
   # plot(ggplot(subDat, aes(x=X1, y=imp, fill = X2))+
   ggplot(subDat, aes(x=reorder(X1, imp), y=imp, fill = X2))+      
      geom_bar(stat="identity",
               position = position_dodge2(preserve = "single"))+
      coord_flip() +
      theme_minimal() +
      theme(axis.text=element_text(hjust = 0.5, size=15))+
      ylab("RF permutation importance") +
      xlab("") +
      labs(fill = "res.")+
      facet_wrap(~col, scales="free_y")
   dev.off()
}

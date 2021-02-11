# load general data #####
library(corrplot)
library(reshape2)
library(ggplot2)
library(png)
library(gridExtra)
library(grid)

sizes = c("bins1Mb" = 1000000, "bins100kb" = 100000,
          "bins10kb" = 10000, "bins1kb" = 1000, "bins100bp" = 100) #  "bins10bp" = 10
tissues = c("luad","skin", "kidney","colon",
            "ovary",  "prostate" ,"breast")
load(paste0("data/procData/bins/bins1Mb/bins1Mb.RData"))
load("data/rdata/binsAnalysis_features.RData")

MutDensitiesWGS = lapply(names(sizes), function(x){
   load(paste0("data/rdata/MutsPerBin_WGS_", x, ".RData"))
   return(WGSMutsPerBins)
})
names(MutDensitiesWGS) = names(sizes)
GCcontents = lapply(names(sizes), function(x){
   load(paste0("data/procData/bins/",x,"/GCcontent.RData"))
   return(GCcontent)
})
names(GCcontents) = names(sizes)
MutDensities = lapply(names(sizes), function(x){
   load(paste0("data/rdata/MutsPerBin_", x, ".RData"))
   return(MutsPerBin)
})
names(MutDensities) = names(sizes)

#####

# create Tables with all binned features for each tissue and resolution  #####
wholeBin = sapply(tissues, function(tissue){
   print(tissue)
   res = lapply(names(sizes), function(size){
      cat(size, ' ')
      load(paste0("data/procData/bins/", size,"/", tissue,
                  "/DNAbinding/DNAbinding.RData"))
      load(paste0("data/procData/bins/", size, "/", tissue,
                  "/UCSC_tracks/UCSC_histones.RData"))
      load(paste0("data/procData/bins/", size,"/",
                  tissue, "/UCSC_tracks/Tfbs.RData"))
      load(paste0("data/procData/bins/", size, "/", tissue,
                  "/DNAaccessibility/DNAaccessibility.RData"))
      colnames(DNAaccessibility) = c("DNAaccessibility",
                                     "UCSC_DNAaccessibility")
      load(paste0("data/procData/bins/", size, "/", tissue,
                  "/Nucleosome.RData"))
      GCcontent = GCcontents[[size]]
      MutsPerBin = MutDensities[[size]]
      if(tissue != "ovary"){
         load(paste0("data/procData/bins/", size, "/", tissue,
                     "/methylation/methylation.RData"))
      } else{
         methylation = rep(NA, nrow(MutsPerBin))
      }
      if(tissue %in%  c("skin","ovary", "kidney", "prostate")){
         WGSMutsPerBins = MutDensitiesWGS[[size]][,tissue]
      } else{
         WGSMutsPerBins = rep(NA, nrow(MutsPerBin))
      }
      dat = cbind(MutsPerBin[,c("chr", "start", "end", tissue)],
                  WGSMutsPerBins,
                  DNAbinding, histones, TF_BS, GCcontent, 
                  DNAaccessibility, methylation, Nucleosome) 
      colnames(dat)[4:5] = c("nMuts_WXS", "nMuts_WGS")
      save(dat, file = paste0("data/rdata/", tissue,
                              "/binsDat_", size, ".RData"))
      return(colnames(dat))
   })
   cat('\n')
   return(res)
})
covered = sapply(tissues, function(tissue){
   print(tissue)
   res = sapply(names(sizes), function(size){
      cat(size, ' ')
      load(paste0("data/procData/bins/", size,"/", tissue,
                  "/DNAbinding/covered_DNAbinding.RData"))
      load(paste0("data/procData/bins/", size, "/", tissue,
                  "/UCSC_tracks/covered_UCSC_histones.RData"))
      load(paste0("data/procData/bins/", size,"/", tissue,
                  "/UCSC_tracks/covered_Tfbs.RData"))
      load(paste0("data/procData/bins/", size, "/", tissue,
                  "/Nucleosome_covered.RData"))
      colnames(TF_BS_covered) = c("ETS_BS", "TF_BS")
      load(paste0("data/procData/bins/", size, "/", tissue,
                  "/DNAaccessibility/DNAaccessibility_covered.RData"))
      colnames(DNAaccessibility_covered) = c("DNAaccessibility",
                                             "UCSC_DNAaccessibility")
      load(paste0("data/procData/bins/", size, "/", tissue,
                  "_covered_GCcontent.RData"))
      MutsPerBin = MutDensities[[size]]
      if(tissue != "ovary"){
         load(paste0("data/procData/bins/", size, "/", tissue,
                     "/methylation/methylation_covered.RData"))
      } else{
         methylation_covered = rep(NA, nrow(MutsPerBin))
      }
      dat = cbind(MutsPerBin[,c("chr", "start", "end", tissue)],
                  DNAbinding_covered, histones_covered, 
                  TF_BS_covered, GCcontent = GCcontent_covered,
                  DNAaccessibility_covered,
                  methylation=methylation_covered, 
                  Nucleosome_covered)
      colnames(dat)[4] = "nMuts"
      save(dat, file = paste0("data/rdata/", tissue,
                              "/binsDat_", size, "_covered.RData"))
      return(colnames(dat))
   })
   cat('\n')
   return(res)
})
features = unique(unlist(wholeBin))[-(1:5)]
save(features, file = "data/rdata/binsAnalysis_features.RData")
#####

# general plot, only covered portions #####
print("only covered portions")
covered = sapply(tissues, function(tissue){
   print(tissue)
   png(paste0("fig/", tissue, "/genomePlots_", tissue, "_covered.png"),
       height = 2500, width = 4000, pointsize=20)
   split.screen(rbind(c(0,1,0.2,1),
                      c(0,1,0,0.2)))
   row1 = split.screen(screen=1,
                       figs=cbind(seq(0,0.8,0.2),
                                  seq(0.2,1,0.2),0,1))
   row2 = split.screen(screen=2,
                       figs=cbind(seq(0,0.8,0.2), 
                                  seq(0.2,1,0.2),0,1))
   assignment = rbind(row1, row2)
   colnames(assignment) = names(sizes)
   res = sapply(names(sizes), function(size){
      load(paste0("data/rdata/", tissue,
                  "/binsDat_", size, "_covered.RData"))
      dat = dat[!is.na(dat$nMuts),]
      cols = c("black", rainbow(ncol(dat)-4,end=0.9))
      names(cols) =  colnames(dat[-(1:3)])
      
      # genome plot
      p = assignment[1,size]
      assignment2 = split.screen(screen=p,
                                 c(ncol(dat)-3,1))
      names(assignment2) = c(colnames(dat)[-(1:3)])
      xlim = c(0,tail(bins$end,1))
      xlabels = seq(0,250000000,by=25000000)
      par(mar = c(0.2,2.5,0.1,0.1), oma = c(3,0,1,0))
      for(i in colnames(dat)[-(1:3)]){
         p2 = assignment2[i]
         screen(p2)
         par(mar = c(0.1,2.5,0.1,0.1))
         if(any(!is.na(dat[,i]))){
            plot(dat$end, dat[,i], type = "h", xlim = xlim,
                 col = rgb(t(col2rgb(cols[i])),alpha=250, maxColorValue=255),
                 xlab = "Mb", ylab = "", xaxt = "n", las = 1,tcl = -0.2,mgp = c(3,0.25,0),
                 yaxs = "i", lab = c(10,2,0), pch = 19, lwd = 1)
         } else{
            plot(dat$end, rep(0, nrow(dat)), type = "h", xlim = xlim,
                 col = rgb(t(col2rgb(cols[i])),alpha=250, maxColorValue=255),
                 xlab = "Mb", ylab = "", xaxt = "n", las = 1,tcl = -0.2,mgp = c(3,0.25,0),
                 yaxs = "i", lab = c(10,2,0), pch = 19, lwd = 1)
         }
         abline(v = xlabels,
                lty = 2, col = "grey")
         legend("topleft", legend= NA, bty = "n", title = i)
         if(i == colnames(dat)[4])
            mtext(size, side = 3, line = 0, cex=1.2)
      }
      axis(1, at=xlabels,labels=xlabels/1000000, line = 0,tcl = -0.2,
           mgp = c(3,0.2,0))
      mtext("genomic position (Mb)", side=1, line=0.7)
      # corrplot
      p = assignment[2,size]
      par(oma = c(0,0,0,0), mar = c(0,0,0,0))
      screen(p)
      par(oma = c(0,0,0,0), mar = c(0,1,0,1))
      mat = cor(dat[,-(1:3)], use = "pair")[(ncol(dat)-3):1,]
      image(mat, xaxt = "n",yaxt = "n", col=cm.colors(50))
      axis(2,at = seq(0,1,length.out=ncol(dat)-3),
           labels=rownames(mat), las = 2, tick=F, line=-1.3, hadj=0)
      # correlations
      cors = cor(dat[,4], dat[,-(1:4)], use="pair")[1,]
      return(abs(cors))
   })
   close.screen(all.screens=TRUE)
   dev.off()
   temp = melt(res)
   temp$Var2 = factor(temp$Var2, levels=c("bins100bp", "bins1kb", 
                                          "bins10kb", "bins100kb", "bins1Mb"))
   p = ggplot(temp, aes(x = Var2, y = Var1)) + 
      geom_raster(aes(fill=value)) + 
      scale_fill_gradient(low="grey90", high="red") +
      ggtitle(tissue) +
      xlab("bin size") + 
      ylab("features") + 
      labs(fill="abs. Pearson's R") +
      theme(plot.title = element_text(hjust = 0.5, size=25),
            axis.title=element_text(size=20),
            axis.text.x=element_text(angle=30, hjust = 1),
            axis.text=element_text(size=18), 
            legend.text=element_text(size=18),
            legend.title=element_text(size=18))
   ggsave(p, filename = paste0(
      "fig/", tissue, "/compareCorrelationPerBinsize_", tissue, "_covered.png"), 
      width=10, height=10)
})
#####


# general plot per whole bin #####
print("whole bins")
wholeBin = sapply(tissues, function(tissue){
   print(tissue)
   png(paste0("fig/", tissue, "/genomePlots_", tissue, ".png"),
       height = 2500, width = 4000, pointsize=20)
   split.screen(rbind(c(0,1,0.2,1),
                      c(0,1,0,0.2)))
   row1 = split.screen(screen=1,
                       figs=cbind(seq(0,0.8,0.2),
                                  seq(0.2,1,0.2),0,1))
   row2 = split.screen(screen=2,
                       figs=cbind(seq(0,0.8,0.2), 
                                  seq(0.2,1,0.2),0,1))
   assignment = rbind(row1, row2)
   colnames(assignment) = names(sizes)
   res = lapply(names(sizes), function(size){
      load(paste0("data/rdata/", tissue,
                  "/binsDat_", size, ".RData"))
      cols = c("black", "grey", rainbow(ncol(dat)-5))
      names(cols) = colnames(dat[-(1:3)])
      
      # genome plot
      p = assignment[1,size]
      assignment2 = split.screen(screen=p,
                                 c(ncol(dat)-3,1))
      names(assignment2) = c(colnames(dat)[-(1:3)])
      xlim = c(0,tail(bins$end,1))
      xlabels = seq(0,250000000,by=25000000)
      par(mar = c(0.2,2.5,0.1,0.1), oma = c(3,0,1,0))
      for(i in colnames(dat)[-(1:3)]){
         p2 = assignment2[i]
         screen(p2)
         par(mar = c(0.2,2.5,0.1,0.1))
         if(any(!is.na(dat[,i]))){
            plot(dat$end, dat[,i], type = "h", xlim = xlim,
                 col = rgb(t(col2rgb(cols[i])),alpha=250, maxColorValue=255),
                 xlab = "Mb", ylab = "", xaxt = "n", las = 1,tcl = -0.2,mgp = c(3,0.25,0),
                 yaxs = "i", lab = c(10,2,0), pch = 19, lwd = 1)
         } else{
            plot(dat$end, rep(0, nrow(dat)), type = "h", xlim = xlim,
                 col = rgb(t(col2rgb(cols[i])),alpha=250, maxColorValue=255),
                 xlab = "Mb", ylab = "", xaxt = "n", las = 1,tcl = -0.2,mgp = c(3,0.25,0),
                 yaxs = "i", lab = c(10,2,0), pch = 19, lwd = 1)
         }
         abline(v = xlabels,
                lty = 2, col = "grey")
         legend("topleft", legend= NA, bty = "n", title = i)
         if(i == colnames(dat)[4])
            mtext(size, side = 3, line = 0)
      }
      axis(1, at=xlabels,labels=xlabels/1000000, line = 0,tcl = -0.2,
           mgp = c(3,0.2,0))
      mtext("genomic position (Mb)", side=1, line=0.5)
      # corrplot
      p = assignment[2,size]
      par(oma = c(0,0,0,0), mar = c(0,0,0,0))
      screen(p)
      par(oma = c(0,0,0,0), mar = c(0,1,0,1))
      mat = cor(dat[,-(1:3)], use = "pair")[(ncol(dat)-3):1,]
      image(mat, xaxt = "n",yaxt = "n", col=cm.colors(50))
      axis(2,at = seq(0,1,length.out=ncol(dat)-3),
           labels=rownames(mat), las = 2, tick=F, line=-1.3, hadj=0)
      # correlations
      cors = cor(dat[,4:5], dat[,-(1:5)], use="pair")
      return(abs(cors))
   })
   close.screen(all.screens=TRUE)
   dev.off()
   
   WGSres = sapply(res, function(x){x["nMuts_WGS",]})
   if(!all(is.na(WGSres))){
      colnames(WGSres) = names(sizes)
      temp = melt(WGSres)
      temp$Var2 = factor(temp$Var2, levels=c("bins100bp", "bins1kb", 
                                             "bins10kb", "bins100kb", "bins1Mb"))
      p = ggplot(temp, aes(x = Var2, y = Var1)) + 
         geom_raster(aes(fill=value)) + 
         scale_fill_gradient(low="grey90", high="red") +
         ggtitle(tissue) +
         xlab("bin size") + 
         ylab("features") + 
         labs(fill="abs. Pearson's R") +
         theme(plot.title = element_text(hjust = 0.5, size=25),
               axis.title=element_text(size=20),
               axis.text=element_text(size=18), 
               axis.text.x=element_text(angle=30, hjust = 1),
               legend.text=element_text(size=18),
               legend.title=element_text(size=18))
      ggsave(p, filename = paste0(
         "fig/", tissue, "/compareCorrelationPerBinsize_", tissue, "_wholebinsWGS.png"), 
         width=10, height=10)
   }
   
   WXSres = sapply(res, function(x){x["nMuts_WXS",]})
   colnames(WXSres) = names(sizes)
   temp = melt(WXSres)
   temp$Var2 = factor(temp$Var2, levels=c("bins100bp", "bins1kb", 
                                          "bins10kb", "bins100kb", "bins1Mb"))
   p = ggplot(temp, aes(x = Var2, y = Var1)) + 
      geom_raster(aes(fill=value)) + 
      scale_fill_gradient(low="grey90", high="red") +
      ggtitle(tissue) +
      xlab("bin size") + 
      ylab("features") + 
      labs(fill="abs. Pearson's R") +
      theme(plot.title = element_text(hjust = 0.5, size=25),
            axis.title=element_text(size=20),
            axis.text=element_text(size=18), 
            axis.text.x=element_text(angle=30, hjust = 1),
            legend.text=element_text(size=18),
            legend.title=element_text(size=18))
   ggsave(p, filename = paste0(
      "fig/", tissue, "/compareCorrelationPerBinsize_", tissue, "_wholebinsWXS.png"), 
      width=10, height=10)
})
#####


# compare predictors between tissues and resolutions whole bins #####
# compare between tissues
print("compare predictors whole bin")
plotDat = sapply(names(sizes), function(size){
   print(size)
   res = array(NA,dim = c(length(tissues),
                          nrow(MutDensities[[size]]),
                          length(features)),
               dimnames=list(tissues, seq(nrow(MutDensities[[size]])),
                             features))
   for(tissue in tissues){
      load(paste0("data/rdata/", tissue, "/binsDat_", size, ".RData"))
      res[tissue,,colnames(dat)[-(1:5)]] <- as.matrix(dat[,-(1:5)])
   }
   MutsPerBin = MutDensities[[size]]
   thePlots = lapply(features, function(factor){
      png("fig/tempgraph.png",  
          width = 800, height = 1200, pointsize = 35)
      par(mfrow = c(length(tissues),1),
          mar = c(0.7,4,0,1),
          oma = c(2.5,1,3,1))
      xlim = c(0,tail(bins$end,1))
      xlabels = seq(0,250000000,by=25000000)
      cols = rainbow(length(tissues))
      names(cols) = tissues
      for(tissue in tissues){
         if(any(!is.na(res[tissue,,factor]))){
            plot(MutsPerBin$end, res[tissue,,factor], type = "h", xlim = xlim,
                 col = rgb(t(col2rgb(cols[tissue])),alpha=250, maxColorValue=255),
                 xlab = "Mb", ylab = "", xaxt = "n", las = 1,
                 yaxs = "i", lab = c(10,2,0), pch = 19, cex = 0.2)
         } else{
            plot(MutsPerBin$end, rep(NA, nrow(MutsPerBin)), ylim = c(0,1),
                 type = "h", xlim = xlim,
                 col = rgb(t(col2rgb(cols[tissue])),alpha=250, maxColorValue=255),
                 xlab = "Mb", ylab = "", xaxt = "n", las = 1,
                 yaxs = "i", lab = c(10,2,0), pch = 19, cex = 0.2)
         }
         legend("topleft", legend=tissue, bty="n")
         abline(v = xlabels,
                lty = 2, col = "grey")
      }
      axis(1, at=xlabels,labels=xlabels/1000000)
      mtext("genomic position (Mb)", side=1, line=2)
      mtext(paste0(factor, " at ", substr(size,5,10) , " scale"), 
            outer = T, line =1)
      dev.off()
      rasterGrob(readPNG("fig/tempgraph.png", native = FALSE),
                 interpolate = FALSE)
   })
   pdf(paste0("fig/genomePlot_",size,"_betweenTissues.pdf"))
   for(i in 1:length(thePlots)){
      do.call(grid.arrange, c(thePlots[i], ncol = 1))
   }
   dev.off()     
})
# compare between  resolutions
print("compare resolutions whole bin")
comp = sapply(tissues, function(tissue){
   print(tissue)
   thePlots = lapply(features, function(factor){
      png("fig/tempgraph.png", width=800, height=1200, pointsize=35)
      par(mfrow = c(length(sizes),1),
          mar = c(0.7,4,0,1),
          oma = c(2.5,1,3,1))
      xlim = c(0,tail(bins$end,1))
      xlabels = seq(0,250000000,by=25000000)
      cols = rainbow(length(sizes))
      names(cols) = names(sizes)
      for(size in names(sizes)){
         MutsPerBin = MutDensities[[size]]
         if(any(!is.na(plotDat[[size]][,factor]))){
            plot(MutsPerBin$end, plotDat[[size]][,factor], type = "h", xlim = xlim,
                 col = rgb(t(col2rgb(cols[size])),alpha=250, maxColorValue=255),
                 xlab = "Mb", ylab = "", xaxt = "n", las = 1,
                 yaxs = "i", lab = c(10,2,0), pch = 19, cex = 0.2)
         } else{
            plot(MutsPerBin$end, rep(NA, nrow(MutsPerBin)), ylim = c(0,1),
                 type = "h", xlim = xlim,
                 col = rgb(t(col2rgb(cols[size])),alpha=250, maxColorValue=255),
                 xlab = "Mb", ylab = "", xaxt = "n", las = 1,
                 yaxs = "i", lab = c(10,2,0), pch = 19, cex = 0.2)
         }
         legend("topleft", legend=size, bty="n")
         abline(v = xlabels,
                lty = 2, col = "grey")
      }
      axis(1, at=xlabels,labels=xlabels/1000000)
      mtext("genomic position (Mb)", side=1, line=2)
      mtext(paste0(factor, " in ", tissue), 
            outer = T, line =1)
      dev.off()
      rasterGrob(readPNG("fig/tempgraph.png", native = FALSE),
                 interpolate = FALSE)
   })
   pdf(paste0("fig/", tissue, "/genomePlot_", tissue,
              "_betweenResolutions.pdf"))
   for(i in 1:length(thePlots)){
      do.call(grid.arrange, c(thePlots[i], ncol = 1))
   }
   dev.off() 
})
#####

# compare predictors between tissues covered portions #####
# compare between tissues
print("compare predictors covered")
plotDat_covered = sapply(names(sizes), function(size){
   print(size)
   res = array(NA,dim = c(length(tissues),
                          nrow(MutDensities[[size]]),
                          length(features)),
               dimnames=list(tissues, seq(nrow(MutDensities[[size]])),
                             features))
   for(tissue in tissues){
      load(paste0("data/rdata/", tissue, "/binsDat_", size, "_covered.RData"))
      res[tissue,,colnames(dat)[-(1:4)]] = as.matrix(dat[,-(1:4)])
   }
   MutsPerBin = MutDensities[[size]]
   thePlots = lapply(features, function(factor){
      png("fig/tempgraph.png", width=800, height=1200, pointsize=35)
      par(mfrow = c(length(tissues),1),
          mar = c(0.7,4,0,1),
          oma = c(2.5,1,3,1))
      xlim = c(0,tail(bins$end,1))
      xlabels = seq(0,250000000,by=25000000)
      cols = rainbow(length(tissues))
      names(cols) = tissues
      for(tissue in tissues){
         if(any(!is.na(res[tissue,,factor]))){
            plot(MutsPerBin$end, res[tissue,,factor], type = "h", xlim = xlim,
                 col = rgb(t(col2rgb(cols[tissue])),alpha=250, maxColorValue=255),
                 xlab = "Mb", ylab = "", xaxt = "n", las = 1,
                 yaxs = "i", lab = c(10,2,0), pch = 19, cex = 0.2)
         } else{
            plot(MutsPerBin$end, rep(NA, nrow(MutsPerBin)), ylim = c(0,1),
                 type = "h", xlim = xlim,
                 col = rgb(t(col2rgb(cols[tissue])),alpha=250, maxColorValue=255),
                 xlab = "Mb", ylab = "", xaxt = "n", las = 1,
                 yaxs = "i", lab = c(10,2,0), pch = 19, cex = 0.2)
         }
         legend("topleft", legend=tissue, bty="n")
         abline(v = xlabels,
                lty = 2, col = "grey")
      }
      axis(1, at=xlabels,labels=xlabels/1000000)
      mtext("genomic position (Mb)", side=1, line=2)
      mtext(paste0(factor, " at ", substr(size,5,10) , " scale"), 
            outer = T, line =1)
      dev.off()
      rasterGrob(readPNG("fig/tempgraph.png", native = FALSE),
                 interpolate = FALSE)
   })
   pdf(paste0("fig/genomePlot_covered_",size,"_betweenTissues.pdf"))
   for(i in 1:length(thePlots)){
      do.call(grid.arrange, c(thePlots[i], ncol = 1))
   }
   dev.off()    
})
# compare different resolutions 
print("compare resolutions covered")
comp_covered = sapply(tissues, function(tissue){
   print(tissue)
   thePlots = lapply(features, function(factor){
      png("fig/tempgraph.png", width=800, height=1200, pointsize=35)
      par(mfrow = c(length(sizes),1),
          mar = c(0.7,4,0,1),
          oma = c(2.5,1,3,1))
      xlim = c(0,tail(bins$end,1))
      xlabels = seq(0,250000000,by=25000000)
      cols = rainbow(length(sizes))
      names(cols) = names(sizes)
      for(size in names(sizes)){
         MutsPerBin = MutDensities[[size]]
         if(any(!is.na(plotDat_covered[[size]][,factor]))){
            plot(MutsPerBin$end, plotDat_covered[[size]][,factor], type = "h", xlim = xlim,
                 col = rgb(t(col2rgb(cols[size])),alpha=250, maxColorValue=255),
                 xlab = "Mb", ylab = "", xaxt = "n", las = 1,
                 yaxs = "i", lab = c(10,2,0), pch = 19, cex = 0.2)
         } else{
            plot(MutsPerBin$end, rep(NA, nrow(MutsPerBin)), ylim = c(0,1),
                 type = "h", xlim = xlim,
                 col = rgb(t(col2rgb(cols[size])),alpha=250, maxColorValue=255),
                 xlab = "Mb", ylab = "", xaxt = "n", las = 1,
                 yaxs = "i", lab = c(10,2,0), pch = 19, cex = 0.2)
         }
         legend("topleft", legend=size, bty="n")
         abline(v = xlabels,
                lty = 2, col = "grey")
      }
      axis(1, at=xlabels,labels=xlabels/1000000)
      mtext("genomic position (Mb)", side=1, line=2)
      mtext(paste0(factor, " in ", tissue), 
            outer = T, line =1)
      dev.off()
      rasterGrob(readPNG("fig/tempgraph.png", native = FALSE),
                 interpolate = FALSE)
   })
   pdf(paste0("fig/", tissue, "/genomePlot_covered_", tissue,
              "_betweenResolutions.pdf"))
   for(i in 1:length(thePlots)){
      do.call(grid.arrange, c(thePlots[i], ncol = 1))
   }
   dev.off() 
})
#####


# compare approxfun, smooth.spline and splinefun on GCcontent #####
png("fig/compareApproxfunSmoothSplinefun_luadGCcontent.png", 
    width=800, height=1200, pointsize=30)
par(mfrow = c(5,1), mar = c(3,4,1,1))
for (size in names(sizes)){
   a = approxfun(MutDensities[[size]]$start, GCcontents[[size]])
   b = smooth.spline(MutDensities[[size]]$start[!is.na(GCcontents[[size]])],
                     GCcontents[[size]][!is.na(GCcontents[[size]])])
   # predict(b,MutDensities[[size]]$start)
   d = splinefun(MutDensities[[size]]$start[!is.na(GCcontents[[size]])],
                 GCcontents[[size]][!is.na(GCcontents[[size]])])
   xmax = min(250000000,sizes[size]*1000)
   plot(MutDensities[[size]]$start, GCcontents[[size]],
        pch = 19, cex = 0.9, xlim = c(0,xmax), col = rgb(0,0,0,0.6),
        xlab = "genomic position", ylab = "GCcontent", main = size)
   curve(a, col = 2, lty = 1, lwd = 3, add = T)
   lines(b, col = 3, lty = 2, lwd = 3)
   curve(d, col = 4, lty = 3, lwd = 3, add = T)
   legend("bottomright", 
          col = 2:4, lty = 1:3, lwd = 3,
          legend=c("approxfun", "smooth.spline", "splinefun"))
}
dev.off()
# splinefun doesnt handle outliers well and it interpolates 
# beyond the ends, otherwise it is the same as approxfun 
# --> don't use it.


# compare approxfun, smooth.spline and splinefun on +H3K27ac+ #####
tissue = "luad"
png("fig/compareApproxfunSmoothSplinefun_luadH3K27ac.png", 
    width=800, height=1200, pointsize=30)
par(mfrow = c(5,1), mar = c(3,4,1,1))
for (size in names(sizes)){
   load(paste0("data/procData/bins/", size,"/", tissue,
               "/DNAbinding/DNAbinding.RData"))
   x = MutDensities[[size]]$start
   y = DNAbinding[,"H3K27ac"]
   xmax = min(250000000,sizes[size]*1000)
   y = y[x<xmax]
   x = x[x<xmax]
   a = approxfun(x, y)
   b = smooth.spline(x[!is.na(y)],y[!is.na(y)])
   # predict(b,MutDensities[[size]]$start)
   d = splinefun(x[!is.na(y)],y[!is.na(y)])
   plot(x, y,  ylab = "H3K27ac",
        pch = 19, cex = 0.9,  col = rgb(0,0,0,0.6),
        xlab = "genomic position",main = size)
   curve(a, col = 2, lty = 1, lwd = 3, add = T)
   lines(b, col = 3, lty = 2, lwd = 3)
   curve(d, col = 4, lty = 3, lwd = 3, add = T)
   legend("topleft", 
          col = 2:4, lty = 1:3, lwd = 3,
          legend=c("approxfun", "smooth.spline", "splinefun"))
}
dev.off()
# here, smooth.spline creates problems because it oscillates 
# unnecessarily at the higher resolutions when lots of datapoints are 0.


# compare approxfun, smooth.spline and splinefun on UCSC_H3k4me3 #####
tissue = "luad"
png("fig/compareApproxfunSmoothSplinefun_luadUCSCH3k4me3.png", 
    width=800, height=1200, pointsize=30)
par(mfrow = c(5,1), mar = c(3,4,1,1))
for (size in names(sizes)){
   load(paste0("data/procData/bins/", size, "/", tissue,
               "/UCSC_tracks/UCSC_histones.RData"))
   x = MutDensities[[size]]$start
   y = histones[,"UCSC_H3k4me3"]
   xmax = min(250000000,sizes[size]*1000)
   y = y[x<xmax]
   x = x[x<xmax]
   a = approxfun(x, y)
   b = smooth.spline(x[!is.na(y)],y[!is.na(y)])
   # predict(b,MutDensities[[size]]$start)
   d = splinefun(x[!is.na(y)],y[!is.na(y)])
   plot(x, y,  ylab = "UCSC_H3k4me3",
        pch = 19, cex = 0.9,  col = rgb(0,0,0,0.6),
        xlab = "genomic position",main = size)
   curve(a, col = 2, lty = 1, lwd = 3, add = T)
   lines(b, col = 3, lty = 2, lwd = 3)
   curve(d, col = 4, lty = 3, lwd = 3, add = T)
   legend("topleft", 
          col = 2:4, lty = 1:3, lwd = 3,
          legend=c("approxfun", "smooth.spline", "splinefun"))
}
dev.off()
##### CONCLUSION: use approxfun!




# use approxfun to represent predictors #####
# iterate through tissues, then sizes, then features to save all approxfunctions
# whole bins
temp = sapply(tissues, function(tissue){
   print(tissue)
   if(!dir.exists(paste0("data/rdata/", tissue, "/approx/"))){
      dir.create(paste0("data/rdata/", tissue, "/approx/"))
   }
   res = lapply(names(sizes), function(size){
      cat(size,' ')
      load(paste0("data/rdata/", tissue, "/binsDat_", size, ".RData"))
      app = sapply(features[features %in% colnames(dat)],function(feature){
         x = c(dat$start[1], 
               rowMeans(dat[,c("start","end")]), 
               tail(dat$start,n=1))
         y = c(dat[1,feature], 
               dat[,feature], 
               dat[nrow(dat),feature])
         if(all(is.na(y))){return(NA)}
         a = approxfun(x,y)
         save(a, file = paste0("data/rdata/", tissue, 
                               "/approx/",feature, "_", size,".RData" ))
         return(a)
      })
   })
   cat('\n')
})

# covered parts
temp = sapply(tissues, function(tissue){
   print(tissue)
   res = lapply(names(sizes), function(size){
      cat(size,' ')
      load(paste0("data/rdata/", tissue, "/binsDat_",
                  size, "_covered.RData"))
      app = sapply(features[features %in% colnames(dat)],function(feature){
         x = c(dat$start[1], 
               rowMeans(dat[,c("start","end")]), 
               tail(dat$start,n=1))
         y = c(dat[1,feature], 
               dat[,feature], 
               dat[nrow(dat),feature])
         if(all(is.na(y))){return(NA)}
         a = approxfun(x,y)
         save(a, file = paste0("data/rdata/", tissue, 
                               "/approx/",feature, "_", size,"_covered.RData" ))
         return(a)
      })
   })
   cat('\n')
})
#####


# complete plot: on page per tissue, one subpanel per feature, #####
# different resolutions as different colors in the same plot
cols = c("bins1Mb" =2, "bins100kb"=3, "bins10kb"=4, 
         "bins1kb"=5, "bins100bp"=6)
# whole bins
temp = sapply(tissues, function(tissue){
   print(tissue)
   pdf(paste0("fig/compareSplineResolutions_", tissue, ".pdf"), 
       height=8, width=12)
   sapply(features, function(feature){
      if(file.exists(paste0("data/rdata/", tissue, "/approx/",
                            feature, "_", "bins1Mb",".RData" ))){
         sapply(names(sizes[5:1]), function(size){
            load(paste0("data/rdata/", tissue, 
                        "/approx/",feature, "_", size,".RData" ))
            curve(a, add = !size==names(sizes[5]), 
                  from=0, to = 250000000, 
                  col = cols[size], n=10000, ylab = feature, 
                  xlab = "position", main= feature)
         })
         legend("topright", legend=names(cols), col=cols,
                lty =1, bg="white")
      }
   })
   dev.off()
})

# covered part
temp = sapply(tissues, function(tissue){
   print(tissue)
   pdf(paste0("fig/compareSplineResolutionsCovered_", tissue, ".pdf"),
       height=8, width=12)
   sapply(features, function(feature){
      if(file.exists(paste0("data/rdata/", tissue, 
                            "/approx/",feature, "_", 
                            "bins1Mb_covered.RData" ))){
         sapply(names(sizes[5:1]), function(size){
            load(paste0("data/rdata/", tissue, "/approx/",
                        feature, "_", size,"_covered.RData" ))
            curve(a, add = !size==names(sizes[5]),
                  from=0, to = 250000000, 
                  col = cols[size], n=10000, ylab = feature, 
                  xlab = "position", main= feature)
         })
         legend("topright", legend=names(cols), col=cols, 
                lty =1, bg="white")
      }
   })
   dev.off()
})
#####


# get binwise values for TN and TP positions ######
# whole bins
temp = sapply(tissues, function(tissue){
   load(paste0("data/rdata/", tissue, "/completeData.RData"))
   data = data[removed$chr == "chr1",]
   removed = removed[removed$chr == "chr1",]
   Muts = data.frame(pos = removed$pos, mutated = as.integer(data$mutated))
   temp2 = sapply(names(sizes), function(size){
      load(paste0("data/rdata/", tissue,
                  "/binsDat_", size, ".RData"))
      fun = approxfun(dat$start, seq(nrow(dat)), method = "constant")
      inds = fun(Muts$pos)
      binwise = dat[inds, -(1:5)]
      save(binwise, file = paste0("data/rdata/", tissue,
                                  "/binwisePreds_", size, ".RData"))
   })
})
# covered parts
temp = sapply(tissues, function(tissue){
   load(paste0("data/rdata/", tissue, "/completeData.RData"))
   data = data[removed$chr == "chr1",]
   removed = removed[removed$chr == "chr1",]
   Muts = data.frame(pos = removed$pos, mutated = as.integer(data$mutated))
   temp2 = sapply(names(sizes), function(size){
      load(paste0("data/rdata/", tissue,
                  "/binsDat_", size, "_covered.RData"))
      fun = approxfun(dat$start, seq(nrow(dat)), method = "constant")
      inds = fun(Muts$pos)
      binwise = dat[inds, -(1:4)]
      save(binwise, file = paste0("data/rdata/", tissue,
                                  "/binwisePreds_", size, "_covered.RData"))
   })
})
#####


# compare correlation #####
pdf("fig/comparePredictorTypes.pdf")
# per feature one plot
sapply(features, function(feature){
   cat(feature, ' ')
   # per tissue one panel
   temp = lapply(tissues, function(tissue){
      load(paste0("data/rdata/", tissue, "/completeData.RData"))
      data = data[removed$chr == "chr1",]
      removed = removed[removed$chr == "chr1",]
      Muts = data.frame(pos = removed$pos, mutated = as.integer(data$mutated))
      # barplot of correlations
      cors = rbind(
         # 1. bin-wise readout
         data.frame(x="binwise",resolution = substr(names(sizes),5,10),
                    z=sapply(names(sizes), function(size){
                       load(paste0("data/rdata/", tissue,
                                   "/binwisePreds_", size, ".RData"))
                       if(feature %in% colnames(binwise)){
                          cor(Muts[,"mutated"],binwise[,feature], use = "pair")
                       } else
                          return(NA)
                    })),
         # 2. bin-wise readout covered
         data.frame(x="binwise_covered",resolution = substr(names(sizes),5,10),
                    z=sapply(names(sizes), function(size){
                       load(paste0("data/rdata/", tissue,
                                   "/binwisePreds_", size, "_covered.RData"))
                       if(feature %in% colnames(binwise)){
                          cor(Muts[,"mutated"],binwise[,feature], use = "pair")
                       } else
                          return(NA)
                    })), 
         # 3. spline
         data.frame(x="approx",resolution = substr(names(sizes),5,10),
                    z=if(file.exists(paste0("data/rdata/", tissue, 
                                            "/approx/",feature, "_", "bins1Mb",".RData" ))){
                       sapply(names(sizes), function(size){
                          print(size)
                          load(paste0("data/rdata/", tissue, 
                                      "/approx/",feature, "_", size,".RData" ))
                          cor(Muts[,"mutated"],a(Muts[,"pos"]), use = "pair")
                       })
                    } else {
                       rep(NA,length(sizes))
                    }),
         # 4. spline covered
         data.frame(x = "approx_covered", resolution = substr(names(sizes),5,10), 
                    z = if(file.exists(paste0("data/rdata/", tissue, 
                                              "/approx/",feature, "_", "bins1Mb_covered.RData" ))){
                       sapply(names(sizes), function(size){
                          load(paste0("data/rdata/", tissue, 
                                      "/approx/",feature, "_", size,"_covered.RData" ))
                          cor(Muts[,"mutated"],a(Muts[,"pos"]), use = "pair")
                       })
                    } else {
                       rep(NA,length(sizes))
                    }),
         # 5. single base readout and, if available, direct readout around feature
         data.frame(x="baseReadout",
                    resolution = if(feature == "GCcontent"){
                       do.call(rbind,
                               strsplit(colnames(data)[grep(feature, colnames(data))],
                                        split="_"))[,2]
                    } else{"1bp"},
                    z=if(length(grep(feature, colnames(data)))>0){
                       as.vector(cor(Muts[,"mutated"], 
                                     data[,grep(feature, colnames(data))],
                                     use = "pair"))}
                    else{NA})
      )
      cors$resolution = factor(cors$resolution,
                               levels = levels(cors$resolution)[c(7,6,1,4,3,2,5)])
      cors$tissue = tissue
      return(cors)
   })
   temp2 = do.call(rbind,lapply(temp,function(x){
      y = unlist(sapply(split(x,x$x), function(i){
         t = rep("",nrow(i))
         t[which.max(abs(i$z))] = "group max"
         t
      }))
      y[which.max(abs(x$z))] = "total max"
      x$max = factor(y)
      x
   }))
   temp2$tissue = factor(temp2$tissue)
   g = ggplot(data=temp2,aes(x=x, y=z,fill = resolution, color = max)) +
      geom_bar(stat="identity", position=position_dodge2(preserve="single"), 
               size = 0.5) +
      scale_color_manual(values = c("transparent","grey33" ,"black")) +
      ggtitle(feature) +
      facet_wrap( ~ tissue, ncol = 1) +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5, size=25))+#,
            # axis.title=element_text(size=15),
            # axis.text=element_text(size=15),
            # legend.text=element_text(size=15),
            # legend.title=element_text(size=15))+
      ylab("correlation with TP/TN labels") +
      xlab("")
      # guides(color=F)   
   ggsave("fig/binComparison/comparePredictorTypes_", feature,".png")
   print(g)
})
dev.off()
cat('\n')
######


# compare RF prediction performance between predictor types #####
library(ranger)
pdf("fig/comparePredictorTypes_RFpredError_scaled.pdf", height=10,width=10)
RFimportances = lapply(tissues, function(tissue){
   print(tissue)
   load(paste0("data/rdata/", tissue, "/completeData.RData"))
   data = data[removed$chr == "chr1",]
   removed = removed[removed$chr == "chr1",]
   Muts = data.frame(pos = removed$pos, mutated = as.integer(data$mutated))
   par(mfrow = c(5,5), mar = c(5,4,2,0), oma = c(2,1,2,0))
   importances = sapply(features, function(feature){
      res = cbind(
         # 1. bin-wise readout
         sapply(names(sizes), function(size){
            load(paste0("data/rdata/", tissue,
                        "/binwisePreds_", size, ".RData"))
            if(feature %in% colnames(binwise)){
               binwise[,feature]
            } else
               rep(NA,nrow(binwise))
         }),
         # 5. single base readout and, if available, direct readout around feature
         if(feature == "GCcontent"){
            data[,grep(feature, colnames(data))]
         } else if (feature == "methylation"){
            data[,grep("meth", colnames(data))]
         } else if (feature == "DNAaccessibility"){
            data[,grep("DNAaccessibility_tissue", colnames(data))]
         } else if (feature == "UCSC_DNAaccessibility"){
            data[,"DNAaccessibility_UCSC"]
         } else if (feature %in% colnames(data)){
            data[,feature]
         } else {NA}
      )
      colnames(res)[colnames(res) == ""] = "signal1bp"
      notMissing = apply(res, 1, function(x){!any(is.na(x))})
      if(sum(notMissing)<=5){
         barplot(NA, ylim = c(0,10), las = 2, main = feature)
         return(rep(NA, length(sizes)))
      }
      dat = data.frame(y =Muts$mutated[notMissing],
                       res[notMissing,apply(res, 2, function(x){!all(is.na(x))})])
      rf = ranger(y~.,data=dat,importance='permutation', holdout=T, 
                  scale.permutation.importance=T,
                  case.weights=sample(0:1,nrow(dat),replace=T,c(0.3,0.7)))
      barplot(rf$variable.importance, main=feature, las=2,
              col = c(rainbow(5),grey.colors(length(rf$variable.importance)-5)))
      return(rf$variable.importance[1:5])
   })   
   mtext(tissue, outer = T, line = -0.5, font=2, cex = 1.3)
   mtext("RF permutation importance", outer = T, line = -0.5, side = 2)
   return(importances)
})
dev.off()
bestResolution_RFimportance = sapply(features, function(feature){
   best = unlist(sapply(RFimportances, function(x){
      which.max(x[,feature])
   }))
   table(factor(names(best),levels=names(sizes)))
})
##### 


# compare lm prediction performance #####
pdf("fig/comparePredictorTypes_LMsingificance.pdf", height=10,width=10)
LMpvals = lapply(tissues, function(tissue){
   print(tissue)
   load(paste0("data/rdata/", tissue, "/completeData.RData"))
   data = data[removed$chr == "chr1",]
   removed = removed[removed$chr == "chr1",]
   Muts = data.frame(pos = removed$pos, mutated = as.integer(data$mutated))
   par(mfrow = c(5,5), mar = c(5,4,2,0), oma = c(2,1,2,0))
   pvals = sapply(features, function(feature){
      res = cbind(
         # 1. bin-wise readout
         sapply(names(sizes), function(size){
            load(paste0("data/rdata/", tissue,
                        "/binwisePreds_", size, ".RData"))
            if(feature %in% colnames(binwise)){
               binwise[,feature]
            } else
               rep(NA,nrow(binwise))
         }),
         # 5. single base readout and, if available, direct readout around feature
         if(feature == "GCcontent"){
            data[,grep(feature, colnames(data))]
         } else if (feature == "methylation"){
            data[,grep("meth", colnames(data))]
         } else if (feature == "DNAaccessibility"){
            data[,grep("DNAaccessibility_tissue", colnames(data))]
         } else if (feature == "UCSC_DNAaccessibility"){
            data[,"DNAaccessibility_UCSC"]
         } else if (feature %in% colnames(data)){
            data[,feature]
         } else {NA}
      )
      colnames(res)[colnames(res) == ""] = "signal1bp"
      notMissing = apply(res, 1, function(x){!any(is.na(x))})
      if(sum(notMissing)<=5){
         barplot(NA, ylim = c(0,10), las = 2, main = feature)
         return(rep(NA, length(sizes)))
      }
      dat = data.frame(y =Muts$mutated[notMissing],
                       res[notMissing,apply(res, 2, function(x){!all(is.na(x))})])
      model = lm(y~.,data=dat)
      barplot(-log10(summary(model)$coefficients[-1,4]), main=feature, las=2,
              col = c(rainbow(5),grey.colors(ncol(dat)-6)))
      return(-log10(summary(model)$coefficients[-1,4])[1:5])
   })
   mtext(tissue, outer = T, line = -0.5, font=2, cex = 1.3)
   mtext("-log10(p-value) in linear model", side=2,
         outer = T, line = -0.5, cex = 1)
   return(pvals)
})
dev.off()
bestResolution_LMpvals = sapply(features, function(feature){
   best = unlist(sapply(LMpvals, function(x){
      which.max(x[,feature])
   }))
   table(factor(names(best),levels=names(sizes)))
})
#####



# input all features at all resolutions into RF #####
library(ranger)
pdf("fig/comparePredictorTypes_RFpredErrorAllFeatures.pdf",
    height=12,width=12)
temp2 = lapply(tissues, function(tissue){
   print(tissue)
   load(paste0("data/rdata/", tissue, "/completeData.RData"))
   data = data[removed$chr == "chr1",]
   removed = removed[removed$chr == "chr1",]
   Muts = data.frame(pos = removed$pos, 
                     mutated = as.integer(data$mutated))
   predictors = sapply(features, function(feature){
      res = cbind(
         # 1. bin-wise readout
         sapply(names(sizes), function(size){ # [c(2,5)]
            load(paste0("data/rdata/", tissue,
                        "/binwisePreds_", size, ".RData"))
            if(feature %in% colnames(binwise)){
               binwise[,feature]
            } else
               return(rep(NA, nrow(data)))
         }),
         # 5. single base readout and, if available, direct readout around feature
         if(feature == "GCcontent"){
            # data[,grep(feature, colnames(data))]
            data[,"GCcontent_5bp"]
         } else if (feature == "methylation" & tissue != "ovary"){
            data[,"methbank_1kbp"]
         } else if (feature == "DNAaccessibility"){
            data[,"DNAaccessibility_tissue"]
         } else if (feature == "UCSC_DNAaccessibility"){
            data[,"DNAaccessibility_UCSC"]
         } else if (feature %in% colnames(data)){
            data[,feature]
         } else{NA}
      )
      colnames(res)[colnames(res) == ""] = "signal1bp"
      colnames(res) = paste(feature, colnames(res), sep=":")
      return(res)
   }, simplify = F)
   predictors = do.call(cbind, predictors)
   predictors = predictors[,apply(predictors,2,function(x){!all(is.na(x))})]
   notMissing = apply(predictors, 1, function(x){!any(is.na(x))})
   dat = data.frame(y = Muts$mutated[notMissing],
                    predictors[notMissing,])
   rf = ranger(y~.,data=dat,importance='permutation',
               holdout=T, scale.permutation.importance=T,
               case.weights=sample(0:1,nrow(dat),replace=T,c(0.3,0.7)))
   labels = sapply(colnames(dat)[-1],function(x){
      strsplit(x,".", fixed = T)[[1]]
   })
   plotData = data.frame(t(labels),rf$variable.importance)
   plotData$X1 = factor(plotData$X1,levels=features)
   plotData$X2 = factor(plotData$X2,levels=c(names(sizes), 
                                             "signal1bp"))
   plotData$column = plotData$X1 %in% levels(plotData$X1)[1:(length(levels(plotData$X1))/2)]
   g = ggplot(data=plotData,aes(y=rf.variable.importance, 
                                x=X1,fill = X2)) +
      geom_bar(stat="identity", 
               position=position_dodge2(preserve="single"), 
               size = 0.5) +
      coord_flip() +
      facet_wrap(~column,scales="free_y", 
                 labeller = as_labeller(c("TRUE" = "", "FALSE" = ""))) +
      ggtitle(tissue) +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5, size=25),
            axis.text=element_text(size=12))+
      labs(fill = "Binsize") +
      ylab("RF permutation importance") +
      xlab("")
   ggsave(paste0("fig/comparePredictorTypes_RFpredErrorAllFeatures_",
                 tissue, ".png"),plot=g, scale=1.2)
   print(g)
})
dev.off()
#####


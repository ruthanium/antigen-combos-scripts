# *******************************************
# Author: Ruth Dannenfelser
# Date: October 2019
#
# make all the scatterplots in the paper
# *******************************************

library(ggrepel)
library(data.table)
library(ggplot2)
library(plotly)


rnaseqmat = readRDS("scripts/clustering-method/rdata/combat-rnaseq-tm-clin-zeroadj.RData")
sketchLAll = rbindlist(lapply(cans,
                              function(can) data.table(can=can,
                                                       dataset=readLines(paste0("scripts/clustering-method/sketches/", gsub(" ", "-", tolower(can)), "-sketch.txt")))))


# Figure 4 - contour plots and scatter plots for supplementary figure
genes = c("PRR7", "CLDN6")
target = "Ovarian Serous Cystadenocarcinoma"
mSP(target, genes, T)
makePlotsWithoutSaving(target, genes, T)
makePlotsWithoutSavingHighlights(target, genes, F, c("Ovary"))

genes = c("B4GALNT1", "MSLN")
target = "Mesothelioma"
makePlotsWithoutSaving(target, genes, T)
mSP(target, genes, T)
mSP.withCutoff(target, genes, 9, T)
#makePlotsWithoutSavingHighlights(target, genes, F, c("Pancreas"))


genes = c("B4GALNT1", "KREMEN2")
target = "Head and Neck Squamous Cell Carcinoma"
mSP.withCutoff(target, genes, 7.5, T)
mSP(target, genes, T)
makePlotsWithoutSaving(target, genes, T)

genes = c("TREM1", "CA9")
target = "Lung Adenocarcinoma"
mSP.withCutoff(target, genes, 9, T)
mSP(target, genes, T)
makePlotsWithoutSaving(target, genes, T)

genes = c("FAP", "MUC16")
target = "Pancreatic Adenocarcinoma"
makePlots(target, genes)

target = "Brain Lower Grade Glioma"
genes = c("MOG", "CACNG7")
makePlotsWithoutSaving(target, genes, T)

target = "Glioblastoma Multiforme"
genes = c("MOG", "IL13RA2")
mSP.withCutoff(target, genes, 7.5, T)
mSP(target, genes, T)
makePlotsWithoutSaving(target, genes, T)
makePlotsWithoutSavingHighlights(target, genes, F, c("Brain"))


target = "Kidney Renal Clear Cell Carcinoma"
genes = c("AXL", "CDH6")
mSP(target, genes, T)
mSP.highlight(target, genes, T, c("Lung"))
makePlotsWithoutSaving(target, genes, T)
makePlotsWithoutSavingHighlights(target, genes, F, c("Lung"))

genes = c("CD70", "AXL")
makePlotsWithoutSaving(target, genes, T)
mSP(target, genes, T)
mSP.highlight(target, genes, T, c("Blood"))
#makePlotsWithoutSavingHighlights(target, genes, F, c("Blood"))

genes = c("CDH6", "CD70")
makePlotsWithoutSaving(target, genes, T)


genes = c("MSLN", "FAP")
target = "Pancreatic Adenocarcinoma"

mSP.highlight(target, genes, T, c("Lung"))
mSP.highlight.withCutoff(target, genes, 7.5, T, c("Lung"))
mSP.withCutoff(target, genes, 7.5, T)
mSP(target, genes, T)

makePlotsWithoutSaving(target, genes, T)
makePlotsWithoutSavingHighlights(target, genes, F, c("Lung"))

target = "Kidney Renal Clear Cell Carcinoma"
genes = c("SEMA5B", "APLN")
mSP(target, genes, T)
makePlotsWithoutSaving(target, genes, T)


#potential bad examples - for wendell
genes = c("EGFR", "MSLN")
target = "Pancreatic Adenocarcinoma"
makePlotsWithoutSaving(target, genes)

genes = c("EPCAM", "MSLN")
target = "Pancreatic Adenocarcinoma"
makePlotsWithoutSaving(target, genes)


# first set that Wendell doesn't really like
# making the best, med, bad plots
genes = c("CD180", "CXCR5") # best overall
target = "Lymphoid Neoplasm Diffuse Large B-cell Lymphoma"
makePlotsWithoutSaving(target, genes)

genes = c("LRCH3", "TTYH2") # worst overall
makePlotsWithoutSaving(target, genes)

genes = c("GKN1", "CXCR5") # medium potentially?
makePlotsWithoutSaving(target, genes)

# second set for the good medium and bad
target = "Kidney Renal Clear Cell Carcinoma"
genes = c("GABRD", "AXL") # best
makePlotsWithoutSaving(target, genes, T)

genes = c("CALY", "AXL") # worst
makePlotsWithoutSaving(target, genes, T)

genes = c("ATP8B3", "AXL") # medium-ish
makePlotsWithoutSaving(target, genes, T)

target = "Pancreatic Adenocarcinoma"
genes = c("PTPRN", "MSLN") # best
makePlotsWithoutSaving(target, genes, T)

genes = c("KRT5", "MSLN") # worst
makePlotsWithoutSaving(target, genes, T)

genes = c("DLK1", "MSLN") # medium-ish
makePlotsWithoutSaving(target, genes, T)

# extras that are H:L
target = "Brain Lower Grade Glioma"
genes = c("EPCAM", "SYT11")
mSP(target, genes)
mSP.withCutoff(target, genes, 9, T)

target = "Cervical Squamous Cell Carcinoma and Endocervical Adenocarcinoma"
genes = c("KREMEN2", "KDR")
makePlotsWithoutSaving(target, genes, T)

target="Kidney Renal Clear Cell Carcinoma"
genes = c("VIPR1", "CA9")
mSP(target, genes)
mSP.withCutoff(target, genes, 9, T)

target="Liver Hepatocellular Carcinoma"
genes = c("MUC1", "PAQR9")
makePlotsWithoutSaving(target, genes, T)

target="Pheochromocytoma and Paraganglioma"
genes = c("CELSR1", "L1CAM")
makePlotsWithoutSaving(target, genes, T)

# experimenting with 3D scatter plots with plotly
target = "Skin Cutaneous Melanoma"
genes = c("ROR1", "MUC1", "EPCAM") # best in test case
make3DPlotsWithoutSaving(target, genes, T)

genes = c("EPHA2","ROR1", "FAP") # worst in test case
make3DPlotsWithoutSaving(target, genes, T)

target = "Acute Myeloid Leukemia"
genes = c("EGFR", "MET", "ERBB2") # best overall
make3DPlotsWithoutSaving(target, genes, T)

target = "Uterine Carcinosarcoma"
genes = c("EPCAM", "FOLH1", "GPC3") # worst overall


genes = c("FAP", "MSLN", "CTRB2")
target = "Pancreatic Adenocarcinoma"
make3DPlotsWithoutSaving(target, genes, T)

make3DPlots(target, genes, T)

genes = c("FAP", "MSLN", "CTRB1")
make3DPlotsWithoutSaving(target, genes, T)

genes = c("FAP", "MSLN", "CLPS")
make3DPlotsWithoutSaving(target, genes, T)

genes = c("FAP", "MSLN", "PTPRN")
make3DPlotsWithoutSaving(target, genes, T)

genes = c("FAP", "MSLN", "SST")
make3DPlotsWithoutSaving(target, genes, T)


genes = c("FAP", "MSLN", "CLDN18")
make3DPlotsWithoutSaving(target, genes, T)

genes = c("FAP", "MSLN", "CEACAM5")
make3DPlotsWithoutSaving(target, genes, T)

genes = c("FAP", "MSLN", "CA9")
make3DPlotsWithoutSaving(target, genes, T)

# supplemental plots
# 3 C:C
target = "Kidney Renal Clear Cell Carcinoma"
genes = c("KDR", "FOLH1")
mSP.withCutoff(target, genes, 7.5, T)
genes = c("MET", "CD70")
mSP(target, genes, T)

target = "Ovarian Serous Cystadenocarcinoma"
genes = c("MSLN", "L1CAM")
mSP.withCutoff(target, genes, 9, T)

# 3 C:Ns
target = "Pheochromocytoma and Paraganglioma" # yes
genes = c("EGFR", "IL13RA2")
mSP(target, genes, T)
genes = c("EPHA8", "B4GALNT1") # YES!
mSP(target, genes, T)

target = "Mesothelioma" # yes
genes = c("CA9", "KISS1R")
mSP.withCutoff(target, genes, 9, T)

# 6 N:Ns
target = "Testicular Germ Cell Tumors"
genes = c("SYT2", "CALHM3")
mSP.withCutoff(target, genes, 6, T)
genes = c("SYT2", "STRA6")
mSP.withCutoff(target, genes, 6, T)

target = "Sarcoma"
genes = c("CLDN4", "MMP14")
mSP(target, genes, T)

target = "Breast Invasive Carcinoma"
genes = c("ADAM12", "GPR133")
mSP(target, genes, T)

target = "Colon Adenocarcinoma"  
genes = c("LIFR", "GRIN2D")
mSP(target, genes, T)

target = "Stomach Adenocarcinoma"
genes = c("PRR7", "MUC17")
mSP.withCutoff(target, genes, 8, T)
genes = c("CHRNG", "KCNK9") # maybe but didn't use yet


target = "Pancreatic Adenocarcinoma"
genes = c("MSLN", "GRIN2D")
mSP.withCutoff(target, genes, 6.5, T)

target = "Mesothelioma"
genes = c("CA9", "KISS1R")
mSP.withCutoff(target, genes, 7.5,T)

target = "Kidney Renal Clear Cell Carcinoma"
genes = c("VIPR1","CA9")
mSP.withCutoff(target, genes, 8, T)

target = "Stomach Adenocarcinoma"
genes = c("CACNG8","CEACAM5")
mSP.withCutoff(target, genes, 8, T)


res[clin==T,]
target = "Pancreatic Adenocarcinoma"
apply (head(res, 5), 1, function(row) { 
  genes = unlist(strsplit(as.character(row["pair"]), ":"))  
  make3DPlotsWithoutSaving(target, genes, T)
})


# make plot out of the melanoma data
mela = fread("data/melanoma-ribas-plot-data.txt", header=T, sep="\t")
mela$MET = log(mela$MET)
mela$MLANA = log(mela$MLANA)
mela[,col:= ifelse(mela$highlight != "blah", "hilite", "blah")]
mela$MET = ifelse(mela$MET < 0, 0, mela$MET)
mela$MLANA = ifelse(mela$MLANA < 0, 0, mela$MLANA)
mela[,col:= ifelse(mela$MET > 1 & mela$MLANA > 2, "hilite", col)]

m = max(mela$MET, mela$MLANA)

mm = ggplot(mela, aes(x=MET, y=MLANA)) + geom_point(aes(color=col), alpha = 1, size=2) +
  scale_color_manual(values=c("#85919e", "#df9c95")) + 
  geom_text_repel(data=mela[highlight != "blah",], aes(label=highlight), size=4, box.padding = 0.25, 
                  point.padding = 0.2, segment.color = 'grey50', segment.size = 0.25) + theme_bw() +
  ylab(paste("MLANA", "log(TPM)")) + xlab(paste("MET", "log(TPM)")) + xlim(0,4) +
  theme(axis.title = element_text(size=18), axis.text= element_text(size=16), panel.grid.minor=element_blank(),
        plot.background=element_blank(), panel.background=element_blank(), panel.grid.major=element_blank()) +
  guides(colour=F)

ggsave(paste0("scripts/clustering-method/res-figures/paper/", genes[1], "-", genes[2], "-", target, "-ribas.pdf"), 
       mm, width=3.5, height=3.5, units="in", device="pdf", useDingbats=FALSE)


# 3d plots we're including in the paper

genes = c("CD70", "AXL", "CDH6") 
make3DPlotsWithoutSaving("Kidney Renal Clear Cell Carcinoma", genes, T) 
makeSplitViolins("Kidney Renal Clear Cell Carcinoma", genes, T)

# pancreatic top triples
genes = c("FAP", "PSCA", "MSLN") 
make3DPlotsWithoutSaving("Pancreatic Adenocarcinoma", genes, T) 
makeSplitViolins("Pancreatic Adenocarcinoma", genes, F)

genes = c("FAP", "CLDN18", "MSLN")
make3DPlotsWithoutSaving("Pancreatic Adenocarcinoma", genes, T) 
makeSplitViolins("Pancreatic Adenocarcinoma", genes, T)

genes = c("FAP", "CEACAM5", "MSLN") 
make3DPlotsWithoutSaving("Pancreatic Adenocarcinoma", genes, T) 
makeSplitViolins("Pancreatic Adenocarcinoma", genes, T)

genes = c("FAP", "KDR", "MSLN") 
make3DPlotsWithoutSaving("Pancreatic Adenocarcinoma", genes, T) 
makeSplitViolins("Pancreatic Adenocarcinoma", genes, F)

genes = c("FAP", "PTPRN", "MSLN")
make3DPlotsWithoutSaving("Pancreatic Adenocarcinoma", genes, T) 
makeSplitViolins("Pancreatic Adenocarcinoma", genes, T)


makeSplitViolins("Pancreatic Adenocarcinoma", c("FAP", "MSLN"))

genes = c("EGFR", "ERBB2", "MET")
make3DPlotsWithoutSaving("Lung Adenocarcinoma", genes, T) 


# make 3D html files for requested plotly files
request = fread("data/triples/greg-plot-request.txt", header=F)
apply(request, 1, function(row){
  genes = strsplit(row[["V1"]], ":")
  genes = unlist(genes)
  print(genes)
  make3DPlots(row[["V2"]], genes, F, T)
})


# -------------------------------
#          functions
# -------------------------------
makePlots <- function(target, genes, jitta=F)
{
  tmp = as.data.frame(rnaseqmat[,c(genes[1], genes[2], "tissue.cancer", "type"),with=F])
  tmp$target = ifelse(tmp$tissue.cancer == target, "tumor", "normal")
  tmp = tmp[tmp$target == "tumor" | (tmp$type == "tissue"),]
  
  tmp[,1] = ifelse(as.numeric(unlist(tmp[,1])) < -5, -5, as.numeric(unlist(tmp[,1])))
  tmp[,2] = ifelse(as.numeric(unlist(tmp[,2])) < -5, -5, as.numeric(unlist(tmp[,2])))
  
  if (jitta)
  {
    tmp[,1] = jitter(ifelse(as.numeric(unlist(tmp[,1])) < 0, 0, as.numeric(unlist(tmp[,1]))), 0.001)
    tmp[,2] = jitter(ifelse(as.numeric(unlist(tmp[,2])) < 0, 0, as.numeric(unlist(tmp[,2]))), 0.001)
  }
  
  tmp2 = as.data.table(tmp)
  names(tmp2) = c("g1", "g2", "tc", "type", "tar")
  tmpagg <- tmp2[, .(g1 = mean(g1), g2 = mean(g2)), by=tc]
  rm(tmp2)
  names(tmpagg) <- c("tc", names(tmp)[1], names(tmp)[2])
  tmpagg$target = ifelse(tmpagg$tc == target, "tumor", "normal")
  
  
  cplot = ggplot(tmp, aes_string(x=names(tmp)[1], y=names(tmp)[2])) +
    stat_density2d(aes(alpha=..level.., fill=target, color=target), size=0.25, bins=6, geom="polygon") +
    geom_point(data=tmpagg, aes_string(x=names(tmp)[1], y=names(tmp)[2], colour="target"), size=1.5, alpha=0.8) +
    geom_text_repel(data=tmpagg, aes(label=tc), size=3, box.padding = 0.25, 
                    point.padding = 0.2, segment.color = 'grey50', segment.size = 0.25) +
    scale_fill_manual(values=c("#34495e", "#e74c3c")) + 
    scale_color_manual(values=c("#34495e", "#c82a19", "#bdc3c7")) +
    scale_alpha_continuous(range = c(0.10, 0.8))  + theme_bw() +
    ylab(paste(names(tmp)[2], "log(TPM)")) + xlab(paste(names(tmp)[1], "log(TPM)")) + 
    theme(axis.title = element_text(size=16), axis.text= element_text(size=12), panel.grid.minor=element_blank(),
          plot.background=element_blank(), panel.background=element_blank(), panel.grid.major=element_blank()) +
    guides(colour=F)
  
  ggsave(paste0("scripts/clustering-method/res-figures/paper/", genes[1], "-", genes[2], "-", target, "-contour.pdf"), 
         cplot, width=6, height=4, device="pdf")
  
  
  tmp[,1] = ifelse(as.numeric(unlist(tmp[,1])) < 0, 0, as.numeric(unlist(tmp[,1])))
  tmp[,2] = ifelse(as.numeric(unlist(tmp[,2])) < 0, 0, as.numeric(unlist(tmp[,2])))
  
  splot = ggplot(tmp, aes_string(x=names(tmp)[1], y=names(tmp)[2])) + geom_point(aes(color=target), alpha = 0.3)  +
    scale_color_manual(values=c("#34495e", "#e74c3c", "#bdc3c7")) + theme_bw() + 
    ylab(paste(names(tmp)[2], "log(TPM)")) + xlab(paste(names(tmp)[1], "log(TPM)")) + 
    theme(legend.position="bottom", axis.title = element_text(size=16), 
          axis.text= element_text(size=12), panel.grid.minor=element_blank(),
          plot.background=element_blank(), panel.background=element_blank(), panel.grid.major=element_blank())
  
  ggsave(paste0("scripts/clustering-method/res-figures/paper/", genes[1], "-", genes[2], "-", target, "-scatter.pdf"), 
         splot, width=5, height=5, device="pdf")
}


makePlotsWithoutSaving <- function(target, genes, jitta=F)
{
  tmp = as.data.frame(rnaseqmat[,c(genes[1], genes[2], "tissue.cancer", "type"),with=F])
  tmp$target = ifelse(tmp$tissue.cancer == target, "tumor", "normal")
  tmp = tmp[tmp$target == "tumor" | (tmp$type == "tissue"),]
  
  tmp[,1] = ifelse(as.numeric(unlist(tmp[,1])) < -5, -5, as.numeric(unlist(tmp[,1])))
  tmp[,2] = ifelse(as.numeric(unlist(tmp[,2])) < -5, -5, as.numeric(unlist(tmp[,2])))
  
  if (jitta)
  {
    tmp[,1] = jitter(ifelse(as.numeric(unlist(tmp[,1])) < 0, 0, as.numeric(unlist(tmp[,1]))), 0.001)
    tmp[,2] = jitter(ifelse(as.numeric(unlist(tmp[,2])) < 0, 0, as.numeric(unlist(tmp[,2]))), 0.001)
  }
    
  tmp2 = as.data.table(tmp)
  names(tmp2) = c("g1", "g2", "tc", "type", "tar")
  tmpagg <- tmp2[, .(g1 = mean(g1), g2 = mean(g2)), by=tc]
  rm(tmp2)
  names(tmpagg) <- c("tc", names(tmp)[1], names(tmp)[2])
  tmpagg$target = ifelse(tmpagg$tc == target, "tumor", "normal2")
  tmpagg = tmpagg[target != "tumor",]

  # cplot = ggplot(tmp, aes_string(x=names(tmp)[1], y=names(tmp)[2])) +
  #   stat_density2d(aes(alpha=..level.., fill=target, color=target), size=0.25, bins=6, geom="polygon") +
  #   geom_point(data=tmpagg, aes_string(x=names(tmp)[1], y=names(tmp)[2], colour="target"), size=1.5, alpha=0.8) +
  #   geom_text_repel(data=tmpagg, aes(label=tc), size=3, box.padding = 0.25, 
  #                   point.padding = 0.2, segment.color = 'grey50', segment.size = 0.25) +
  #   scale_fill_manual(values=c("#34495e", "#e74c3c")) + 
  #   scale_color_manual(values=c("#34495e", "#c82a19", "#bdc3c7")) +
  #   scale_alpha_continuous(range = c(0.10, 0.8))  + theme_bw() + 
  #   ylab(paste(names(tmp)[2], "log(TPM)")) + xlab(paste(names(tmp)[1], "log(TPM)")) + 
  #   theme(axis.title = element_text(size=18), axis.text= element_text(size=16), panel.grid.minor=element_blank(),
  #         plot.background=element_blank(), panel.background=element_blank(), panel.grid.major=element_blank()) +
  #   guides(colour=F)
  # 
  # print(cplot)
  
  tmp[,1] = ifelse(as.numeric(unlist(tmp[,1])) < 0, 0, as.numeric(unlist(tmp[,1])))
  tmp[,2] = ifelse(as.numeric(unlist(tmp[,2])) < 0, 0, as.numeric(unlist(tmp[,2])))
  
  m = max(tmp[,1], tmp[,2])
  # splot = ggplot(tmp, aes_string(x=names(tmp)[1], y=names(tmp)[2])) + geom_point(aes(color=target), alpha = 0.3)  +
  #   scale_color_manual(values=c("#34495e", "#e74c3c", "#bdc3c7")) + theme_bw() + 
  #   ylab(paste(names(tmp)[2], "log(TPM)")) + xlab(paste(names(tmp)[1], "log(TPM)")) + 
  #   theme(legend.position="bottom", axis.title = element_text(size=18), 
  #         axis.text= element_text(size=16), panel.grid.minor=element_blank(),
  #         plot.background=element_blank(), panel.background=element_blank(), panel.grid.major=element_blank())
  # 
  # print(splot)
  
  splot = ggplot(tmp, aes_string(x=names(tmp)[1], y=names(tmp)[2])) +
    geom_point(aes(color=target), alpha = 0.7) +
    scale_color_manual(values=c("#eaecee", "#34495e", "#e74c3c", "#bdc3c7")) + 
    geom_point(data=tmpagg, aes_string(x=names(tmp)[1], y=names(tmp)[2], colour="target"), size=2.5, alpha=1) +
    geom_text_repel(data=tmpagg, aes(label=tc), size=4, box.padding = 0.25, 
                    point.padding = 0.2, segment.color = 'grey50', segment.size = 0.25) +
    scale_fill_manual(values=c("#34495e", "#e74c3c")) + 
    #scale_color_manual(values=c("#34495e", "#c82a19", "#bdc3c7")) + 
    ylim(0, m) + xlim(0, m) +
    ylab(paste(names(tmp)[2], "log(TPM)")) + xlab(paste(names(tmp)[1], "log(TPM)")) + 
    theme(axis.title = element_text(size=18), axis.text= element_text(size=16), panel.grid.minor=element_blank(),
          plot.background=element_blank(), panel.background=element_blank(), panel.grid.major=element_blank()) +
    guides(colour=F)
  
  print(splot) # scatter plot with all centroids and cancer
  
  # tmp = tmp[tmp$target == "tumor",]
  # 
  # splot = ggplot(tmp, aes_string(x=names(tmp)[1], y=names(tmp)[2])) +
  #   geom_point(aes(color=target), alpha = 0.2) +
  #   scale_color_manual(values=c("#34495e", "#e74c3c", "#bdc3c7")) + 
  #   geom_point(data=tmpagg, aes_string(x=names(tmp)[1], y=names(tmp)[2], colour="target"), size=2, alpha=1) +
  #   geom_text_repel(data=tmpagg, aes(label=tc), size=4, box.padding = 0.25, 
  #                   point.padding = 0.2, segment.color = 'grey50', segment.size = 0.25) +
  #   scale_fill_manual(values=c("#34495e", "#e74c3c")) + ylim(0, m) + xlim(0, m) +
  #   #scale_color_manual(values=c("#34495e", "#c82a19", "#bdc3c7")) + theme_bw() + 
  #   ylab(paste(names(tmp)[2], "log(TPM)")) + xlab(paste(names(tmp)[1], "log(TPM)")) + 
  #   theme(axis.title = element_text(size=18), axis.text= element_text(size=16), panel.grid.minor=element_blank(),
  #         plot.background=element_blank(), panel.background=element_blank(), panel.grid.major=element_blank()) +
  #   guides(colour=F)
  # 
  # print(splot) # scatter plot with all centroids shown over faded selected only
  
}


mSP <- function(target, genes, printToFile=F, jitta=T)
{
  tmp = as.data.frame(rnaseqmat[,c(genes[1], genes[2], "tissue.cancer", "type"),with=F])
  tmp$target = ifelse(tmp$tissue.cancer == target, "tumor", "normal")
  tmp = tmp[tmp$target == "tumor" | (tmp$type == "tissue"),]
  
  tmp[,1] = ifelse(as.numeric(unlist(tmp[,1])) < -5, -5, as.numeric(unlist(tmp[,1])))
  tmp[,2] = ifelse(as.numeric(unlist(tmp[,2])) < -5, -5, as.numeric(unlist(tmp[,2])))
  
  if (jitta)
  {
    tmp[,1] = jitter(ifelse(as.numeric(unlist(tmp[,1])) < 0, 0, as.numeric(unlist(tmp[,1]))), 0.001)
    tmp[,2] = jitter(ifelse(as.numeric(unlist(tmp[,2])) < 0, 0, as.numeric(unlist(tmp[,2]))), 0.001)
  }
  
  tmp2 = as.data.table(tmp)
  names(tmp2) = c("g1", "g2", "tc", "type", "tar")
  tmpagg <- tmp2[, .(g1 = mean(g1), g2 = mean(g2)), by=tc]
  rm(tmp2)
  names(tmpagg) <- c("tc", names(tmp)[1], names(tmp)[2])
  tmpagg$target = ifelse(tmpagg$tc == target, "tumor", "normal2")
  tmpagg = tmpagg[target != "tumor",]

  tmp[,1] = ifelse(as.numeric(unlist(tmp[,1])) < 0, 0, as.numeric(unlist(tmp[,1])))
  tmp[,2] = ifelse(as.numeric(unlist(tmp[,2])) < 0, 0, as.numeric(unlist(tmp[,2])))
  
  m = max(tmp[,1], tmp[,2])

  splot = ggplot(tmp, aes_string(x=names(tmp)[1], y=names(tmp)[2])) +
    geom_point(aes(color=target), alpha = 0.7) +
    scale_color_manual(values=c("#c2c8ce", "#34495e", "#e74c3c", "#bdc3c7")) + 
    geom_point(data=tmpagg, aes_string(x=names(tmp)[1], y=names(tmp)[2], colour="target"), size=2.5, alpha=1) +
    geom_text_repel(data=tmpagg, aes(label=tc), size=3, box.padding = 0.25, 
                    point.padding = 0.2, segment.color = 'grey50', segment.size = 0.25) +
    scale_fill_manual(values=c("#34495e", "#e74c3c")) + 
    ylim(0, m) + xlim(0, m) + theme_bw() + 
    ylab(paste(names(tmp)[2], "log(TPM)")) + xlab(paste(names(tmp)[1], "log(TPM)")) + 
    theme(axis.title = element_text(size=18), axis.text= element_text(size=16), panel.grid.minor=element_blank(),
          plot.background=element_blank(), panel.background=element_blank(), panel.grid.major=element_blank()) +
    guides(colour=F)
  
  if (printToFile)
  {
    ggsave(paste0("scripts/clustering-method/res-figures/paper/", genes[1], "-", genes[2], "-", target, "-update-scatter.pdf"), 
         splot, width=3.5, height=3.5, units="in", device="pdf", useDingbats=FALSE)
  }
  else
  {
    print(splot)
  }
}

mSP.single <- function(target, gene)
{
  tmp = as.data.frame(rnaseqmat[,c(gene, "tissue.cancer", "type"),with=F])
  tmp$target = ifelse(tmp$tissue.cancer == target, "tumor", "normal")
  tmp = tmp[tmp$target == "tumor" | (tmp$type == "tissue"),]
  
  tmp[,1] = ifelse(as.numeric(unlist(tmp[,1])) < 0, 0, as.numeric(unlist(tmp[,1])))
  
  tmp2 = as.data.table(tmp)
  names(tmp2) = c("g1", "tc", "type", "tar")
  tmpagg <- tmp2[, .(g1 = mean(g1)), by=tc]
  rm(tmp2)
  names(tmpagg) <- c("tc", names(tmp)[1])
  tmpagg$target = ifelse(tmpagg$tc == target, "tumor", "normal2")
  tmpagg = tmpagg[target != "tumor",]
  
  m = max(tmp[,1])
  
  ggplot(tmp, aes_string(x=names(tmp)[1], y="1")) +
    geom_point(aes(color=target), alpha = 0.7) +
    scale_color_manual(values=c("#c2c8ce", "#34495e", "#e74c3c", "#bdc3c7")) + 
    geom_point(data=tmpagg, aes_string(x=names(tmp)[1], y="1", colour="target"), size=2.5, alpha=1) +
    geom_text_repel(data=tmpagg, aes(label=tc), size=3, box.padding = 0.25, 
                    point.padding = 0.2, segment.color = 'grey50', segment.size = 0.25) +
    scale_fill_manual(values=c("#34495e", "#e74c3c")) + 
    ylim(0, m) + xlim(0, m) + theme_bw() + 
    ylab(paste(names(tmp)[2], "log(TPM)")) + xlab(paste(names(tmp)[1], "log(TPM)")) + 
    theme(axis.title = element_text(size=18), axis.text= element_text(size=16), panel.grid.minor=element_blank(),
          plot.background=element_blank(), panel.background=element_blank(), panel.grid.major=element_blank()) +
    guides(colour=F)
  
}



mSP.withCutoff <- function(target, genes, cutoff, printToFile=F, jitta=T)
{
  tmp = as.data.frame(rnaseqmat[,c(genes[1], genes[2], "tissue.cancer", "type"),with=F])
  tmp$target = ifelse(tmp$tissue.cancer == target, "tumor", "normal")
  tmp = tmp[tmp$target == "tumor" | (tmp$type == "tissue"),]
  
  tmp[,1] = ifelse(as.numeric(unlist(tmp[,1])) < -5, -5, as.numeric(unlist(tmp[,1])))
  tmp[,2] = ifelse(as.numeric(unlist(tmp[,2])) < -5, -5, as.numeric(unlist(tmp[,2])))
  
  if (jitta)
  {
    tmp[,1] = jitter(ifelse(as.numeric(unlist(tmp[,1])) < 0, 0, as.numeric(unlist(tmp[,1]))), 0.001)
    tmp[,2] = jitter(ifelse(as.numeric(unlist(tmp[,2])) < 0, 0, as.numeric(unlist(tmp[,2]))), 0.001)
  }
  
  tmp2 = as.data.table(tmp)
  names(tmp2) = c("g1", "g2", "tc", "type", "tar")
  tmpagg <- tmp2[, .(g1 = mean(g1), g2 = mean(g2)), by=tc]
  rm(tmp2)
  names(tmpagg) <- c("tc", names(tmp)[1], names(tmp)[2])
  tmpagg$target = ifelse(tmpagg$tc == target, "tumor", "normal2")
  tmpagg = tmpagg[target != "tumor",]
  
  tmp[,1] = ifelse(as.numeric(unlist(tmp[,1])) < 0, 0, as.numeric(unlist(tmp[,1])))
  tmp[,2] = ifelse(as.numeric(unlist(tmp[,2])) < 0, 0, as.numeric(unlist(tmp[,2])))
  
  m = max(tmp[,1], tmp[,2])
  
  splot = ggplot(tmp, aes_string(x=names(tmp)[1], y=names(tmp)[2])) +
    geom_point(aes(color=target), alpha = 0.7) +
    scale_color_manual(values=c("#c2c8ce", "#34495e", "#e74c3c", "#bdc3c7")) + 
    geom_point(data=tmpagg, aes_string(x=names(tmp)[1], y=names(tmp)[2], colour="target"), size=2.5, alpha=1) +
    geom_text_repel(data=tmpagg, aes(label=tc), size=3, box.padding = 0.25, 
                    point.padding = 0.2, segment.color = 'grey50', segment.size = 0.25) +
    scale_fill_manual(values=c("#34495e", "#e74c3c")) + 
    ylim(0, cutoff) + xlim(0, cutoff) + theme_bw() + 
    ylab(paste(names(tmp)[2], "log(TPM)")) + xlab(paste(names(tmp)[1], "log(TPM)")) + 
    theme(axis.title = element_text(size=18), axis.text= element_text(size=16), panel.grid.minor=element_blank(),
          plot.background=element_blank(), panel.background=element_blank(), panel.grid.major=element_blank()) +
    guides(colour=F)
  
  if (printToFile)
  {
    ggsave(paste0("scripts/clustering-method/res-figures/paper/", genes[1], "-", genes[2], "-", target, "-update-scatter.pdf"), 
           splot, width=3.5, height=3.5, units="in", device="pdf", useDingbats=FALSE)
  }
  else
  {
    print(splot)
  }
}

mSP.highlight <- function(target, genes, jitta=F, hilite=c())
{
  tmp = as.data.frame(rnaseqmat[,c(genes[1], genes[2], "tissue.cancer", "type"),with=F])
  tmp$target = ifelse(tmp$tissue.cancer == target, "tumor", "normal")
  tmp = tmp[tmp$target == "tumor" | (tmp$type == "tissue"),]
  
  if (jitta)
  {
    tmp[,1] = jitter(ifelse(as.numeric(unlist(tmp[,1])) < 0, 0, as.numeric(unlist(tmp[,1]))), 0.001)
    tmp[,2] = jitter(ifelse(as.numeric(unlist(tmp[,2])) < 0, 0, as.numeric(unlist(tmp[,2]))), 0.001)
  }
  else
  {
    tmp[,1] = ifelse(as.numeric(unlist(tmp[,1])) < 0, 0, as.numeric(unlist(tmp[,1])))
    tmp[,2] = ifelse(as.numeric(unlist(tmp[,2])) < 0, 0, as.numeric(unlist(tmp[,2])))
    
  }
  
  tmp2 = as.data.table(tmp)
  names(tmp2) = c("g1", "g2", "tc", "type", "tar")
  tmpagg <- tmp2[, .(g1 = mean(g1), g2 = mean(g2)), by=tc]
  rm(tmp2)
  names(tmpagg) <- c("tc", names(tmp)[1], names(tmp)[2])
  tmpagg$target = ifelse(tmpagg$tc == target, "tumor", "normal2")
  tmpagg = tmpagg[target != "tumor",]
  
  tmp[,1] = ifelse(as.numeric(unlist(tmp[,1])) < 0, 0, as.numeric(unlist(tmp[,1])))
  tmp[,2] = ifelse(as.numeric(unlist(tmp[,2])) < 0, 0, as.numeric(unlist(tmp[,2])))
  
  m = max(tmp[,1], tmp[,2])
  
  tmp$target = ifelse(tmp$tissue.cancer %in% hilite, "hilite", tmp$target)
  tmpagg$target = ifelse(tmpagg$tc %in% hilite, "hilite", tmpagg$target)
  
  #print(tmp[tmp$target == "hilite",])
  
  splot = ggplot(tmp[tmp$target != "hilite",], aes_string(x=names(tmp)[1], y=names(tmp)[2])) +
    geom_point(aes(color=target), alpha = 0.7) +
    scale_color_manual(values=c("#16a085", "#c2c8ce", "#34495e", "#e74c3c", "#bdc3c7")) + 
    geom_point(data=tmp[tmp$target=="hilite",], aes_string(x=names(tmp)[1], y=names(tmp)[2], colour="target"), alpha=0.7) +
    geom_point(data=tmpagg, aes_string(x=names(tmp)[1], y=names(tmp)[2], colour="target"), size=2.5, alpha=1) +
    # geom_text_repel(data=tmpagg, aes(label=tc), size=3, box.padding = 0.25, 
    #                 point.padding = 0.2, segment.color = 'grey50', segment.size = 0.25) +
    scale_fill_manual(values=c("#34495e", "#e74c3c")) + 
    ylim(0, m) + xlim(0, m) + theme_bw() + 
    ylab(paste(names(tmp)[2], "log(TPM)")) + xlab(paste(names(tmp)[1], "log(TPM)")) + 
    theme(axis.title = element_text(size=18), axis.text= element_text(size=16), panel.grid.minor=element_blank(),
          plot.background=element_blank(), panel.background=element_blank(), panel.grid.major=element_blank()) +
    guides(colour=F)
  
#  print(splot)
  ggsave(paste0("scripts/clustering-method/res-figures/paper/", genes[1], "-", genes[2], "-", target, "-update-hilite-scatter.pdf"), 
         splot, width=3.5, height=3.5, units="in", device="pdf", useDingbats=FALSE)
  
}  

mSP.highlight.withCutoff <- function(target, genes, cutoff, jitta=F, hilite=c())
{
  tmp = as.data.frame(rnaseqmat[,c(genes[1], genes[2], "tissue.cancer", "type"),with=F])
  tmp$target = ifelse(tmp$tissue.cancer == target, "tumor", "normal")
  tmp = tmp[tmp$target == "tumor" | (tmp$type == "tissue"),]
  
  if (jitta)
  {
    tmp[,1] = jitter(ifelse(as.numeric(unlist(tmp[,1])) < 0, 0, as.numeric(unlist(tmp[,1]))), 0.001)
    tmp[,2] = jitter(ifelse(as.numeric(unlist(tmp[,2])) < 0, 0, as.numeric(unlist(tmp[,2]))), 0.001)
  }
  else
  {
    tmp[,1] = ifelse(as.numeric(unlist(tmp[,1])) < 0, 0, as.numeric(unlist(tmp[,1])))
    tmp[,2] = ifelse(as.numeric(unlist(tmp[,2])) < 0, 0, as.numeric(unlist(tmp[,2])))
    
  }
  
  tmp2 = as.data.table(tmp)
  names(tmp2) = c("g1", "g2", "tc", "type", "tar")
  tmpagg <- tmp2[, .(g1 = mean(g1), g2 = mean(g2)), by=tc]
  rm(tmp2)
  names(tmpagg) <- c("tc", names(tmp)[1], names(tmp)[2])
  tmpagg$target = ifelse(tmpagg$tc == target, "tumor", "normal2")
  tmpagg = tmpagg[target != "tumor",]
  
  tmp[,1] = ifelse(as.numeric(unlist(tmp[,1])) < 0, 0, as.numeric(unlist(tmp[,1])))
  tmp[,2] = ifelse(as.numeric(unlist(tmp[,2])) < 0, 0, as.numeric(unlist(tmp[,2])))
  
  m = max(tmp[,1], tmp[,2])
  
  tmp$target = ifelse(tmp$tissue.cancer %in% hilite, "hilite", tmp$target)
  tmpagg$target = ifelse(tmpagg$tc %in% hilite, "hilite", tmpagg$target)
  
  #print(tmp[tmp$target == "hilite",])
  
  splot = ggplot(tmp[tmp$target != "hilite",], aes_string(x=names(tmp)[1], y=names(tmp)[2])) +
    geom_point(aes(color=target), alpha = 0.7) +
    scale_color_manual(values=c("#16a085", "#c2c8ce", "#34495e", "#e74c3c", "#bdc3c7")) + 
    geom_point(data=tmp[tmp$target=="hilite",], aes_string(x=names(tmp)[1], y=names(tmp)[2], colour="target"), alpha=0.7) +
    geom_point(data=tmpagg, aes_string(x=names(tmp)[1], y=names(tmp)[2], colour="target"), size=2.5, alpha=1) +
    # geom_text_repel(data=tmpagg, aes(label=tc), size=3, box.padding = 0.25, 
    #                 point.padding = 0.2, segment.color = 'grey50', segment.size = 0.25) +
    scale_fill_manual(values=c("#34495e", "#e74c3c")) + 
    ylim(0, cutoff) + xlim(0, cutoff) + theme_bw() + 
    ylab(paste(names(tmp)[2], "log(TPM)")) + xlab(paste(names(tmp)[1], "log(TPM)")) + 
    theme(axis.title = element_text(size=18), axis.text= element_text(size=16), panel.grid.minor=element_blank(),
          plot.background=element_blank(), panel.background=element_blank(), panel.grid.major=element_blank()) +
    guides(colour=F)
  
  #  print(splot)
  ggsave(paste0("scripts/clustering-method/res-figures/paper/", genes[1], "-", genes[2], "-", target, "-update-hilite-scatter.pdf"), 
         splot, width=3.5, height=3.5, units="in", device="pdf", useDingbats=FALSE)
  
} 

make3DPlots <- function(target, genes, jitta=F, save=F)
{
  tmp = as.data.frame(rnaseqmat[,c(genes[1], genes[2], genes[3], "tissue.cancer", "type"),with=F])
  tmp$target = ifelse(tmp$tissue.cancer == target, "tumor", "normal")
  tmp = tmp[tmp$target == "tumor" | (tmp$type == "tissue"),]
  
  tmp[,1] = ifelse(as.numeric(unlist(tmp[,1])) < 0, 0, as.numeric(unlist(tmp[,1])))
  tmp[,2] = ifelse(as.numeric(unlist(tmp[,2])) < 0, 0, as.numeric(unlist(tmp[,2])))
  tmp[,3] = ifelse(as.numeric(unlist(tmp[,3])) < 0, 0, as.numeric(unlist(tmp[,3])))
  
  if (jitta)
  {
    tmp[,1] = jitter(ifelse(as.numeric(unlist(tmp[,1])) < 0, 0, as.numeric(unlist(tmp[,1]))), 0.001)
    tmp[,2] = jitter(ifelse(as.numeric(unlist(tmp[,2])) < 0, 0, as.numeric(unlist(tmp[,2]))), 0.001)
    tmp[,3] = jitter(ifelse(as.numeric(unlist(tmp[,3])) < 0, 0, as.numeric(unlist(tmp[,3]))), 0.001)
  }
  
  tmp2 = as.data.table(tmp)
  names(tmp2) = c("g1", "g2", "g3", "tc", "type", "tar")
  tmpagg <- tmp2[, .(g1 = mean(g1), g2 = mean(g2), g3 = mean(g3)), by=tc]
  rm(tmp2)
  tmpagg$tar = ifelse(tmpagg$tc == target, "tumor", "normal.centroid")
  tmpagg = tmpagg[tar != "tumor",]
  tmpagg[, sz := 3.5]
  tmpagg[, opa := 1]
  
  names(tmp) = c("g1", "g2", "g3", "tc", "type", "tar")
  tmp = as.data.table(tmp)
  tmp[, type := NULL]
  tmp[, sz := 3]
  tmp[, opa := ifelse(target == "normal", 0.2, 0.7)]
  #print(head(tmp))
  #print(head(tmpagg))
  
  tmp = rbind(tmp, tmpagg)
  
  p <- plot_ly(tmp, x=~g1, y=~g2, z=~g3, type="scatter3d", mode="markers", 
          color=~tar, colors=c("#c2c8ce","#34495e", "#e74c3c"), 
          marker=list(size=5, opacity=0.5)) %>%
  layout(
    title = target,
    scene = list(
      xaxis = list(title = genes[1], tickfont = list(size=16), titlefont=list(size=24)),
      yaxis = list(title = genes[2], tickfont = list(size=16), titlefont=list(size=24)),
      zaxis = list(title = genes[3], tickfont = list(size=16), titlefont=list(size=24))
    ))
  
  if(save)
  {
    htmlwidgets::saveWidget(p, paste0(getwd(),"/scripts/clustering-method/res-figures/html/",gsub(" ", "-", target), "-",genes[1], "-", genes[2], "-", genes[3],".html"))  
  }
  else
  {
    print(p)
  }
  
}

makePlotsWithoutSavingHighlights <- function(target, genes, jitta=F, hilite=c())
{
  tmp = as.data.frame(rnaseqmat[,c(genes[1], genes[2], "tissue.cancer", "type"),with=F])
  tmp$target = ifelse(tmp$tissue.cancer == target, "tumor", "normal")
  tmp = tmp[tmp$target == "tumor" | (tmp$type == "tissue"),]
  
  if (jitta)
  {
    tmp[,1] = jitter(ifelse(as.numeric(unlist(tmp[,1])) < 0, 0, as.numeric(unlist(tmp[,1]))), 0.001)
    tmp[,2] = jitter(ifelse(as.numeric(unlist(tmp[,2])) < 0, 0, as.numeric(unlist(tmp[,2]))), 0.001)
  }
  else
  {
    tmp[,1] = ifelse(as.numeric(unlist(tmp[,1])) < 0, 0, as.numeric(unlist(tmp[,1])))
    tmp[,2] = ifelse(as.numeric(unlist(tmp[,2])) < 0, 0, as.numeric(unlist(tmp[,2])))
    
  }
  
  tmp2 = as.data.table(tmp)
  names(tmp2) = c("g1", "g2", "tc", "type", "tar")
  tmpagg <- tmp2[, .(g1 = mean(g1), g2 = mean(g2)), by=tc]
  rm(tmp2)
  names(tmpagg) <- c("tc", names(tmp)[1], names(tmp)[2])
  tmpagg$target = ifelse(tmpagg$tc == target, "tumor", "normal")
  
  
  cplot = ggplot(tmp, aes_string(x=names(tmp)[1], y=names(tmp)[2])) +
    geom_point(aes(color=target), alpha = 0.10) +
    scale_color_manual(values=c("#34495e", "#ff3a24", "#bdc3c7")) + 
    geom_point(data=tmpagg, aes_string(x=names(tmp)[1], y=names(tmp)[2], colour="target"), size=2, alpha=1) +
    geom_text_repel(data=tmpagg, aes(label=tc), size=4, box.padding = 0.25, 
                    point.padding = 0.2, segment.color = 'grey50', segment.size = 0.25) +
    scale_fill_manual(values=c("#34495e", "#e74c3c")) + 
    #scale_color_manual(values=c("#34495e", "#c82a19", "#bdc3c7")) + theme_bw() + 
    ylab(paste(names(tmp)[2], "log(TPM)")) + xlab(paste(names(tmp)[1], "log(TPM)")) + coord_equal(ratio=1) +
    theme(axis.title = element_text(size=18), axis.text= element_text(size=16), panel.grid.minor=element_blank(),
          plot.background=element_blank(), panel.background=element_blank(), panel.grid.major=element_blank()) +
    guides(colour=F)
  
  print(cplot) # scatter plot with all centroids shown over faded samples
  

  tmp$target = ifelse(tmp$tissue.cancer %in% hilite, "hilite", tmp$target)
  tmpagg$target = ifelse(tmpagg$tc %in% hilite, "hilite", tmpagg$target)
  
  splot = ggplot(tmp, aes_string(x=names(tmp)[1], y=names(tmp)[2])) +
    geom_point(aes(color=target), alpha = 0.10) +
    scale_color_manual(values=c("#16a085", "#34495e", "#ff3a24", "#8e44ad")) + 
    geom_point(data=tmpagg, aes_string(x=names(tmp)[1], y=names(tmp)[2], colour="target"), size=2, alpha=1) +
    geom_text_repel(data=tmpagg, aes(label=tc), size=4, box.padding = 0.25, 
                    point.padding = 0.2, segment.color = 'grey50', segment.size = 0.25) +
    scale_fill_manual(values=c("#34495e", "#e74c3c")) + 
    #scale_color_manual(values=c("#34495e", "#c82a19", "#bdc3c7")) + theme_bw() + 
    ylab(paste(names(tmp)[2], "log(TPM)")) + xlab(paste(names(tmp)[1], "log(TPM)")) + coord_equal(ratio=1) +
    theme(axis.title = element_text(size=18), axis.text= element_text(size=16), panel.grid.minor=element_blank(),
          plot.background=element_blank(), panel.background=element_blank(), panel.grid.major=element_blank()) +
    guides(colour=F)
  
  print(splot) # scatter plot with all centroids and highlighted tissue
  
  tmp = tmp[tmp$target == "tumor" | (tmp$target == "hilite"),]
  
  splot = ggplot(tmp, aes_string(x=names(tmp)[1], y=names(tmp)[2])) +
    geom_point(aes(color=target), alpha = 0.2) +
    scale_color_manual(values=c("#16a085", "#34495e", "#ff3a24", "#8e44ad")) + 
    geom_point(data=tmpagg, aes_string(x=names(tmp)[1], y=names(tmp)[2], colour="target"), size=2, alpha=1) +
    geom_text_repel(data=tmpagg, aes(label=tc), size=4, box.padding = 0.25, 
                    point.padding = 0.2, segment.color = 'grey50', segment.size = 0.25) +
    scale_fill_manual(values=c("#34495e", "#e74c3c")) + 
    #scale_color_manual(values=c("#34495e", "#c82a19", "#bdc3c7")) + theme_bw() + 
    ylab(paste(names(tmp)[2], "log(TPM)")) + xlab(paste(names(tmp)[1], "log(TPM)")) + coord_equal(ratio=1) +
    theme(axis.title = element_text(size=18), axis.text= element_text(size=16), panel.grid.minor=element_blank(),
          plot.background=element_blank(), panel.background=element_blank(), panel.grid.major=element_blank()) +
    guides(colour=F)
  
  print(splot) # scatter plot with all centroids shown over faded selected only
  
}

make3DPlotsBLAH <- function(target, genes, jitta=F)
{
  tmp = as.data.frame(rnaseqmat[,c(genes[1], genes[2], genes[3], "tissue.cancer", "type"),with=F])
  tmp$target = ifelse(tmp$tissue.cancer == target, "tumor", "normal")
  tmp = tmp[tmp$target == "tumor" | (tmp$type == "tissue"),]
  
  tmp[,1] = ifelse(as.numeric(unlist(tmp[,1])) < -5, -5, as.numeric(unlist(tmp[,1])))
  tmp[,2] = ifelse(as.numeric(unlist(tmp[,2])) < -5, -5, as.numeric(unlist(tmp[,2])))
  tmp[,3] = ifelse(as.numeric(unlist(tmp[,3])) < -5, -5, as.numeric(unlist(tmp[,3])))
  
  if (jitta)
  {
    tmp[,1] = jitter(ifelse(as.numeric(unlist(tmp[,1])) < 0, 0, as.numeric(unlist(tmp[,1]))), 0.001)
    tmp[,2] = jitter(ifelse(as.numeric(unlist(tmp[,2])) < 0, 0, as.numeric(unlist(tmp[,2]))), 0.001)
    tmp[,3] = jitter(ifelse(as.numeric(unlist(tmp[,3])) < 0, 0, as.numeric(unlist(tmp[,3]))), 0.001)
  }
  
  names(tmp) = c("g1", "g2", "g3", "tc", "type", "tar")
  
  p = plot_ly(x=tmp$g1, y=tmp$g2, z=tmp$g3, type="scatter3d", mode="markers", color=tmp$tar, colors=c("#34495e", "#e74c3c"), size=5) %>%
    layout(
      title = target,
      scene = list(
        xaxis = list(title = genes[1]),
        yaxis = list(title = genes[2]),
        zaxis = list(title = genes[3])
      ))
  
  orca(p, paste0("scripts/clustering-method/res-figures/paper/", genes[1], "-", genes[2], "-", genes[3], "-", target, "-scatter.png"))
  
}

makeSplitViolins <- function(target, genes, printToFile=F)
{
  if (length(genes) > 2)
  {
    tmp = as.data.frame(rnaseqmat[,c(genes[1], genes[2], genes[3], "tissue.cancer", "type"),with=F])
    tmp$target = ifelse(tmp$tissue.cancer == target, "tumor", "normal")
    tmp = tmp[tmp$target == "tumor" | (tmp$type == "tissue"),]
    
    tmp[,1] = ifelse(as.numeric(unlist(tmp[,1])) < 0, 0, as.numeric(unlist(tmp[,1])))
    tmp[,2] = ifelse(as.numeric(unlist(tmp[,2])) < 0, 0, as.numeric(unlist(tmp[,2])))
    tmp[,3] = ifelse(as.numeric(unlist(tmp[,3])) < 0, 0, as.numeric(unlist(tmp[,3])))
    tmp = as.data.table(tmp)
    tmp[, ANDgate := pmin(tmp[[1]], tmp[[2]], tmp[[3]])]
  }
  else
  {
    tmp = as.data.frame(rnaseqmat[,c(genes[1], genes[2], "tissue.cancer", "type"),with=F])
    tmp$target = ifelse(tmp$tissue.cancer == target, "tumor", "normal")
    tmp = tmp[tmp$target == "tumor" | (tmp$type == "tissue"),]
    
    tmp[,1] = ifelse(as.numeric(unlist(tmp[,1])) < 0, 0, as.numeric(unlist(tmp[,1])))
    tmp[,2] = ifelse(as.numeric(unlist(tmp[,2])) < 0, 0, as.numeric(unlist(tmp[,2])))
    tmp = as.data.table(tmp)
    tmp[, ANDgate := pmin(tmp[[1]], tmp[[2]])]
    
  }
  ttmp = tmp
  ttmp[,sample := 1:.N]
  ttmp[, tissue.cancer := NULL]
  ttmp[, target := NULL]
  ttmp = melt(ttmp, id=c("sample", "type"))
  ttmp[, type := ifelse(type == "cancer" | type == "tumor", "tumor", "normal")]
  gg = ggplot(ttmp, aes(variable, value, fill = type)) + geom_split_violin() + xlab("") +
    scale_fill_manual(values=c("#34495e", "#e74c3c", "#bdc3c7")) + 
    ylab("log(TPM)") + coord_flip() + theme_bw() + theme(legend.position="none", axis.title = element_text(size=16), 
                                                         axis.text= element_text(size=16))
  
  
  if (printToFile)
  {
    if (length(genes) > 2)
    {
      ggsave(paste0("scripts/clustering-method/res-figures/paper/", genes[1], "-", genes[2], "-", genes[3], "-", target, "-splitViolin.pdf"), 
           gg, width=3.5, height=3.5, units="in", device="pdf", useDingbats=FALSE)
    }
    else
    {
      ggsave(paste0("scripts/clustering-method/res-figures/paper/", genes[1], "-", genes[2], "-", target, "-splitViolin.pdf"), 
             gg, width=3, height=9.5, units="in", device="pdf", useDingbats=FALSE)
      
    }
  }
  else
  {
    print(gg)
  }
  
}



makeSplitViolins.facet <- function(target, genes, printToFile=F)
{
  if (length(genes) > 2)
  {
    tmp = as.data.frame(rnaseqmat[,c(genes[1], genes[2], genes[3], "tissue.cancer", "type"),with=F])
    tmp$target = ifelse(tmp$tissue.cancer == target, "tumor", "normal")
    tmp = tmp[tmp$target == "tumor" | (tmp$type == "tissue"),]
    
    tmp[,1] = ifelse(as.numeric(unlist(tmp[,1])) < 0, 0, as.numeric(unlist(tmp[,1])))
    tmp[,2] = ifelse(as.numeric(unlist(tmp[,2])) < 0, 0, as.numeric(unlist(tmp[,2])))
    tmp[,3] = ifelse(as.numeric(unlist(tmp[,3])) < 0, 0, as.numeric(unlist(tmp[,3])))
    tmp = as.data.table(tmp)
    tmp[, ANDgate := pmin(tmp[[1]], tmp[[2]], tmp[[3]])]
  }
  else
  {
    tmp = as.data.frame(rnaseqmat[,c(genes[1], genes[2], "tissue.cancer", "type"),with=F])
    tmp$target = ifelse(tmp$tissue.cancer == target, "tumor", "normal")
    tmp = tmp[tmp$target == "tumor" | (tmp$type == "tissue"),]
    
    tmp[,1] = ifelse(as.numeric(unlist(tmp[,1])) < 0, 0, as.numeric(unlist(tmp[,1])))
    tmp[,2] = ifelse(as.numeric(unlist(tmp[,2])) < 0, 0, as.numeric(unlist(tmp[,2])))
    tmp = as.data.table(tmp)
    tmp[, ANDgate := pmin(tmp[[1]], tmp[[2]])]
    
  }
  ttmp = tmp
  ttmp[,sample := 1:.N]
  ttmp[, tissue.cancer := NULL]
  ttmp[, target := NULL]
  ttmp = melt(ttmp, id=c("sample", "type"))
  ttmp[, type := ifelse(type == "cancer" | type == "tumor", "tumor", "normal")]
  gg = ggplot(ttmp, aes("antigen", value, fill = type)) + geom_split_violin() + xlab("") +
    scale_fill_manual(values=c("#34495e", "#e74c3c", "#bdc3c7")) + facet_wrap(~variable, nrow=1) + 
    ylab("log(TPM)") + theme_bw() + theme(legend.position="none", axis.title = element_text(size=16), 
                                                         axis.text= element_text(size=16), axis.title.x=element_blank(),
                                                         axis.text.x=element_blank(), axis.ticks.x=element_blank())
  
  gg2 = ggplot(ttmp[variable == "ANDgate",], aes(x=type, y=log2(value + 1))) + geom_violin()
  gg3 = ggplot(ttmp[variable == "ANDgate",], aes(x="antigen", y=value)) + geom_violin() + facet_wrap(~type)
  
  if (printToFile)
  {
    if (length(genes) > 2)
    {
      ggsave(paste0("scripts/clustering-method/res-figures/paper/", genes[1], "-", genes[2], "-", genes[3], "-", target, "-splitViolin.pdf"), 
             gg, width=3.5, height=3.5, units="in", device="pdf", useDingbats=FALSE)
    }
    else
    {
      ggsave(paste0("scripts/clustering-method/res-figures/paper/", genes[1], "-", genes[2], "-", target, "-splitViolin.pdf"), 
             gg, width=3, height=9.5, units="in", device="pdf", useDingbats=FALSE)
      
    }
  }
  else
  {
    print(gg)
    print(gg2)
    print(gg3)
  }
  
}

makeBeeswarms <- function(target, genes, printToFile=F)
{
  if (length(genes) > 2)
  {
    tmp = as.data.frame(rnaseqmat[,c(genes[1], genes[2], genes[3], "tissue.cancer", "type"),with=F])
    tmp$target = ifelse(tmp$tissue.cancer == target, "tumor", "normal")
    tmp = tmp[tmp$target == "tumor" | (tmp$type == "tissue"),]
    
    tmp[,1] = ifelse(as.numeric(unlist(tmp[,1])) < 0, 0, as.numeric(unlist(tmp[,1])))
    tmp[,2] = ifelse(as.numeric(unlist(tmp[,2])) < 0, 0, as.numeric(unlist(tmp[,2])))
    tmp[,3] = ifelse(as.numeric(unlist(tmp[,3])) < 0, 0, as.numeric(unlist(tmp[,3])))
    tmp = as.data.table(tmp)
    tmp[, ANDgate := pmin(tmp[[1]], tmp[[2]], tmp[[3]])]
  }
  else
  {
    tmp = as.data.frame(rnaseqmat[,c(genes[1], genes[2], "tissue.cancer", "type"),with=F])
    tmp$target = ifelse(tmp$tissue.cancer == target, "tumor", "normal")
    tmp = tmp[tmp$target == "tumor" | (tmp$type == "tissue"),]
    
    tmp[,1] = ifelse(as.numeric(unlist(tmp[,1])) < 0, 0, as.numeric(unlist(tmp[,1])))
    tmp[,2] = ifelse(as.numeric(unlist(tmp[,2])) < 0, 0, as.numeric(unlist(tmp[,2])))
    tmp = as.data.table(tmp)
    tmp[, ANDgate := pmin(tmp[[1]], tmp[[2]])]
    
  }
  ttmp = tmp
  ttmp[,sample := 1:.N]
  ttmp[, tissue.cancer := NULL]
  ttmp[, target := NULL]
  ttmp = melt(ttmp, id=c("sample", "type"))
  ttmp[, type := ifelse(type == "cancer" | type == "tumor", "tumor", "normal")]
  
  p <- ggplot(ttmp, aes(variable,value,color=factor(type), group=variable)) + 
    geom_quasirandom(method = 'tukeyDense', alpha=0.5) + scale_color_manual(values=c("#34495e", "#e74c3c", "#bdc3c7"))  +
    ylab("log(TPM)") + xlab("") + theme_bw() + theme(legend.position="none", axis.title = element_text(size=16), 
                                                   axis.text= element_text(size=16))
  gg <- sb(p, "#34495e")
  
  if (printToFile)
  {
    print("no save option yet")
  }
  else
  {
    plot(gg)
  }
  
}

plotDensity <- function(target, gene)
{
  # tmp = as.data.frame(rnaseqmat[,c(gene, "tissue.cancer", "type"),with=F])
  # tmp$target = ifelse(tmp$tissue.cancer == target, "tumor", "normal")
  # tmp = tmp[tmp$target == "tumor" | (tmp$type == "tissue"),]

  # calculate it on the sketch instead of the whole dataset
  sketchL = sketchLAll[can==target]$dataset
  s2 = rnaseqmat$dataset
  s2 = s2[grepl("^TCGA.*", s2)] # keep all TCGA
  sketchL = append(sketchL, s2)
  tmp = rnaseqmat[dataset %in% sketchL,]
  tmp = as.data.frame(tmp[,c(gene, "tissue.cancer", "type"),with=F])
  tmp$target = ifelse(tmp$tissue.cancer == target, "tumor", "normal")
  tmp = tmp[tmp$target == "tumor" | (tmp$type == "tissue"),]
  
  tmp[,1] = ifelse(as.numeric(unlist(tmp[,1])) < 0, 0, as.numeric(unlist(tmp[,1])))
  
  ggplot(tmp, aes_string(x=names(tmp)[1], fill="target")) + geom_density(alpha=0.75) + 
    scale_color_manual(values=c("#c2c8ce", "#34495e", "#e74c3c", "#bdc3c7")) + 
    scale_fill_manual(values=c("#34495e", "#e74c3c")) +  theme_bw() + 
    ylab("") + xlab("") + 
    theme(axis.title = element_text(size=18), axis.text= element_text(size=16), panel.grid.minor=element_blank(),
          plot.background=element_blank(), panel.background=element_blank(), panel.grid.major=element_blank()) +
    guides(fill=F)
  
}

plotDensityAll <- function(gene)
{
  # tmp = as.data.frame(rnaseqmat[,c(gene, "tissue.cancer", "type"),with=F])
  # tmp$target = ifelse(tmp$tissue.cancer == target, "tumor", "normal")
  # tmp = tmp[tmp$target == "tumor" | (tmp$type == "tissue"),]
  
  # calculate it on the sketch instead of the whole dataset
  sketchL = sketchLAll[can=="Brain Lower Grade Glioma"]$dataset # picked a random cancer here to get subset of normals d
  s2 = rnaseqmat$dataset
  s2 = s2[grepl("^TCGA.*", s2)] # keep all TCGA
  sketchL = append(sketchL, s2)
  tmp = rnaseqmat[dataset %in% sketchL,]
  tmp = as.data.frame(tmp[,c(gene, "tissue.cancer", "type"),with=F])
  tmp[,1] = ifelse(as.numeric(unlist(tmp[,1])) < 0, 0, as.numeric(unlist(tmp[,1])))
  
  #m = max(tmp[,1])
  
  gg1 = ggplot(tmp, aes_string(x=names(tmp)[1], fill="type")) + geom_density(alpha=0.75) + 
    scale_color_manual(values=c("#c2c8ce", "#34495e", "#e74c3c", "#bdc3c7")) + 
    scale_fill_manual(values=c("#34495e", "#e74c3c")) +  theme_bw() + 
    ylab("density") + xlab(paste(names(tmp)[1], "log(TPM)")) +
    theme(axis.title = element_text(size=18), axis.text= element_text(size=16), panel.grid.minor=element_blank(),
          plot.background=element_blank(), panel.background=element_blank(), panel.grid.major=element_blank()) +
    guides(fill=F)
  
  # gg2 = ggplot(tmp[tmp$type=="tumor",], aes_string(x=names(tmp)[1], fill="type")) + geom_density(alpha=1) +
  #   scale_fill_manual(values=c("#e74c3c")) +  theme_bw() +
  #   ylab("density") + xlab(paste(names(tmp)[1], "log(TPM)")) + facet_wrap(~tissue.cancer) +
  # theme(axis.title = element_text(size=18), axis.text= element_text(size=16), panel.grid.minor=element_blank(),
  #       plot.background=element_blank(), panel.background=element_blank(), panel.grid.major=element_blank()) +
  #   guides(fill=F)
  gg3 = ggplot(tmp[tmp$type=="tumor",], aes_string(x=paste0("reorder(tissue.cancer, ",names(tmp)[1],", FUN=mean)"), y= names(tmp)[1])) + 
    geom_boxplot(alpha=1) + theme_bw() +
    xlab("") + ylab(paste(names(tmp)[1], "log(TPM)")) + scale_x_discrete(label=function(x) abbreviate(x, minlength=15)) + coord_flip() + 
    theme(axis.title = element_text(size=18), axis.text= element_text(size=16),
          panel.grid.minor=element_blank(),
          plot.background=element_blank(), panel.background=element_blank(), panel.grid.major=element_blank()) +
    guides(fill=F) + ggtitle(names(tmp)[1])
  

  
  print(gg1)
  print(gg3)
  
}


# playing around with idea of 3 split violin plots
# https://stackoverflow.com/questions/35717353/split-violin-plot-with-ggplot2
GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, 
                           draw_group = function(self, data, ..., draw_quantiles = NULL) {
                             data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
                             grp <- data[1, "group"]
                             newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
                             newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
                             newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
                             
                             if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                               stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
                                                                         1))
                               quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
                               aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
                               aesthetics$alpha <- rep(1, nrow(quantiles))
                               both <- cbind(quantiles, aesthetics)
                               quantile_grob <- GeomPath$draw_panel(both, ...)
                               ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
                             }
                             else {
                               ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
                             }
                           })

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, 
        position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}

# split beeswarm
sb <- function(plot, clr){
  p <- ggplot_build(plot)
  # update the layer data
  #ms <- p$layout$panel_params[[1]]$x.major_source
  p$data[[1]] <- p$data[[1]] %>%
    mutate(x=case_when(
      colour==clr ~ as.numeric(as.character(group)) + abs(as.numeric(as.character(group)) - x),
      TRUE ~ as.numeric(as.character(group)) - abs(as.numeric(as.character(group)) - x))
    )
  # plot the update
  return(ggplot_gtable(p))
}


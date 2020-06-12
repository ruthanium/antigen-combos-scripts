# **************************************
# Author: Ruth Dannenfelser
# Date:   December 4, 2019
#
# cleaning up the aggregate stats
# script - just makes the plots
# messy grunt work done in 
# DensitySketchedBasedAggregateStats
#
# **************************************
library(ggplot2)
library(data.table)
library(sp)
library(parallel)
library(e1071)

setwd("/Genomics/ogtr04/rd6/BioGPS")

set.seed(1066)

# get the markers set up
tm = fread("data/transmembrane-markers/filteredTransmembraneList-updated-12-19.txt", header=F)
rnaseqmat = readRDS("scripts/clustering-method/rdata/combat-rnaseq-tm-clin-backcompatible.RData")
cans = unique(rnaseqmat[rnaseqmat$type == "tumor",]$tissue.cancer)

clinSmall = fread("data/transmembrane-markers/clin-paper-dec-19.txt", header=F)
clinSmall = clinSmall$V1
clinSmall = clinSmall[clinSmall %in% names(rnaseqmat)]

tm = tm$V1
tm = tm[tm %in% names(rnaseqmat)]


# work with all the cancer types - the full set of results (for clinicals)
bigres = readRDS("scripts/clustering-method/rdata/all-rnaseq-combat-sketch-0.2-top-fixed-backcompatible.RData")

bigres.clean = bigres[(gene1 %in% tm | gene1 %in% clinSmall) & (gene2 %in% tm | gene2 %in% clinSmall),]
bigres = bigres.clean
rm(bigres.clean)

# load sketches
sketchLAll = rbindlist(lapply(cans,
                              function(can) data.table(can=can,
                                                       dataset=readLines(paste0("scripts/clustering-method/sketches/", gsub(" ", "-", tolower(can)), "-sketch.txt")))))

# load triples
#tripleres = readRDS("scripts/clustering-method/rdata/all-rnaseq-combat-sketch-0.2-clin-triples.RData")
topTriples = readRDS("scripts/clustering-method/rdata/topCCCTriples-10.RData")

# load singles
singles = readRDS("scripts/clustering-method/rdata/singles-complete.RData")

#singTmp = readRDS("scripts/clustering-method/rdata/singles.RData") 
#names(singTmp) = c("can", "gene", "db", "dist.man", "ttype")
#singles = merge(singles, singTmp, by=c("can", "gene"))
#rm(singTmp)
#singles[, glab := ifelse(gene %in% clinSmall, "C", "N") ] # label C or S
#saveRDS(singles, "scripts/clustering-method/rdata/singles-complete.RData")
singles[,ttype := glab] # remove all singles that are just low
singles[, glab := NULL]
singles = singles[!is.na(db),]# remove the one NA - seems like it was completely stacked on each other

#TO DO NEED TO order singles
singles[,dbscaled := 1-((log(db) - min(log(db)))/(max(log(db)) - min(log(db))))] # make score between 0 and 1
singles[,distscaled := (log(dist.man) - min(log(dist.man)))/(max(log(dist.man)) - min(log(dist.man)))]
singles[,score := apply(singles[,c("dbscaled", "distscaled"), with=F], 1, min)]


# need to reduce the singles to top k per tumor type
k = 50
singlesN = singles[ttype=="N",]
singlesN = singlesN[order(-score),]
singlesN = singlesN[, head(.SD, k), by=can]

topSingles = rbind(singlesN, singles[ttype == "C",])

maxClin = singles[singles[ttype=="C", .I[which.max(f1)], by=can]$V1]

# load doubles
doubles = readRDS("scripts/clustering-method/rdata/doubles-topPerClin.RData")

# plotssss
# reduce the doubles down to the top k - doubles are top 100 orginally
k = 10
doublesTmp = doubles[order(-score),]
doublesTmp = doublesTmp[, head(.SD, k), by=c("type", "ttype")]

df = topSingles
names(df) <- c("type", "pair", "f1", "prec", "rec", "uth", "acc", "db", "dist.man", "ttype", "dbscaled", "distscaled", "score")
df = df[,c("pair", "type", "f1", "prec", "rec", "ttype", "db", "dist.man", "score")]


df = rbind(df, doublesTmp[,c("pair", "type", "f1", "prec", "rec", "ttype", "db", "dist.man", "score")], fill=T)

ggplot(df, aes(x=ttype, y=f1)) + geom_boxplot(notch = T) + 
  theme_bw() + ylim(0,1) + xlab("") + ylab("F1 score") + 
  theme(legend.position="bottom", axis.title = element_text(size=16), 
        axis.text= element_text(size=16), panel.grid.minor=element_blank(),
        plot.background=element_blank(), panel.background=element_blank(), panel.grid.major=element_blank())

wilcox.test(f1 ~ ttype, data = df[ttype == "C" | ttype == "C:C",])
wilcox.test(f1 ~ ttype, data = df[ttype == "C" | ttype == "C:N",])
wilcox.test(f1 ~ ttype, data = df[ttype == "C:N" | ttype == "N:N",])

ggplot(df, aes(x=ttype, y=score)) + geom_boxplot(notch = T) + 
  theme_bw()  + xlab("") + ylab("score") + 
  theme(legend.position="bottom", axis.title = element_text(size=16), 
        axis.text= element_text(size=16), panel.grid.minor=element_blank(),
        plot.background=element_blank(), panel.background=element_blank(), panel.grid.major=element_blank())


# combine with triples
df = rbind(df, topTriples[,c("pair", "type", "f1", "prec", "rec", "ttype", "db", "dist.man", "score")])

# plot of C, C:C, C:C:C
gg <- ggplot(df[ttype %in% c("C", "C:C", "C:C:C")], aes(x=ttype, y=f1, fill=ttype)) + 
  geom_boxplot(notch = T, outlier.alpha = 0) + 
  theme_bw() + ylim(0,1) + xlab("") + ylab("F1 score") + scale_fill_manual(values=c("#B2DEDB", "#029688", "#004D41")) +
  theme(legend.position="none", axis.title = element_text(size=16), 
        axis.text= element_text(size=18), panel.grid.minor=element_blank(),
        plot.background=element_blank(), panel.background=element_blank(), panel.grid.major=element_blank())

ggsave(paste0("scripts/clustering-method/res-figures/paper/one-two-three-clins-boxplot.pdf"), 
       gg, width=1.75, height=2.25, units="in", device="pdf", useDingbats=FALSE)


# plot of C, N, C:C, C:N, N:N
df[["ttype"]] <- factor(df$ttype, levels = c("C", "N", "C:C", "C:N", "N:N", "C:C:C"))
gg <- ggplot(df[ttype == "C" | ttype == "N" | ttype == "C:C" | ttype == "C:N" | ttype == "N:N"], aes(x=ttype, y=f1)) + 
  geom_boxplot(notch = T, outlier.alpha = 0) + 
  theme_bw() + ylim(0,1) + xlab("") + ylab("F1 score") + scale_fill_manual(values=c("#B2DEDB", "#FFD550", "#029688", "#8BC24B", "#FF8F00")) +
  theme(legend.position="none", axis.title = element_text(size=16), 
        axis.text= element_text(size=18), panel.grid.minor=element_blank(),
        plot.background=element_blank(), panel.background=element_blank(), panel.grid.major=element_blank())

ggsave(paste0("scripts/clustering-method/res-figures/paper/n-increase-performance-boxplot.pdf"), 
       gg, width=2, height=2.25, units="in", device="pdf", useDingbats=FALSE)


gg <- ggplot(df[ttype == "C" | ttype == "C:C" | ttype == "C:N" | ttype == "N:N"], aes(x=ttype, y=f1)) + 
  geom_boxplot(notch = T, outlier.alpha = 0) + 
  theme_bw() + ylim(0,1) + xlab("") + ylab("F1 score") + scale_fill_manual(values=c("#B2DEDB", "#FFD550", "#029688", "#8BC24B", "#FF8F00")) +
  theme(legend.position="none", axis.title = element_text(size=16), 
        axis.text= element_text(size=18), panel.grid.minor=element_blank(),
        plot.background=element_blank(), panel.background=element_blank(), panel.grid.major=element_blank())

ggsave(paste0("scripts/clustering-method/res-figures/paper/withoutn-increase-performance-boxplot.pdf"), 
       gg, width=2, height=2.25, units="in", device="pdf", useDingbats=FALSE)


# plot of single, to dual
df[, count := ifelse((ttype == "C" | ttype == "N"), "single", "dual")]
df[, count := ifelse(ttype == "C:C:C", "triple", count)]

df.tmp = df[ttype == "N",]
df.tmp = df.tmp[order(-score),]
df.tmp = df.tmp[, head(.SD, 30), by=type]
df.tmp = rbind(df.tmp, df[ttype == "C",])
df.tmp = rbind(df.tmp, df[count=="dual",])

df.tmp[["count"]] <- factor(df.tmp$count, levels = c("single", "dual"))

gg <- ggplot(df.tmp, aes(x=count, y=f1)) + 
  geom_boxplot(notch = T, outlier.alpha = 0) + 
  theme_bw() + ylim(0,1) + xlab("") + ylab("F1 score") + scale_fill_manual(values=c("#B2DEDB", "#029688", "#004D41")) +
  theme(legend.position="none", axis.title = element_text(size=16), 
        axis.text= element_text(size=18), panel.grid.minor=element_blank(),
        plot.background=element_blank(), panel.background=element_blank(), panel.grid.major=element_blank())

ggsave(paste0("scripts/clustering-method/res-figures/paper/one-two-overall-boxplot.pdf"), 
       gg, width=1.75, height=2.25, units="in", device="pdf", useDingbats=FALSE)



# lollipop plot best single to best double
# simple just show best C pair to best double

# TEMP CLEAN UP!
dbl = df[count == "dual",]
sin = df[count == "single",]
maxDPair = dbl[dbl[, .I[which.max(f1)], by=type]$V1]
maxClin = sin[sin[ttype=="C", .I[which.max(f1)], by=type]$V1]
maxPairs = merge(maxDPair, maxClin, by="type")

maxPairs$type <- factor(maxPairs$type, levels=rev(unique(maxPairs[order(-f1.x), ]$type)))
ggplot(maxPairs) + geom_segment(aes(x=type, xend=type, y=f1.y, yend=f1.x), size=1.5, color="#dadfe3") +
  geom_point(aes(x=type, y=f1.x), color="#16a085", size=2 ) +
  geom_point(aes(x=type, y=f1.y), color="#8e44ad", size=2 ) + xlab("") + ylab("F1") + 
  scale_x_discrete(label=function(x) abbreviate(x, minlength=15)) + coord_flip() + theme_minimal()

# more complicated - show exactly which type of pair is the best double
gglolli = ggplot(maxPairs) + geom_segment(aes(x=type, xend=type, y=f1.y, yend=f1.x), size=1.5, color="#dadfe3") +
  geom_point(aes(x=type, y=f1.x, color=ttype.x), size=2 ) + scale_color_manual(values=c("#e67e22", "#16a085", "#2980b9")) +
  geom_text(aes(label=pair.x, x=type, y=f1.x),hjust=1.3, vjust=0.5)  +
  geom_point(aes(x=type, y=f1.y), color="#8e44ad", size=2) +
  geom_text(aes(label=pair.y, x=type, y=f1.y),hjust=-1.1, vjust=0.5) + xlab("") + ylab("F1") + 
  scale_x_discrete(label=function(x) abbreviate(x, minlength=15)) + coord_flip() + theme_bw()

ggsave(paste0("scripts/clustering-method/res-figures/paper/maxF1Lollipop.pdf"), 
       gglolli, width=5, height=6, units="in", device="pdf", useDingbats=FALSE)

maxPairs[, dF1 := f1.x - f1.y]
maxPairs[, dprec := prec.x - prec.y]
maxPairs[, drec := rec.x - rec.y]
maxPairs[, du := uth.y - uth.x]

# improvement 1 vs 2 antigens
df[, antigens := ifelse(ttype == "C", "single", "double")]
df$antigens <- factor(df$antigens, levels=c("single", "double"))
ggplot(df, aes(x=antigens, y=f1)) + geom_boxplot(notch = T) + 
  theme_bw() + ylim(0,1) + xlab("") + ylab("F1 score") +
  theme(legend.position="bottom", axis.title = element_text(size=16), 
        axis.text= element_text(size=16), panel.grid.minor=element_blank(),
        plot.background=element_blank(), panel.background=element_blank(), panel.grid.major=element_blank())
wilcox.test(f1 ~ antigens, data=df)




# add triples - since limited set of C:C:Cs don't enforce the maximum of 2 pairings per gene b/c that seems weird
# but still take only the top 10 per cancer - 4000 per tumor type
tripleres[,can := type]
memm = data.table(pair="dummy", val="d")


# testttttt to toss!!
testset = df[ttype=="N:N",]
testset = df


#temp updated clins
cs = merge(cSingles, singles2, by=c("can", "gene")) 


summary(testdf.tmp7$f1)
testdf.tmp77 = merge(df[,c("type", "pair", "db", "dist.man", "score", "ttype")], testdf.tmp7, by.x=c("pair", "type"), by.y=c("pair", "can"))

#df.tmp = radial with proportional weights
#df.tmp2 = radial no weights
#df.tmp3 = radial 0.7 / 0.3
ggplot(testdf.tmp77[ttype %in% c("C", "C:C", "C:C:C")], aes(x=ttype, y=f1, fill=ttype)) + 
  geom_boxplot(notch = T, outlier.alpha = 0) + 
  theme_bw() + ylim(0,1) + xlab("") + ylab("F1 score") + scale_fill_manual(values=c("#B2DEDB", "#029688", "#004D41")) +
  theme(legend.position="none", axis.title = element_text(size=16), 
        axis.text= element_text(size=18), panel.grid.minor=element_blank(),
        plot.background=element_blank(), panel.background=element_blank(), panel.grid.major=element_blank())


#df[["ttype"]] <- factor(df$ttype, levels = c("C", "N", "C:C", "C:N", "N:N", "C:C:C"))
gg <- ggplot(testdf.tmp77[ttype == "C" | ttype == "N" | ttype == "C:C" | ttype == "C:N" | ttype == "N:N"], aes(x=ttype, y=f1)) + 
  geom_boxplot(notch = T, outlier.alpha = 0) + 
  theme_bw() + ylim(0,1) + xlab("") + ylab("F1 score") + scale_fill_manual(values=c("#B2DEDB", "#FFD550", "#029688", "#8BC24B", "#FF8F00")) +
  theme(legend.position="none", axis.title = element_text(size=16), 
        axis.text= element_text(size=18), panel.grid.minor=element_blank(),
        plot.background=element_blank(), panel.background=element_blank(), panel.grid.major=element_blank())


# LETS SEE MORE PLOTS BLAH BLAH BLAH - PRECSION
# plot of C, C:C, C:C:C
gg <- ggplot(df[ttype %in% c("C", "C:C", "C:C:C")], aes(x=ttype, y=prec, fill=ttype)) + 
  geom_boxplot(notch = T, outlier.alpha = 0) + 
  theme_bw() + ylim(0,1) + xlab("") + ylab("precision") + scale_fill_manual(values=c("#B2DEDB", "#029688", "#004D41")) +
  theme(legend.position="none", axis.title = element_text(size=13), 
        axis.text= element_text(size=13), panel.grid.minor=element_blank(),
        plot.background=element_blank(), panel.background=element_blank(), panel.grid.major=element_blank())

ggsave(paste0("scripts/clustering-method/res-figures/paper/one-two-three-clins-boxplot-prec.pdf"), 
       gg, width=1.75, height=2.25, units="in", device="pdf", useDingbats=FALSE)


# plot of C, N, C:C, C:N, N:N
df[["ttype"]] <- factor(df$ttype, levels = c("C", "N", "C:C", "C:N", "N:N", "C:C:C"))
gg <- ggplot(df[ttype == "C" | ttype == "N" | ttype == "C:C" | ttype == "C:N" | ttype == "N:N"], aes(x=ttype, y=prec)) + 
  geom_boxplot(notch = T, outlier.alpha = 0) + 
  theme_bw() + ylim(0,1) + xlab("") + ylab("precision") + scale_fill_manual(values=c("#B2DEDB", "#FFD550", "#029688", "#8BC24B", "#FF8F00")) +
  theme(legend.position="none", axis.title = element_text(size=10), 
        axis.text= element_text(size=10), panel.grid.minor=element_blank(),
        plot.background=element_blank(), panel.background=element_blank(), panel.grid.major=element_blank())

ggsave(paste0("scripts/clustering-method/res-figures/paper/n-increase-performance-boxplot-prec.pdf"), 
       gg, width=2, height=2.25, units="in", device="pdf", useDingbats=FALSE)


gg <- ggplot(df[ttype == "C" | ttype == "C:C" | ttype == "C:N" | ttype == "N:N"], aes(x=ttype, y=prec)) + 
  geom_boxplot(notch = T, outlier.alpha = 0) + 
  theme_bw() + ylim(0,1) + xlab("") + ylab("precison") + scale_fill_manual(values=c("#B2DEDB", "#FFD550", "#029688", "#8BC24B", "#FF8F00")) +
  theme(legend.position="none", axis.title = element_text(size=13), 
        axis.text= element_text(size=13), panel.grid.minor=element_blank(),
        plot.background=element_blank(), panel.background=element_blank(), panel.grid.major=element_blank())

ggsave(paste0("scripts/clustering-method/res-figures/paper/withoutn-increase-performance-boxplot-prec.pdf"), 
       gg, width=2, height=2.25, units="in", device="pdf", useDingbats=FALSE)


# plot of single, to dual
gg <- ggplot(df.tmp, aes(x=count, y=prec)) + 
  geom_boxplot(notch = T, outlier.alpha = 0) + 
  theme_bw() + ylim(0,1) + xlab("") + ylab("precision") + scale_fill_manual(values=c("#B2DEDB", "#029688", "#004D41")) +
  theme(legend.position="none", axis.title = element_text(size=16), 
        axis.text= element_text(size=18), panel.grid.minor=element_blank(),
        plot.background=element_blank(), panel.background=element_blank(), panel.grid.major=element_blank())

ggsave(paste0("scripts/clustering-method/res-figures/paper/one-two-overall-boxplot-prec.pdf"), 
       gg, width=1.75, height=2.25, units="in", device="pdf", useDingbats=FALSE)


# more complicated - show exactly which type of pair is the best double
maxDPair = dbl[dbl[, .I[which.max(prec)], by=type]$V1]
maxClin = sin[sin[ttype=="C", .I[which.max(prec)], by=type]$V1]
maxPairs = merge(maxDPair, maxClin, by="type")

maxPairs$type <- factor(maxPairs$type, levels=rev(unique(maxPairs[order(-prec.x), ]$type)))
gglolli = ggplot(maxPairs) + geom_segment(aes(x=type, xend=type, y=prec.y, yend=prec.x), size=1.5, color="#dadfe3") +
  geom_point(aes(x=type, y=prec.x, color=ttype.x), size=2 ) + scale_color_manual(values=c("#e67e22", "#16a085", "#2980b9")) +
  #geom_text(aes(label=pair.x, x=type, y=prec.x),hjust=1.3, vjust=0.5)  +
  geom_point(aes(x=type, y=prec.y), color="#8e44ad", size=2) +
  #geom_text(aes(label=pair.y, x=type, y=prec.y),hjust=-1.1, vjust=0.5)
  xlab("") + ylab("precision") + 
  scale_x_discrete(label=function(x) abbreviate(x, minlength=15)) + coord_flip() + theme_bw()
gglolli
ggsave(paste0("scripts/clustering-method/res-figures/paper/maxPrecLollipop.pdf"), 
       gglolli, width=5, height=6, units="in", device="pdf", useDingbats=FALSE)



# LETS SEE MORE PLOTS BLAH BLAH BLAH - RECALL
# plot of C, C:C, C:C:C
gg <- ggplot(df[ttype %in% c("C", "C:C", "C:C:C")], aes(x=ttype, y=rec, fill=ttype)) + 
  geom_boxplot(notch = T, outlier.alpha = 0) + 
  theme_bw() + ylim(0,1) + xlab("") + ylab("recall") + scale_fill_manual(values=c("#B2DEDB", "#029688", "#004D41")) +
  theme(legend.position="none", axis.title = element_text(size=13), 
        axis.text= element_text(size=13), panel.grid.minor=element_blank(),
        plot.background=element_blank(), panel.background=element_blank(), panel.grid.major=element_blank())
gg
ggsave(paste0("scripts/clustering-method/res-figures/paper/one-two-three-clins-boxplot-recall.pdf"), 
       gg, width=1.75, height=2.25, units="in", device="pdf", useDingbats=FALSE)


# plot of C, N, C:C, C:N, N:N
df[["ttype"]] <- factor(df$ttype, levels = c("C", "N", "C:C", "C:N", "N:N", "C:C:C"))
gg <- ggplot(df[ttype == "C" | ttype == "N" | ttype == "C:C" | ttype == "C:N" | ttype == "N:N"], aes(x=ttype, y=rec)) + 
  geom_boxplot(notch = T, outlier.alpha = 0) + 
  theme_bw() + ylim(0,1) + xlab("") + ylab("recall") + scale_fill_manual(values=c("#B2DEDB", "#FFD550", "#029688", "#8BC24B", "#FF8F00")) +
  theme(legend.position="none", axis.title = element_text(size=13), 
        axis.text= element_text(size=13), panel.grid.minor=element_blank(),
        plot.background=element_blank(), panel.background=element_blank(), panel.grid.major=element_blank())
gg
ggsave(paste0("scripts/clustering-method/res-figures/paper/n-increase-performance-boxplot-recall.pdf"), 
       gg, width=2, height=2.25, units="in", device="pdf", useDingbats=FALSE)


gg <- ggplot(df[ttype == "C" | ttype == "C:C" | ttype == "C:N" | ttype == "N:N"], aes(x=ttype, y=rec)) + 
  geom_boxplot(notch = T, outlier.alpha = 0) + 
  theme_bw() + ylim(0,1) + xlab("") + ylab("recall") + scale_fill_manual(values=c("#B2DEDB", "#FFD550", "#029688", "#8BC24B", "#FF8F00")) +
  theme(legend.position="none", axis.title = element_text(size=13), 
        axis.text= element_text(size=13), panel.grid.minor=element_blank(),
        plot.background=element_blank(), panel.background=element_blank(), panel.grid.major=element_blank())

gg
ggsave(paste0("scripts/clustering-method/res-figures/paper/withoutn-increase-performance-boxplot-recall.pdf"), 
       gg, width=2, height=2.25, units="in", device="pdf", useDingbats=FALSE)


# plot of single, to dual
gg <- ggplot(df.tmp, aes(x=count, y=rec)) + 
  geom_boxplot(notch = T, outlier.alpha = 0) + 
  theme_bw() + ylim(0,1) + xlab("") + ylab("recall") + scale_fill_manual(values=c("#B2DEDB", "#029688", "#004D41")) +
  theme(legend.position="none", axis.title = element_text(size=16), 
        axis.text= element_text(size=18), panel.grid.minor=element_blank(),
        plot.background=element_blank(), panel.background=element_blank(), panel.grid.major=element_blank())
gg
ggsave(paste0("scripts/clustering-method/res-figures/paper/one-two-overall-boxplot-recall.pdf"), 
       gg, width=1.75, height=2.25, units="in", device="pdf", useDingbats=FALSE)


# more complicated - show exactly which type of pair is the best double
maxDPair = dbl[dbl[, .I[which.max(f1)], by=type]$V1]
maxClin = sin[sin[ttype=="C", .I[which.max(f1)], by=type]$V1]
maxPairs = merge(maxDPair, maxClin, by="type")

maxPairs$type <- factor(maxPairs$type, levels=rev(unique(maxPairs[order(-rec.x), ]$type)))
gglolli = ggplot(maxPairs) + geom_segment(aes(x=type, xend=type, y=rec.y, yend=rec.x), size=1.5, color="#dadfe3") +
  geom_point(aes(x=type, y=rec.x, color=ttype.x), size=2 ) + scale_color_manual(values=c("#e67e22", "#16a085", "#2980b9")) +
  #geom_text(aes(label=pair.x, x=type, y=rec.x),hjust=1.3, vjust=0.5)  +
  geom_point(aes(x=type, y=rec.y), color="#8e44ad", size=2) +
  #geom_text(aes(label=pair.y, x=type, y=rec.y),hjust=-1.1, vjust=0.5) + 
  xlab("") + ylab("recall") + 
  scale_x_discrete(label=function(x) abbreviate(x, minlength=15)) + coord_flip() + theme_bw()
gglolli
ggsave(paste0("scripts/clustering-method/res-figures/paper/maxRecallLollipop.pdf"), 
       gglolli, width=5, height=6, units="in", device="pdf", useDingbats=FALSE)

# LETS SEE MORE PLOTS BLAH BLAH BLAH - Unique Tissues Hit
# plot of C, C:C, C:C:C
gg <- ggplot(df[ttype %in% c("C", "C:C", "C:C:C")], aes(x=ttype, y=uth, fill=ttype)) + 
  geom_boxplot(notch = T, outlier.alpha = 0) + 
  theme_bw() + ylim(0,33) + xlab("") + ylab("unique tissues hit") + scale_fill_manual(values=c("#B2DEDB", "#029688", "#004D41")) +
  theme(legend.position="none", axis.title = element_text(size=13), 
        axis.text= element_text(size=13), panel.grid.minor=element_blank(),
        plot.background=element_blank(), panel.background=element_blank(), panel.grid.major=element_blank())
gg
ggsave(paste0("scripts/clustering-method/res-figures/paper/one-two-three-clins-boxplot-uth.pdf"), 
       gg, width=1.75, height=2.25, units="in", device="pdf", useDingbats=FALSE)


# plot of C, N, C:C, C:N, N:N
df[["ttype"]] <- factor(df$ttype, levels = c("C", "N", "C:C", "C:N", "N:N", "C:C:C"))
gg <- ggplot(df[ttype == "C" | ttype == "N" | ttype == "C:C" | ttype == "C:N" | ttype == "N:N"], aes(x=ttype, y=uth)) + 
  geom_boxplot(notch = T, outlier.alpha = 0) + 
  theme_bw() + ylim(0,33) + xlab("") + ylab("unique tissues hit") + scale_fill_manual(values=c("#B2DEDB", "#FFD550", "#029688", "#8BC24B", "#FF8F00")) +
  theme(legend.position="none", axis.title = element_text(size=13), 
        axis.text= element_text(size=13), panel.grid.minor=element_blank(),
        plot.background=element_blank(), panel.background=element_blank(), panel.grid.major=element_blank())
gg
ggsave(paste0("scripts/clustering-method/res-figures/paper/n-increase-performance-boxplot-uth.pdf"), 
       gg, width=2, height=2.25, units="in", device="pdf", useDingbats=FALSE)


gg <- ggplot(df[ttype == "C" | ttype == "C:C" | ttype == "C:N" | ttype == "N:N"], aes(x=ttype, y=uth)) + 
  geom_boxplot(notch = T, outlier.alpha = 0) + 
  theme_bw() + ylim(0,33) + xlab("") + ylab("unique tissues hit") + scale_fill_manual(values=c("#B2DEDB", "#FFD550", "#029688", "#8BC24B", "#FF8F00")) +
  theme(legend.position="none", axis.title = element_text(size=13), 
        axis.text= element_text(size=13), panel.grid.minor=element_blank(),
        plot.background=element_blank(), panel.background=element_blank(), panel.grid.major=element_blank())
gg
ggsave(paste0("scripts/clustering-method/res-figures/paper/withoutn-increase-performance-boxplot-uth.pdf"), 
       gg, width=2, height=2.25, units="in", device="pdf", useDingbats=FALSE)


# plot of single, to dual
gg <- ggplot(df.tmp, aes(x=count, y=uth)) + 
  geom_boxplot(notch = T, outlier.alpha = 0) + 
  theme_bw() + ylim(0,33) + xlab("") + ylab("unique tissues hit") + scale_fill_manual(values=c("#B2DEDB", "#029688", "#004D41")) +
  theme(legend.position="none", axis.title = element_text(size=16), 
        axis.text= element_text(size=18), panel.grid.minor=element_blank(),
        plot.background=element_blank(), panel.background=element_blank(), panel.grid.major=element_blank())
gg
ggsave(paste0("scripts/clustering-method/res-figures/paper/one-two-overall-boxplot-uth.pdf"), 
       gg, width=1.75, height=2.25, units="in", device="pdf", useDingbats=FALSE)


maxDPair = dbl[dbl[, .I[which.max(f1)], by=type]$V1]
maxClin = sin[sin[ttype=="C", .I[which.max(f1)], by=type]$V1]
maxPairs = merge(maxDPair, maxClin, by="type")

maxPairs$type <- factor(maxPairs$type, levels=rev(unique(maxPairs[order(uth.x), ]$type)))
gglolli = ggplot(maxPairs) + geom_segment(aes(x=type, xend=type, y=uth.y, yend=uth.x), size=1.5, color="#dadfe3") +
  geom_point(aes(x=type, y=uth.x, color=ttype.x), size=2 ) + scale_color_manual(values=c("#e67e22", "#16a085", "#2980b9")) +
  #geom_text(aes(label=pair.x, x=type, y=uth.x),hjust=1.3, vjust=0.5)  +
  geom_point(aes(x=type, y=uth.y), color="#8e44ad", size=2) +
  #geom_text(aes(label=pair.y, x=type, y=uth.y),hjust=-1.1, vjust=0.5) + xlab("") + ylab("unique tissues hit") + 
  scale_x_discrete(label=function(x) abbreviate(x, minlength=15)) + coord_flip() + theme_bw()
gglolli
ggsave(paste0("scripts/clustering-method/res-figures/paper/maxUthLollipop.pdf"), 
       gglolli, width=5, height=6, units="in", device="pdf", useDingbats=FALSE)

# LETS SEE MORE PLOTS BLAH BLAH BLAH - Accuracy
# plot of C, C:C, C:C:C
gg <- ggplot(df[ttype %in% c("C", "C:C", "C:C:C")], aes(x=ttype, y=acc, fill=ttype)) + 
  geom_boxplot(notch = T, outlier.alpha = 0) + 
  theme_bw() + ylim(0,1) + xlab("") + ylab("accuracy") + scale_fill_manual(values=c("#B2DEDB", "#029688", "#004D41")) +
  theme(legend.position="none", axis.title = element_text(size=13), 
        axis.text= element_text(size=13), panel.grid.minor=element_blank(),
        plot.background=element_blank(), panel.background=element_blank(), panel.grid.major=element_blank())
gg
ggsave(paste0("scripts/clustering-method/res-figures/paper/one-two-three-clins-boxplot-acc.pdf"), 
       gg, width=1.75, height=2.25, units="in", device="pdf", useDingbats=FALSE)


# plot of C, N, C:C, C:N, N:N
df[["ttype"]] <- factor(df$ttype, levels = c("C", "N", "C:C", "C:N", "N:N", "C:C:C"))
gg <- ggplot(df[ttype == "C" | ttype == "N" | ttype == "C:C" | ttype == "C:N" | ttype == "N:N"], aes(x=ttype, y=acc)) + 
  geom_boxplot(notch = T, outlier.alpha = 0) + 
  theme_bw() + ylim(0,1) + xlab("") + ylab("accuracy") + scale_fill_manual(values=c("#B2DEDB", "#FFD550", "#029688", "#8BC24B", "#FF8F00")) +
  theme(legend.position="none", axis.title = element_text(size=13), 
        axis.text= element_text(size=13), panel.grid.minor=element_blank(),
        plot.background=element_blank(), panel.background=element_blank(), panel.grid.major=element_blank())
gg
ggsave(paste0("scripts/clustering-method/res-figures/paper/n-increase-performance-boxplot-acc.pdf"), 
       gg, width=2, height=2.25, units="in", device="pdf", useDingbats=FALSE)

# plot of single, to dual
gg <- ggplot(df.tmp, aes(x=count, y=acc)) + 
  geom_boxplot(notch = T, outlier.alpha = 0) + 
  theme_bw() + ylim(0,1) + xlab("") + ylab("accuracy") + scale_fill_manual(values=c("#B2DEDB", "#029688", "#004D41")) +
  theme(legend.position="none", axis.title = element_text(size=13), 
        axis.text= element_text(size=13), panel.grid.minor=element_blank(),
        plot.background=element_blank(), panel.background=element_blank(), panel.grid.major=element_blank())
gg
ggsave(paste0("scripts/clustering-method/res-figures/paper/one-two-overall-boxplot-acc.pdf"), 
       gg, width=1.75, height=2.25, units="in", device="pdf", useDingbats=FALSE)


# LETS SEE MORE PLOTS BLAH BLAH BLAH - clustering scores - DB
# plot of C, C:C, C:C:C
gg <- ggplot(df[ttype %in% c("C", "C:C", "C:C:C")], aes(x=ttype, y=log(db), fill=ttype)) + 
  geom_boxplot(notch = T, outlier.alpha = 0) + 
  theme_bw()  + xlab("") + ylab("log(Davies-Bouldin)") + scale_fill_manual(values=c("#B2DEDB", "#029688", "#004D41")) +
  theme(legend.position="none", axis.title = element_text(size=13), 
        axis.text= element_text(size=13), panel.grid.minor=element_blank(),
        plot.background=element_blank(), panel.background=element_blank(), panel.grid.major=element_blank())
gg
ggsave(paste0("scripts/clustering-method/res-figures/paper/one-two-three-clins-boxplot-db.pdf"), 
       gg, width=1.75, height=2.25, units="in", device="pdf", useDingbats=FALSE)


# plot of C, N, C:C, C:N, N:N
df[["ttype"]] <- factor(df$ttype, levels = c("C", "N", "C:C", "C:N", "N:N", "C:C:C"))
gg <- ggplot(df[ttype == "C" | ttype == "N" | ttype == "C:C" | ttype == "C:N" | ttype == "N:N"], aes(x=ttype, y=log(db))) + 
  geom_boxplot(notch = T, outlier.alpha = 0) + 
  theme_bw()  + xlab("") + ylab("log(Davies-Bouldin)") + scale_fill_manual(values=c("#B2DEDB", "#FFD550", "#029688", "#8BC24B", "#FF8F00")) +
  theme(legend.position="none", axis.title = element_text(size=13), 
        axis.text= element_text(size=13), panel.grid.minor=element_blank(),
        plot.background=element_blank(), panel.background=element_blank(), panel.grid.major=element_blank())
gg
ggsave(paste0("scripts/clustering-method/res-figures/paper/n-increase-performance-boxplot-db.pdf"), 
       gg, width=2, height=2.25, units="in", device="pdf", useDingbats=FALSE)


# plot of single, to dual
gg <- ggplot(df.tmp, aes(x=count, y=log(db))) + 
  geom_boxplot(notch = T, outlier.alpha = 0) + 
  theme_bw()  + xlab("") + ylab("log(Davies-Bouldin)") + scale_fill_manual(values=c("#B2DEDB", "#029688", "#004D41")) +
  theme(legend.position="none", axis.title = element_text(size=16), 
        axis.text= element_text(size=18), panel.grid.minor=element_blank(),
        plot.background=element_blank(), panel.background=element_blank(), panel.grid.major=element_blank())
gg
ggsave(paste0("scripts/clustering-method/res-figures/paper/one-two-overall-boxplot-db.pdf"), 
       gg, width=1.75, height=2.25, units="in", device="pdf", useDingbats=FALSE)


# LETS SEE MORE PLOTS BLAH BLAH BLAH - clustering scores - distance
# plot of C, C:C, C:C:C
gg <- ggplot(df[ttype %in% c("C", "C:C", "C:C:C")], aes(x=ttype, y=dist.man, fill=ttype)) + 
  geom_boxplot(notch = T, outlier.alpha = 0) + 
  theme_bw()  + xlab("") + ylab("distance") + scale_fill_manual(values=c("#B2DEDB", "#029688", "#004D41")) +
  theme(legend.position="none", axis.title = element_text(size=13), 
        axis.text= element_text(size=13), panel.grid.minor=element_blank(),
        plot.background=element_blank(), panel.background=element_blank(), panel.grid.major=element_blank())
gg
ggsave(paste0("scripts/clustering-method/res-figures/paper/one-two-three-clins-boxplot-dist.pdf"), 
       gg, width=1.75, height=2.25, units="in", device="pdf", useDingbats=FALSE)


# plot of C, N, C:C, C:N, N:N
df[["ttype"]] <- factor(df$ttype, levels = c("C", "N", "C:C", "C:N", "N:N", "C:C:C"))
gg <- ggplot(df[ttype == "C" | ttype == "N" | ttype == "C:C" | ttype == "C:N" | ttype == "N:N"], aes(x=ttype, y=dist.man)) + 
  geom_boxplot(notch = T, outlier.alpha = 0) + 
  theme_bw()  + xlab("") + ylab("distance") + scale_fill_manual(values=c("#B2DEDB", "#FFD550", "#029688", "#8BC24B", "#FF8F00")) +
  theme(legend.position="none", axis.title = element_text(size=13), 
        axis.text= element_text(size=13), panel.grid.minor=element_blank(),
        plot.background=element_blank(), panel.background=element_blank(), panel.grid.major=element_blank())
gg
ggsave(paste0("scripts/clustering-method/res-figures/paper/n-increase-performance-boxplot-dist.pdf"), 
       gg, width=2, height=2.25, units="in", device="pdf", useDingbats=FALSE)


# plot of single, to dual
gg <- ggplot(df.tmp, aes(x=count, y=dist.man)) + 
  geom_boxplot(notch = T, outlier.alpha = 0) + 
  theme_bw()  + xlab("") + ylab("distance") + scale_fill_manual(values=c("#B2DEDB", "#029688", "#004D41")) +
  theme(legend.position="none", axis.title = element_text(size=16), 
        axis.text= element_text(size=18), panel.grid.minor=element_blank(),
        plot.background=element_blank(), panel.background=element_blank(), panel.grid.major=element_blank())
gg
ggsave(paste0("scripts/clustering-method/res-figures/paper/one-two-overall-boxplot-dist.pdf"), 
       gg, width=1.75, height=2.25, units="in", device="pdf", useDingbats=FALSE)

# save table of the top Cs and the top k (10) pairs per type per tumor
write.table(df[ttype %in% c("C", "C:C", "C:N", "N:N"),], paste0("scripts/clustering-method/res-tables/top-",k,"-clins-and-pairs.txt"), sep="\t", quote=F, row.names=F)

# check how many pairs fall above x
bigres.clean=bigres.clean[score > 0.85,] # CLEAN - use bigres
brc = rbindlist(mclapply(cans, function(cn) { 
  sketchL = sketchLAll[can==cn]$dataset
  sketch = rnaseqmat[dataset %in% sketchL,]
  sketch = sketch[sketch$type == "tissue" | sketch$tissue.cancer == cn,]
  sketch[, "target" := ifelse(sketch$tissue.cancer == cn, "target", "other")]
  tres = bigres.clean[type==cn]
  ctmp = apply(tres, 1, GetClassType4, tmp=sketch)
  tres[, ctype := ctmp]
  tres
}, mc.cores=15))

# remove all LL pairs from bigres for sara
br = rbindlist(mclapply(cans, function(cn) { 
  sketchL = sketchLAll[can==cn]$dataset
  sketch = rnaseqmat[dataset %in% sketchL,]
  sketch = sketch[sketch$type == "tissue" | sketch$tissue.cancer == cn,]
  sketch[, "target" := ifelse(sketch$tissue.cancer == cn, "target", "other")]
  tres = bigres[type==cn]
  ctmp = apply(tres, 1, GetClassType4, tmp=sketch)
  tres[, ctype := ctmp]
  tres
}, mc.cores=15))
brdone = br[ctype != "LL",]
saveRDS(brdone, "scripts/clustering-method/rdata/bigres-no-LL-for-sara.RData")


bigres.clean = brc[ctype != "LL",]
bigres.clean[, ctype := NULL]
bigres.clean = fintersect(bigres.clean[,c("pair", "type")], bigres[score > 0.85, c("pair", "type")])
bigres.clean = merge(bigres.clean, bigres, by=c("pair", "type"))
bigres.clean[,.N, by=c("type")]

tab85 = bigres.clean[score > 0.85,.N, by=c("type")]
tab90 = bigres.clean[score > 0.9,.N, by=c("type")]
tab = merge(tab85, tab90, by="type")
names(tab) = c("cancer", "N85", "N90") # BUSTEDDDD b/c not everything is mergedddd
tab = melt(tab, id.vars = c("cancer"))
tab$value = log10(tab$value)

tab$variable = factor(tab$variable, levels = c("N90", "N85"))
gg <- ggplot(tab, aes(reorder(cancer, value), value, fill=variable)) + geom_bar(stat="identity", position="dodge") + 
  facet_wrap(~variable) +
  theme_bw() + coord_flip() + xlab("") + ylab("log10(number of pairs)") + 
  #scale_fill_manual(values = c("#172869", "#1BB6AF")) +
  theme(legend.position="none", axis.text= element_text(size=9), panel.grid.minor=element_blank(),
      plot.background=element_blank(), panel.background=element_blank(), panel.grid.major=element_blank()) #+ 
  #scale_y_continuous(labels = function(x) format(x, scientific = TRUE))
gg  
ggsave(paste0("scripts/clustering-method/res-figures/paper/numOfPairsPerThreshold.pdf"), 
       gg, width=7, height=5.5, units="in", device="pdf", useDingbats=FALSE)



tab85$N = log10(tab85$N)
gg <- ggplot(tab85, aes(reorder(type, N), N)) + geom_bar(stat="identity", position="dodge") + 
  theme_bw() + coord_flip() + xlab("") + ylab("log10(number of pairs)") + 
  #scale_fill_manual(values = c("#172869", "#1BB6AF")) +
  theme(legend.position="none", axis.text= element_text(size=9), panel.grid.minor=element_blank(),
        plot.background=element_blank(), panel.background=element_blank(), panel.grid.major=element_blank()) #+ 
#scale_y_continuous(labels = function(x) format(x, scientific = TRUE))
gg  
ggsave(paste0("scripts/clustering-method/res-figures/paper/numOfPairsPerThreshold85.pdf"), 
       gg, width=7, height=5.5, units="in", device="pdf", useDingbats=FALSE)


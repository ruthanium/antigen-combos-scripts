# **************************************
# Author: Ruth Dannenfelser
# Date:   August 26, 2019
#
# generate statistics for clinical
# improvement using radius based metric
#
# **************************************
library(ggplot2)
library(data.table)
library(sp)
library(parallel)
library(e1071)

setwd("/Genomics/ogtr04/rd6/BioGPS")

set.seed(1066)

tm = fread("data/transmembrane-markers/filteredTransmembraneList-updated-12-19.txt", header=F)
rnaseqmat = readRDS("scripts/clustering-method/rdata/combat-rnaseq-tm-clin.RData")
#saveRDS(rnaseqmat, "scripts/clustering-method/rdata/combat-rnaseq-tm-clin-backcompatible.RData", version=2)

cans = unique(rnaseqmat[rnaseqmat$type == "tumor",]$tissue.cancer)

clinSmall = fread("data/transmembrane-markers/clin-paper-dec-19.txt", header=F)
clinSmall = clinSmall$V1
clinSmall = clinSmall[clinSmall %in% names(rnaseqmat)]

tm = tm$V1
tm = tm[tm %in% names(rnaseqmat)]


# # work with all the cancer types - the full set of results (for clinicals)
# #bigres = readRDS("scripts/clustering-method/rdata/all-rnaseq-combat-sketch-0.2.RData") # accio
# bigres = readRDS("scripts/clustering-method/rdata/all-rnaseq-combat-sketch-0.2-top-fixed-backcompatible.RData") # lumos
# #saveRDS(bigres, "scripts/clustering-method/rdata/all-rnaseq-combat-sketch-0.2-top-fixed-backcompatible.RData", version=2)
# 
# # filter out all the genes not in TM from bigres
# bigres.clean = bigres[gene1 %in% tm & gene2 %in% tm,]
# bigres = bigres.clean
# rm(bigres.clean)
# 
# # rescale the scores across all tumor types - only have to do once b/c its saved
# bigres[,dbscaled := 1-((log(db) - min(log(db)))/(max(log(db)) - min(log(db))))] # make score between 0 and 1
# bigres[,distscaled := (log(dist.man) - min(log(dist.man)))/(max(log(dist.man)) - min(log(dist.man)))]
# bigres[,score := apply(bigres[,c("dbscaled", "distscaled"), with=F], 1, min)]
# 
# #saveRDS(bigres, "scripts/clustering-method/rdata/all-rnaseq-combat-sketch-0.2-top-fixed-scaled-over-all-backcompatible.RData", version=2)

bigres = readRDS("scripts/clustering-method/rdata/all-rnaseq-combat-sketch-0.2-top-fixed-scaled-over-all-backcompatible.RData")

sketchLAll = rbindlist(lapply(cans,
                              function(can) data.table(can=can,
                                                       dataset=readLines(paste0("scripts/clustering-method/sketches/", gsub(" ", "-", tolower(can)), "-sketch.txt")))))


tripleres = readRDS("scripts/clustering-method/rdata/all-rnaseq-combat-sketch-0.2-clin-triples.RData")

# calculate metrics for singles
singles
ss = rbindlist(mclapply(cans, function(cann)
{
  f1=numeric()
  prec = numeric()
  rec = numeric()
  uth = numeric()
  tcan = character()
  glist = character()
  acc = numeric()
  for (g in tm)
  {
    acl = assignClusterLabels(cann, c(g), rnaseqmat)
    f1 = append(f1, acl$F1)
    prec = append(prec, acl$prec)
    rec = append(rec, acl$recall)
    uth = append(uth, acl$uniqTisHit)
    tcan = append(tcan, cann)
    glist = append(glist, g)
    acc = append(acc, acl$accuracy)
    
  }
  data.table(gene=glist, can=tcan, f1, prec, rec, uth, acc)
}, mc.cores=15))


singTmp = readRDS("scripts/clustering-method/rdata/singles.RData") 
names(singTmp) = c("can", "gene", "db", "dist.man", "ttype")
singles = merge(singles, singTmp, by=c("can", "gene"))

# label C or S
singles[, glab := ifelse(gene %in% clinSmall, "C", "N") ]





# update bigres typing 
bigres[, g1lab := ifelse(gene1 %in% clinSmall, "C", "N")]
bigres[, g2lab := ifelse(gene2 %in% clinSmall, "C", "N")]
bigres[, ttype := paste0(g1lab, ":", g2lab)]

# get the top k CC, CN, NN pairs per cancer
k = 100 # number of top and bottom to grab
treeTypes = c("C:C", "C:N", "N:N")

topPerClin = rbindlist(lapply(treeTypes, function(tt) {
  rbindlist(mclapply(cans, function(can,tt)
  {
    print(can)
    one = character()
    two = character()
    count = 0
    
    res = bigres[type == can,]
    res = res[ttype == tt]
    
    tpair = character()
    tscore = numeric()
    tcan = character()
    tttype = character()
    tg1 = character()
    tg2 = character()
    tg1lab = character()
    tg2lab = character()
    tctype = character()
    
    i = 1;
    
    sketchL = sketchLAll[can==can]$dataset
    sketch = rnaseqmat[dataset %in% sketchL,]
    
    # pair a maximum of 2 times
    while ((count < k)  & (i < nrow(res)))
    {
      row = res[i,]
      i = i + 1;
      genes = c(row$gene1, row$gene2)
      
      if (!(genes[1] %in% two) & !(genes[2] %in% two)) # if antigen has come up twice print another one for top
      {
        if (genes[1] %in% one)
        {
          two = append(two, genes[1])
          one = one[!one %in% c(genes[1])]
        }
        else
        {
          one = append(one, genes[1])
        }
        if (genes[2] %in% one)
        {
          two = append(two, genes[2])
          one = one[!one %in% c(genes[2])]
        }
        else
        {
          one = append(one, genes[2])
        }
        
        gct = GetClassType2(genes, can, sketch)
        if(gct != "LL")
        {
          count = count + 1;
          tpair = append(tpair, as.character(row$pair))
          tscore = append(tscore, row$score)
          tcan = append(tcan, can)
          tttype = append(tttype, tt)
          tg1 = append(tg1, row$gene1)
          tg2 = append(tg2, row$gene2)
          tg1lab = append(tg1lab, row$g1lab)
          tg2lab = append(tg2lab, row$g2lab)
          tctype = append(tctype, gct)
        }
      }
    }
    data.table(pair=tpair, score=tscore, can=tcan, gene1=tg1, gene2=tg2, g1lab=tg1lab, g2lab=tg2lab, ttype=tttype, ctype=tctype)
    
  }, tt, mc.cores=20))
}))

#topPerClin[,.N,by=.(can,ttype)] #fancy
tpcMetrics = rbindlist(mclapply(cans, function(cann)
{
  tpc = topPerClin[can==cann,]
  f1=numeric()
  prec = numeric()
  rec = numeric()
  uth = numeric()
  tcan = character()
  clist = character()
  acc = numeric()
  for (i in 1:nrow(tpc))
  {
    row = tpc[i,]
    acl = assignClusterLabels(cann, c(row$gene1, row$gene2), rnaseqmat)
    f1 = append(f1, acl$F1)
    prec = append(prec, acl$prec)
    rec = append(rec, acl$recall)
    uth = append(uth, acl$uniqTisHit)
    tcan = append(tcan, cann)
    clist = append(clist, row$pair)
    acc = append(acc, acl$accuracy)
    
  }
  data.table(pair=clist, can=tcan, f1, prec, rec, uth, acc)
}, mc.cores=33))


#tpcMetrics = data.table(pair=clist, can=tcan, f1, prec, rec, uth, acc)
#tpcMetrics = unique(tpcMetrics)
#topPerClin2 = merge(topPerClin, tpcMetrics, by.x=c("pair", "can"), by.y=c("pair", "can"))
#topPerClin2 = unique(topPerClin2)
#topPerClin2 = topPerClin2[order(-f1),]

topPerClin3 = merge(bigres, tpcMetrics, by.x=c("pair", "type"), by.y=c("pair", "can"))
topPerClin3 = topPerClin3[order(-f1),]
saveRDS(topPerClin3, "scripts/clustering-method/rdata/doubles-topPerClin.RData")


# # eliminate L clinicals as those aren't actionable
# tmp = clinMetrics
# names(tmp) = c("pair", "can", "f1", "prec", "rec", "uth", "acc")
# tmp[, ctype := apply(tmp, 1, GetClassType, dt=rnaseqmat)]
# names(tmp) = c("clin", "can", "f1", "prec", "rec", "uth", "acc", "ctype")
# tmp = tmp[ctype != "L",]
# clinMetrics = tmp
# rm(tmp)

maxClin = clinMetrics[clinMetrics[, .I[which.max(f1)], by=can]$V1]

# boxplots of top 5 pairs per clin and top clin
names(topPerClin2) <- c("pair", "can", "db", "dist.man", "dbscaled", "distscaled", "score", "gene1", "gene2", "g1lab", "g2lab", "ttype", "f1", "prec", "rec", "uth", "acc")
topPerClin2$can <- factor(topPerClin2$can, levels=rev(unique(topPerClin2[order(-f1), ]$can)))
ggplot(topPerClin2, aes(y=f1, x=can)) + geom_boxplot(outlier.shape = NA) + theme_minimal() + ylim(0,1) + 
  xlab("") + geom_point(data=maxClin, color="purple") +
  scale_x_discrete(label=function(x) abbreviate(x, minlength=15)) + coord_flip()

topPerClin2$can <- factor(topPerClin2$can, levels=rev(unique(topPerClin2[order(-prec), ]$can)))
ggplot(topPerClin2, aes(y=prec, x=can)) + geom_boxplot(outlier.shape = NA) + theme_minimal() + ylim(0,1) + 
  xlab("") + geom_point(data=maxClin, color="purple") +
  scale_x_discrete(label=function(x) abbreviate(x, minlength=15)) + coord_flip()

topPerClin2$can <- factor(topPerClin2$can, levels=rev(unique(topPerClin2[order(-prec), ]$can)))
ggplot(topPerClin2, aes(y=uth, x=can)) + geom_boxplot(outlier.shape = NA) + theme_minimal() + ylim(0,34) + 
  xlab("") + ylab("unique tissues hit") + geom_point(data=maxClin, color="purple") +
  scale_x_discrete(label=function(x) abbreviate(x, minlength=15)) + coord_flip()


# just show top and top
maxClinPair = topPerClin2[topPerClin2[, .I[which.max(f1)], by=can]$V1]
maxClin = clinMetrics[clinMetrics[, .I[which.max(f1)], by=can]$V1]
maxPairs = merge(maxClinPair, maxClin, by="can")

# simple just show best C pair to best double
maxPairs$can <- factor(maxPairs$can, levels=rev(unique(maxPairs[order(-f1.x), ]$can)))
ggplot(maxPairs) + geom_segment(aes(x=can, xend=can, y=f1.y, yend=f1.x), size=1.5, color="#dadfe3") +
  geom_point(aes(x=can, y=f1.x), color="#16a085", size=2 ) +
  geom_point(aes(x=can, y=f1.y), color="#8e44ad", size=2 ) + xlab("") + ylab("F1") + 
  scale_x_discrete(label=function(x) abbreviate(x, minlength=15)) + coord_flip() + theme_minimal()

# more complicated - show exactly which type of pair is the best double
ggplot(maxPairs) + geom_segment(aes(x=can, xend=can, y=f1.y, yend=f1.x), size=1.5, color="#dadfe3") +
  geom_point(aes(x=can, y=f1.x, color=ttype), size=2 ) + scale_color_manual(values=c("#e67e22", "#16a085", "#2980b9")) +
  geom_text(aes(label=pair, x=can, y=f1.x),hjust=1.3, vjust=0.5)  +
  geom_point(aes(x=can, y=f1.y), color="#8e44ad", size=2) +
  geom_text(aes(label=clin, x=can, y=f1.y),hjust=-1.1, vjust=0.5) + xlab("") + ylab("F1") + 
  scale_x_discrete(label=function(x) abbreviate(x, minlength=15)) + coord_flip() + theme_minimal()

# # instead of showing F1 show prec, accuracy, or unique tissues hit
# maxPairs$can <- factor(maxPairs$can, levels=rev(unique(maxPairs[order(-prec.x), ]$can)))
# ggplot(maxPairs) + geom_segment(aes(x=can, xend=can, y=prec.y, yend=prec.x), size=1.5, color="#dadfe3") +
#   geom_point(aes(x=can, y=prec.x), color="#16a085", size=2 ) +
#   geom_point(aes(x=can, y=prec.y), color="#8e44ad", size=2 ) + xlab("") + ylab("precision") + 
#   scale_x_discrete(label=function(x) abbreviate(x, minlength=15)) + coord_flip() + theme_minimal()
# 
# ggplot(maxPairs) + geom_segment(aes(x=can, xend=can, y=acc.y, yend=acc.x), size=1.5, color="#dadfe3") +
#   geom_point(aes(x=can, y=acc.x), color="#16a085", size=2 ) +
#   geom_point(aes(x=can, y=acc.y), color="#8e44ad", size=2 ) + xlab("") + ylab("accuracy") + ylim(0,1) +
#   scale_x_discrete(label=function(x) abbreviate(x, minlength=15)) + coord_flip() + theme_minimal()
# 
# ggplot(maxPairs) + geom_segment(aes(x=can, xend=can, y=uth.y, yend=uth.x), size=1.5, color="#dadfe3") +
#   geom_point(aes(x=can, y=uth.x), color="#16a085", size=2 ) +
#   geom_point(aes(x=can, y=uth.y), color="#8e44ad", size=2 ) + xlab("") + ylab("unique tissues hit") + ylim(0,34) +
#   scale_x_discrete(label=function(x) abbreviate(x, minlength=15)) + coord_flip() + theme_minimal()


# generate aggregate plots for C, CC, CN, NN
df = clinMetrics
names(df) <- c("pair", "can", "f1", "prec", "rec", "uth", "acc", "ctype")
df[, ttype := "C"]
df = rbind(df, topPerClin2[,c("pair", "can", "f1", "prec", "rec", "uth", "acc", "ttype")], fill=T)

ggplot(df, aes(x=ttype, y=f1)) + geom_boxplot(notch = T) + 
  theme_bw() + ylim(0,1) + xlab("") + ylab("F1 score") + 
  theme(legend.position="bottom", axis.title = element_text(size=16), 
        axis.text= element_text(size=16), panel.grid.minor=element_blank(),
        plot.background=element_blank(), panel.background=element_blank(), panel.grid.major=element_blank())

wilcox.test(f1 ~ ttype, data = df[ttype == "C" | ttype == "C:C",])
wilcox.test(f1 ~ ttype, data = df[ttype == "C:C" | ttype == "C:N",])
wilcox.test(f1 ~ ttype, data = df[ttype == "C:N" | ttype == "N:N",])

# # more specific metrics like prec, recall, accuracy, unique tissue types
# ggplot(df, aes(x=ttype, y=prec)) + geom_boxplot(notch = T) + 
#   theme_minimal() + ylim(0,1) + xlab("") + ylab("Precision")
# wilcox.test(prec ~ ttype, data = df[ttype == "C:C" | ttype == "C:N",])
# 
# ggplot(df, aes(x=ttype, y=rec)) + geom_boxplot(notch = T) + 
#   theme_minimal() + ylim(0,1) + xlab("") + ylab("Recall")
# wilcox.test(rec ~ ttype, data = df[ttype == "CC" | ttype == "CN",])
# 
# ggplot(df, aes(x=ttype, y=acc)) + geom_boxplot(notch = T) + 
#   theme_minimal() + ylim(0,1) + xlab("") + ylab("Accuracy")
# wilcox.test(acc ~ ttype, data = df[ttype == "C" | ttype == "CC",])
# 
# ggplot(df, aes(x=ttype, y=uth)) + geom_boxplot(notch = T) + 
#   theme_minimal() + ylim(0,34) + xlab("") + ylab("Unique Tissue Types Hit")
# wilcox.test(uth ~ ttype, data = df[ttype == "C" | ttype == "CN",], exact=T)

# improvement 1 vs 2 antigens
df[, antigens := ifelse(ttype == "C", "single", "double")]
df$antigens <- factor(df$antigens, levels=c("single", "double"))
ggplot(df, aes(x=antigens, y=f1)) + geom_boxplot(notch = T) + 
  theme_bw() + ylim(0,1) + xlab("") + ylab("F1 score") +
  theme(legend.position="bottom", axis.title = element_text(size=16), 
        axis.text= element_text(size=16), panel.grid.minor=element_blank(),
        plot.background=element_blank(), panel.background=element_blank(), panel.grid.major=element_blank())
wilcox.test(f1 ~ antigens, data=df)


# get distribution of H:H, H:L, and L:L within the top pairs
ctypes = apply(df, 1, GetClassType, dt=rnaseqmat)
df[, ctype := ctypes]
ctypes = gsub("LH", "HL", ctypes)
dfc = as.data.table(table(ctypes))
names(dfc) <- c("gate", "count")
dfc[, gate := ifelse(gate == "LH", "HL", gate)]
dfc[, percent := round((count / sum(dfc$count)), 2) * 100]

ggplot(dfc, aes(x="", y=count, fill=gate)) + 
  geom_bar(width = 1, stat = "identity") + coord_polar("y") +
  scale_fill_brewer(palette="PuBuGn") + theme_classic() + xlab("") + ylab("") +
  geom_text(aes(label = paste0(percent, "%")), position = position_stack(vjust = 0.5)) +
  theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())

# show without the singles 
dfc = dfc[(gate != "L" & gate != "H"),]
dfc[, percent := round((count / sum(dfc$count)), 2) * 100]

ggplot(dfc, aes(x="", y=count, fill=gate)) + 
  geom_bar(width = 1, stat = "identity") + coord_polar("y") +
  scale_fill_manual(values =c("#e74c3c", "#8e44ad", "#3498db")) + theme_classic() + xlab("") + ylab("") +
  geom_text(aes(label = paste0(percent, "%")), position = position_stack(vjust = 0.5)) +
  theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())


# add triples - since limited set of C:C:Cs don't enforce the maximum of 2 pairings per gene b/c that seems weird
# but still take only the top 10 per cancer - 4000 per tumor type
#tripleres = readRDS("scripts/clustering-method/rdata/all-rnaseq-combat-sketch-0.2-clin-triples.RData") #this is updated
tripleres = readRDS("scripts/clustering-method/rdata/all-rnaseq-combat-sketch-0.2-triples-S-reduction.RData")

tripleres[,can := type]
memm = data.table(pair="dummy", val="d")

tripleres[, c("g1", "g2", "g3") := tstrsplit(pair, ":", fixed=TRUE)]
tripleres[, g1lab := ifelse(g1 %in% clinSmall, "C", "N")]
tripleres[, g2lab := ifelse(g2 %in% clinSmall, "C", "N")]
tripleres[, g3lab := ifelse(g3 %in% clinSmall, "C", "N")]
tripleres[, ttype := paste(g1lab, g2lab, g3lab, sep=":")]

tripleres[,ttype := ifelse(ttype == "N:N:C", "C:N:N", ttype)]
tripleres[,ttype := ifelse(ttype == "N:C:N", "C:N:N", ttype)]
tripleres[,ttype := ifelse(ttype == "C:N:C", "C:C:N", ttype)]
tripleres[,ttype := ifelse(ttype == "N:C:C", "C:C:N", ttype)]

# triples with the S reduction should never be LLL
# therefore lets get the top k values before calculating
# their tree type
k = 100 # number of top and bottom to grab
treeTypes = c("C:C:C", "N:N:N", "C:N:N", "C:C:N")

topTriples = rbindlist(lapply(treeTypes, function(tt) {
  rbindlist(mclapply(cans, function(cn,tt)
  {
    print(cn)
    one = character()
    two = character()
    count = 0
    
    res = tripleres[can == cn,]
    res = res[ttype == tt,]
    
    tpair = character()
    tscore = numeric()
    tcan = character()
    tttype = character()
    tdb = numeric()
    tdist = numeric()
    tg1 = character()
    tg2 = character()
    tg3 = character()
    tg1lab = character()
    tg2lab = character()
    tg3lab = character()
    tctype = character()
    
    i = 1;
    sketchL = sketchLAll[can==cn]$dataset
    sketch = rnaseqmat[dataset %in% sketchL,]
    
    # pair a maximum of 2 times
    while ((count < k)  & (i < nrow(res)))
    {
      row = res[i,]
      i = i + 1;
      genes = c(row$g1, row$g2, row$g3)
      
      if (!(genes[1] %in% two) & !(genes[2] %in% two) & !(genes[3] %in% two)) # if antigen has come up twice print another one for top
      {
        if (genes[1] %in% one)
        {
          two = append(two, genes[1])
          one = one[!one %in% c(genes[1])]
        }
        else
        {
          one = append(one, genes[1])
        }
        if (genes[2] %in% one)
        {
          two = append(two, genes[2])
          one = one[!one %in% c(genes[2])]
        }
        else
        {
          one = append(one, genes[2])
        }
        if (genes[3] %in% one)
        {
          two = append(two, genes[3])
          one = one[!one %in% c(genes[3])]
        }
        else
        {
          one = append(one, genes[3])
        }
        
        gct = GetClassType2(genes, cn, sketch)
        
        if (gct != "LLL")
        {
          count = count + 1;
          tpair = append(tpair, as.character(row$pair))
          tscore = append(tscore, row$score)
          tdb = append(tdb, row$db)
          tdist = append(tdist, row$dist.man)
          tcan = append(tcan, cn)
          tttype = append(tttype, tt)
          tg1 = append(tg1, row$g1)
          tg2 = append(tg2, row$g2)
          tg3 = append(tg3, row$g3)
          tg1lab = append(tg1lab, row$g1lab)
          tg2lab = append(tg2lab, row$g2lab)
          tg3lab = append(tg3lab, row$g3lab)
          tctype = append(tctype, gct)
        }

      }
    }
    data.table(pair=tpair, score=tscore, db=tdb, dist.man=tdist, can=tcan, gene1=tg1, gene2=tg2, gene3=tg3, g1lab=tg1lab, g2lab=tg2lab, g3lab=tg3lab, ttype=tttype, ctype=tctype)
    
  }, tt, mc.cores=33))
}))

# # add label for tree type (H or L)
# fcans = rbindlist(mclapply(cans, function(cn) { 
#   sketchL = sketchLAll[can==cn]$dataset
#   sketch = rnaseqmat[dataset %in% sketchL,]
#   sketch = sketch[sketch$type == "tissue" | sketch$tissue.cancer == cn,]
#   sketch[, "target" := ifelse(sketch$tissue.cancer == cn, "target", "other")]
#   tres = topTriples[can==cn]
#   ctmp = apply(tres, 1, GetClassType3, tmp=sketch)
#   tres[, ctype := ctmp]
#   tres
# }, mc.cores=33))
# 
# topTriples = fcans
# rm(fcans)

# get the f1 scores for the top k
tmpMetrics = rbindlist(mclapply(cans, function(cann)
{
  tpc = topTriples[can==cann,]
  f1=numeric()
  prec = numeric()
  rec = numeric()
  uth = numeric()
  tcan = character()
  clist = character()
  acc = numeric()
  for (i in 1:nrow(tpc))
  {
    row = tpc[i,]
    genes = unlist(strsplit(as.character(row$pair), ":"))
    acl = assignClusterLabels(cann, genes, rnaseqmat)
    f1 = append(f1, acl$F1)
    prec = append(prec, acl$prec)
    rec = append(rec, acl$recall)
    uth = append(uth, acl$uniqTisHit)
    tcan = append(tcan, cann)
    clist = append(clist, as.character(row$pair))
    acc = append(acc, acl$accuracy)
    
  }
  data.table(pair=clist, can=tcan, f1, prec, rec, uth, acc)
}, mc.cores=33))

topTriples = merge(topTriples, tmpMetrics, by.x=c("pair", "can"), by.y=c("pair", "can"))
topTriples = topTriples[order(-f1),]
rm(tmpMetrics)

#saveRDS(topTriples,"scripts/clustering-method/rdata/topAllTriples-100.RData") # where 100 is standing in for K

#output to table for viewing
write.table(topTriples, paste0("scripts/clustering-method/res-tables/topAllTriples-", k, ".txt"), quote=F, sep="\t", row.names=F)

# quick performance plot
gg <- ggplot(topTriples, aes(x=ttype, y=f1)) + geom_boxplot(notch = T) + 
  theme_bw() + ylim(0,1) + xlab("") + ylab("F1 score") + 
  theme(legend.position="bottom", axis.title = element_text(size=12), 
        axis.text= element_text(size=12), panel.grid.minor=element_blank(),
        plot.background=element_blank(), panel.background=element_blank(), panel.grid.major=element_blank())

ggsave(paste0("scripts/clustering-method/res-figures/triples-per-type-f1.pdf"), 
       gg, width=3.5, height=4.5, units="in", device="pdf", useDingbats=FALSE)

gg <- ggplot(topTriples, aes(x=can, y=f1)) + geom_boxplot(notch = T) + 
  theme_bw() + ylim(0,1) + xlab("") + ylab("F1 score") + coord_flip() +
  theme(legend.position="bottom", axis.title = element_text(size=12), 
        axis.text= element_text(size=12), panel.grid.minor=element_blank(),
        plot.background=element_blank(), panel.background=element_blank(), panel.grid.major=element_blank())

ggsave(paste0("scripts/clustering-method/res-figures/triples-per-cancer-f1.pdf"), 
       gg, width=7, height=8.5, units="in", device="pdf", useDingbats=FALSE)




#saveRDS(topTriples,"scripts/clustering-method/rdata/topCCCTriples-10.RData") # where 10 is standing in for K

df = clinMetrics
names(df) <- c("pair", "can", "f1", "prec", "rec", "uth", "acc", "ctype", "ttype")
df[, ttype := "C"]
df = rbind(df, topPerClin2[,c("pair", "can", "f1", "prec", "rec", "uth", "acc", "ttype")], fill=T)
df = rbind(df, topTriples[,c("pair", "can", "f1", "prec", "rec", "uth", "acc", "ctype", "ttype")])

ggplot(df, aes(x=ttype, y=f1)) + geom_boxplot(notch = T) + 
  theme_bw() + ylim(0,1) + xlab("") + ylab("F1 score") + 
  theme(legend.position="bottom", axis.title = element_text(size=16), 
        axis.text= element_text(size=16), panel.grid.minor=element_blank(),
        plot.background=element_blank(), panel.background=element_blank(), panel.grid.major=element_blank())


df[, antigens := ifelse(ttype == "C", "single", "double")]
df[, antigens := ifelse(ttype == "C:C:C", "triple", antigens)]
df$antigens <- factor(df$antigens, levels=c("single", "double", "triple"))
ggplot(df, aes(x=antigens, y=f1)) + geom_boxplot(notch = T) + 
  theme_bw() + ylim(0,1) + xlab("") + ylab("F1 score") +
  theme(legend.position="bottom", axis.title = element_text(size=16), 
        axis.text= element_text(size=16), panel.grid.minor=element_blank(),
        plot.background=element_blank(), panel.background=element_blank(), panel.grid.major=element_blank())
wilcox.test(f1 ~ antigens, data=df)



# --------------------------------------
#            functions
# --------------------------------------

sketchDensityAssign <- function(pts, sketch, single, rm.outliers) # takes the gene coordinates and the density estimate and returns cluster labels
{
  if(single)
  {
    ctmp = sketch[target =="target",]
    ctmp[,1] = ifelse(as.numeric(unlist(ctmp[,1])) < 0, 0, as.numeric(unlist(ctmp[,1])))
    
    if (rm.outliers)
    {
      x = as.numeric(unlist(ctmp[,1]))
      x = ifelse(x < 0, 0, x)
      outliers = boxplot.stats(x)$out
      ctmp = ctmp[!(which(x %in% outliers)),]
    }
    cmin = min(ctmp[,1]) 
    #cmax = max(ctmp[,1])
    
    # all singles should be Hs
    labs = ifelse(as.numeric(unlist(pts[,1])) >= cmin, 1, 0)
  }
  else
  {
    ctmp = sketch[target =="target",]
    ctmp[,1] = ifelse(as.numeric(unlist(ctmp[,1])) < 0, 0, as.numeric(unlist(ctmp[,1])))
    ctmp[,2] = ifelse(as.numeric(unlist(ctmp[,2])) < 0, 0, as.numeric(unlist(ctmp[,2])))
    names(ctmp) = c("g1", "g2", "tissue.cancer", "target")
    
    if (rm.outliers)
    {
      x = as.numeric(unlist(ctmp[,1]))
      outliers = boxplot.stats(x)$out
      y = as.numeric(unlist(ctmp[,2]))
      ctmp = ctmp[!(g1 %in% outliers),]
      #ctmp = ctmp[!(which(x %in% outliers)),]
      outliers = boxplot.stats(y)$out
      #ctmp = ctmp[!(which(y %in% outliers)),]
      ctmp = ctmp[!(g2 %in% outliers),]
    }
    
    chp = chull(ctmp[,1:2])
    ch = ctmp[chp,1:2]
    # corner case if the hull only two points add 2 more to close the polygon
    if (length(ch$x) < 3)
    {
      ch = rbind(ch, c(ch[1,1] + 0.0001, ch[1,2] + 0.0001))
      ch = rbind(ch, c(ch[2,1] + 0.0001, ch[2,2] + 0.0001))
    }
    names(ch) <- c("x", "y")
    
    
    # check if it falls within the polygon
    labs = point.in.polygon(as.numeric(unlist(pts[,1])), as.numeric(unlist(pts[,2])), ch$x, ch$y, mode.checked=FALSE)
  }
  
  clusters <- ifelse(labs >= 1, "target", "other")
  return (clusters)
}

# same as sketchDensityAssign but removes data points according to the stdev
sketchDensityAssignStd <- function(pts, sketch, single, rm.outliers) # takes the gene coordinates and the density estimate and returns cluster labels
{
  stdevs = 1
  if(single)
  {
    ctmp = sketch[target =="target",]
    ctmp[,1] = ifelse(as.numeric(unlist(ctmp[,1])) < 0, 0, as.numeric(unlist(ctmp[,1])))
    
    if (rm.outliers)
    {
      x = as.numeric(unlist(ctmp[,1]))
      m = mean(x)
      s = sd(x)
      mx = m + (stdevs*s)
      mn = m - (stdevs*s)
      ctmp = ctmp[(g1 >= mn & g1 <= mx),]
    }
    cmin = min(ctmp[,1]) 
    #cmax = max(ctmp[,1])
    
    # all singles should be Hs
    labs = ifelse(as.numeric(unlist(pts[,1])) >= cmin, 1, 0)
  }
  else
  {
    ctmp = sketch[target =="target",]
    ctmp[,1] = ifelse(as.numeric(unlist(ctmp[,1])) < 0, 0, as.numeric(unlist(ctmp[,1])))
    ctmp[,2] = ifelse(as.numeric(unlist(ctmp[,2])) < 0, 0, as.numeric(unlist(ctmp[,2])))
    names(ctmp) = c("g1", "g2", "tissue.cancer", "target")
    
    if (rm.outliers)
    {
      x = as.numeric(unlist(ctmp[,1]))
      m = mean(x)
      s = sd(x)
      mx = m + (stdevs*s)
      mn = m - (stdevs*s)
      y = as.numeric(unlist(ctmp[,2]))
      ctmp = ctmp[(g1 >= mn & g1 <= mx),]
      m = mean(y)
      s = sd(y)
      mx = m + (stdevs*s)
      mn = m - (stdevs*s)
      ctmp = ctmp[(g2 >= mn & g2 <= mx),]
    }
    
    chp = chull(ctmp[,1:2])
    ch = ctmp[chp,1:2]
    # corner case if the hull only two points add 2 more to close the polygon
    if (length(ch$x) < 3)
    {
      ch = rbind(ch, c(ch[1,1] + 0.0001, ch[1,2] + 0.0001))
      ch = rbind(ch, c(ch[2,1] + 0.0001, ch[2,2] + 0.0001))
    }
    names(ch) <- c("x", "y")
    
    
    # check if it falls within the polygon
    labs = point.in.polygon(as.numeric(unlist(pts[,1])), as.numeric(unlist(pts[,2])), ch$x, ch$y, mode.checked=FALSE)
  }
  
  clusters <- ifelse(labs >= 1, "target", "other")
  return (clusters)
}

# assign based on 1D DTs
sketchDensityAssignDT <- function(pts, sketch, single, printPlot=F) # takes the gene coordinates and the density estimate and returns cluster labels
{
  
  if(single)
  {
    ctmp = sketch
    ctmp[,1] = ifelse(as.numeric(unlist(ctmp[,1])) < 0, 0, as.numeric(unlist(ctmp[,1])))
    names(ctmp) = c("g1", "tissue.cancer", "target")
    ctmp[, y := ifelse(target == "target", 1, -1)]
    ctmp$y = as.factor(ctmp$y)
    
    dt = rpart(y ~ g1, method="class", data=ctmp[,c("g1", "y")], cp=-1, maxdepth=1)
    split = dt$splits[4]
    
    names(pts) = c("g1")
    #pts$prediction = as.numeric(as.character(predict(dt, pts[,c("g1")], type="class")))
    pts$prediction = predictByRuth(pts$g1, dt$splits[4], dt$splits[2])   
  }
  else
  {
    ctmp = sketch
    if (length(names(ctmp)) == 4) # two genes
    {
      ctmp[,1] = ifelse(as.numeric(unlist(ctmp[,1])) < 0, 0, as.numeric(unlist(ctmp[,1])))
      ctmp[,2] = ifelse(as.numeric(unlist(ctmp[,2])) < 0, 0, as.numeric(unlist(ctmp[,2])))
      names(ctmp) = c("g1", "g2", "tissue.cancer", "target")
      ctmp[, y := ifelse(target == "target", 1, -1)]
      ctmp$y = as.factor(ctmp$y)
      
      names(pts) = c("g1", "g2")
      
      dt = rpart(y ~ g1, method="class", data=ctmp[,c("g1", "y")], cp=-1, maxdepth=1)
      #pts$pg1 = as.numeric(as.character(predict(dt, pts[,c("g1")], type="class")))
      pts$pg1 = predictByRuth(pts$g1, dt$splits[4], dt$splits[2])
      
      dt = rpart(y ~ g2, method="class", data=ctmp[,c("g2", "y")], cp=-1, maxdepth=1)
      #pts$pg2 = as.numeric(as.character(predict(dt, pts[,c("g2")], type="class")))
      pts$pg2 = predictByRuth(pts$g2, dt$splits[4], dt$splits[2])
      pts[, prediction := pmin(pg1, pg2)]
      
    }
    else
    {
      ctmp[,1] = ifelse(as.numeric(unlist(ctmp[,1])) < 0, 0, as.numeric(unlist(ctmp[,1])))
      ctmp[,2] = ifelse(as.numeric(unlist(ctmp[,2])) < 0, 0, as.numeric(unlist(ctmp[,2])))
      ctmp[,3] = ifelse(as.numeric(unlist(ctmp[,3])) < 0, 0, as.numeric(unlist(ctmp[,3])))
      names(ctmp) = c("g1", "g2", "g3", "tissue.cancer", "target")
      ctmp[, y := ifelse(target == "target", 1, -1)]
      ctmp$y = as.factor(ctmp$y)
      
      
      names(pts) = c("g1", "g2", "g3")
      
      dt = rpart(y ~ g1, method="class", data=ctmp[,c("g1", "y")], cp=-1, maxdepth=1)
      #pts$pg1 = as.numeric(as.character(predict(dt, pts[,c("g1")], type="class")))
      pts$pg1 = predictByRuth(pts$g1, dt$splits[4], dt$splits[2])
      
      dt = rpart(y ~ g2, method="class", data=ctmp[,c("g2", "y")], cp=-1, maxdepth=1)
      #pts$pg2 = as.numeric(as.character(predict(dt, pts[,c("g2")], type="class")))
      pts$pg2 = predictByRuth(pts$g2, dt$splits[4], dt$splits[2])
      
      dt = rpart(y ~ g3, method="class", data=ctmp[,c("g3", "y")], cp=-1, maxdepth=1)
      #pts$pg3 = as.numeric(as.character(predict(dt, pts[,c("g3")], type="class")))
      pts$pg3 = predictByRuth(pts$g3, dt$splits[4], dt$splits[2])
      
      pts[, prediction := pmin(pg1, pg2, pg3)]
      
    }
  }
#  print(table(pts$pg1))
#  print(table(pts$pg2))
#  print(table(pts$prediction))
  clusters <- ifelse(pts$prediction >= 1, "target", "other")
  return (clusters)
}

predictByRuth <- function(x, cutoff, ncat)
{
  pred = ifelse(x < cutoff, 1, -1)
  if (ncat == -1) # greater than
  {
    pred = ifelse(x > cutoff, 1, -1)
  }
  return (pred)
}


# same as sketchDensityAssign but removes data points according to the maximal margin
#tc = trainControl(method="none", number = 1)
sketchDensityAssignSVM <- function(pts, sketch, single, printPlot=F) # takes the gene coordinates and the density estimate and returns cluster labels
{
  
  if(single)
  {
    ctmp = sketch
    ctmp[,1] = ifelse(as.numeric(unlist(ctmp[,1])) < 0, 0, as.numeric(unlist(ctmp[,1])))
    names(ctmp) = c("g1", "tissue.cancer", "target")
    ctmp[, y := ifelse(target == "target", 1, -1)]
    ctmp$y = as.factor(ctmp$y)
    
    wts = 1/table(ctmp$y) # assign class weights by proportion such that the tumor class is highly weighted
    #wts["-1"] =  length(ctmp[y==1,]$tissue.cancer) / length(ctmp$tissue.cancer) 
    #wts["1"] =  length(ctmp[y==-1,]$tissue.cancer) / length(ctmp$tissue.cancer) 
    #wts["-1"] = 0.5
    #wts["1"] = 0.5
    #wts["-1"] =  1/(2 * length(ctmp[y==-1,]$tissue.cancer)) 
    #wts["1"] =  1/(2 * length(ctmp[y==1,]$tissue.cancer))
    wts["-1"] =  29/(30 * length(ctmp[y==-1,]$tissue.cancer)) 
    wts["1"] =  1/(30 * length(ctmp[y==1,]$tissue.cancer))
    
    svmfit = svm(y ~ ., data = ctmp[,c("g1", "y")], kernel="radial", class.weights=1/table(y), scale = F)
    #svmfit = svm(y ~ ., data = ctmp[,c("g1", "y")], kernel="radial", scale = F)
    
    names(pts) = c("g1")
    pts$prediction = as.numeric(predict(svmfit, pts[,c("g1")]))
    
    # if prediction is all tissue class - help it out by taking min with outliers removed
    if (length(unique(pts$prediction)) == 1)
    {
      x = as.numeric(unlist(ctmp[,1]))
      x = ifelse(x < 0, 0, x)
      outliers = boxplot.stats(x)$out
      ctmp = ctmp[!(which(x %in% outliers)),]
      cmin = min(ctmp[,1]) 
      
      # all singles should be Hs
      pts$prediction = ifelse(as.numeric(unlist(pts[,1])) >= cmin, 2, 1)
    }
    
    
  }
  else
  {
    ctmp = sketch
    if (length(names(ctmp)) == 4) # two genes
    {
      ctmp[,1] = ifelse(as.numeric(unlist(ctmp[,1])) < 0, 0, as.numeric(unlist(ctmp[,1])))
      ctmp[,2] = ifelse(as.numeric(unlist(ctmp[,2])) < 0, 0, as.numeric(unlist(ctmp[,2])))
      names(ctmp) = c("g1", "g2", "tissue.cancer", "target")
      ctmp[, y := ifelse(target == "target", 1, -1)]
      ctmp$y = as.factor(ctmp$y)
      
      wts = 1/(table(ctmp$tissue.cancer)) # assign class weights by proportion such that the tumor class is highly weighted
      #wts["-1"] =  length(ctmp[y==1,]$tissue.cancer) / length(ctmp$tissue.cancer) 
      #wts["1"] =  length(ctmp[y==-1,]$tissue.cancer) / length(ctmp$tissue.cancer) 
      #wts["-1"] =  1/(2 * length(ctmp[y==-1,]$tissue.cancer)) 
      #wts["1"] =  1/(2 * length(ctmp[y==1,]$tissue.cancer))
      #wts["-1"] =  29/(30 * length(ctmp[y==-1,]$tissue.cancer)) 
      #wts["1"] =  1/(30 * length(ctmp[y==1,]$tissue.cancer))
      #wts["-1"] = 0.5
      #wts["1"] = 0.5
      
      svmfit = svm(y ~ ., data = ctmp[,c("g1", "g2", "y")], kernel = "radial", class.weights=wts, scale = FALSE)
      #ctmp$tissue.cancer= as.factor(ctmp$tissue.cancer)
      #svmfit = svm(tissue.cancer ~ ., data = ctmp[,c("g1", "g2", "tissue.cancer")], kernel = "linear", scale = FALSE)
      
      if (printPlot)
      {
        plot(svmfit, ctmp[,c("g1", "g2", "tissue.cancer")])
      }
      
      names(pts) = c("g1", "g2")
      #pts$prediction = as.numeric(predict(svmfit, pts))
      pts$prediction = as.numeric(predict(dt, pts, type="class"))
      #pts$prediction = as.character(predict(svmfit, pts))
      #print(table(ctmp$tissue.cancer, pts$prediction))
      #pts[,prediction := ifelse(prediction == target, 2, 0)]
    }
    else
    {
      ctmp[,1] = ifelse(as.numeric(unlist(ctmp[,1])) < 0, 0, as.numeric(unlist(ctmp[,1])))
      ctmp[,2] = ifelse(as.numeric(unlist(ctmp[,2])) < 0, 0, as.numeric(unlist(ctmp[,2])))
      ctmp[,3] = ifelse(as.numeric(unlist(ctmp[,3])) < 0, 0, as.numeric(unlist(ctmp[,3])))
      names(ctmp) = c("g1", "g2", "g3", "tissue.cancer", "target")
      ctmp[, y := ifelse(target == "target", 1, -1)]
      ctmp$y = as.factor(ctmp$y)
      
      wts = 1/table(ctmp$y) # assign class weights by proportion such that the tumor class is highly weighted
      #wts["-1"] =  length(ctmp[y==1,]$tissue.cancer) / length(ctmp$tissue.cancer) 
      #wts["1"] =  length(ctmp[y==-1,]$tissue.cancer) / length(ctmp$tissue.cancer) 
      wts["-1"] =  29/(30 * length(ctmp[y==-1,]$tissue.cancer)) 
      wts["1"] =  1/(30 * length(ctmp[y==1,]$tissue.cancer))
      #wts["-1"] = 0.5
      #wts["1"] = 0.5
      
      svmfit = svm(y ~ ., data = ctmp[,c("g1", "g2","g3", "y")], kernel = "radial", class.weights=wts, scale = FALSE)
      #svmfit = svm(y ~ ., data = ctmp[,c("g1", "g2","g3", "y")], kernel = "radial", scale = FALSE)
      
      
      names(pts) = c("g1", "g2", "g3")
      pts$prediction = as.numeric(predict(svmfit, pts))
    }
  }
  
  clusters <- ifelse(pts$prediction > 1, "target", "other")
  return (clusters)
}


assignClusterLabels <- function(cn, genes, dt, pp=F)
{
  # get the sketch subset
  #sketchL = readLines(paste0("scripts/clustering-method/sketches/", gsub(" ", "-", tolower(can)), "-sketch.txt"))
  #sketchL = fread(paste0("scripts/clustering-method/sketches/", gsub(" ", "-", tolower(can)), "-sketch.txt"), header=F)
  #sketchL = sketchL$V1
  sketchL = sketchLAll[can==cn]$dataset
  sketch = dt[dataset %in% sketchL,]
  tmp = dt[!(dataset %in% sketchL),]
  
  if (length(genes) == 2)
  {
    sketch = sketch[sketch$type == "tissue" | sketch$tissue.cancer == cn,]
    sketch = sketch[,c(genes[1], genes[2], "tissue.cancer"), with =F]
    sketch[, "target" := ifelse(sketch$tissue.cancer == cn, "target", "other")]
    
    tmp = tmp[tmp$type == "tissue" | tmp$tissue.cancer == cn,]
    tmp = tmp[,c(genes[1], genes[2], "tissue.cancer"), with=F]
    tmp[, "target" := ifelse(tmp$tissue.cancer == cn, "target", "other")]
    tmp[,1] = ifelse(as.numeric(unlist(tmp[,1])) < 0, 0, as.numeric(unlist(tmp[,1])))
    tmp[,2] = ifelse(as.numeric(unlist(tmp[,2])) < 0, 0, as.numeric(unlist(tmp[,2])))
    
    #tmp[, clustlab := sketchDensityAssign(tmp[,1:2, with=F], sketch, F, T)]
    #tmp[, clustlab := sketchDensityAssignStd(tmp[,1:2, with=F], sketch, F, T)]
    tmp[, clustlab := sketchDensityAssignDT(tmp[,1:2, with=F], sketch, F, pp)]
    
    # gg1 = ggplot(tmp, aes_string(x=names(tmp)[1], y=names(tmp)[2])) + geom_point(aes(color=target), alpha = 0.4) + #facet_wrap(~target) +
    #   scale_color_manual(values=c("#34495e", "#e74c3c", "#bdc3c7")) +
    #   theme_minimal() + theme(legend.position="bottom", plot.title = element_text(size=12))
    # 
    # gg2 = ggplot(tmp, aes_string(x=names(tmp)[1], y=names(tmp)[2])) + geom_point(aes(color=clustlab), alpha = 0.4) + #facet_wrap(~target) +
    #   scale_color_manual(values=c("#34495e", "#e74c3c", "#bdc3c7")) +
    #   theme_minimal() + theme(legend.position="bottom", plot.title = element_text(size=12))
    # 
    # plt = grid.arrange(gg1, gg2, ncol = 2)
    # plot(plt)
    
    tp = length(tmp[target == "target" & clustlab == "target",]$clustlab)
    tn = length(tmp[target == "other" & clustlab == "other",]$clustlab)
    fp = length(tmp[target == "other" & clustlab == "target",]$clustlab)
    #precision = tp / (tp + fp)
    prec = tp / (tp + fp)
    #recall = tp / total number of positives (fn) 
    rec = tp / length(tmp[target=="target",]$clustlab)
    f1 = (2 * prec * rec) / (prec + rec)
    #Accuracy = (TP+TN)/(TP+TN+FP+FN)
    acc = (tp + tn) / length(tmp$target)
    
    uniqTisHit = length(unique(tmp[clustlab == "target" & target != "target",]$tissue.cancer))
    
    return(list(F1=f1, prec=prec, recall=rec, uniqTisHit=uniqTisHit, accuracy=acc))
  }
  if (length(genes) == 3)
  {
    sketch = sketch[sketch$type == "tissue" | sketch$tissue.cancer == cn,]
    sketch = sketch[,c(genes[1], genes[2], genes[3], "tissue.cancer"), with =F]
    sketch[, "target" := ifelse(sketch$tissue.cancer == cn, "target", "other")]
    
    tmp = tmp[tmp$type == "tissue" | tmp$tissue.cancer == cn,]
    tmp = tmp[,c(genes[1], genes[2], genes[3], "tissue.cancer"), with=F]
    tmp[, "target" := ifelse(tmp$tissue.cancer == cn, "target", "other")]
    tmp[,1] = ifelse(as.numeric(unlist(tmp[,1])) < 0, 0, as.numeric(unlist(tmp[,1])))
    tmp[,2] = ifelse(as.numeric(unlist(tmp[,2])) < 0, 0, as.numeric(unlist(tmp[,2])))
    tmp[,3] = ifelse(as.numeric(unlist(tmp[,3])) < 0, 0, as.numeric(unlist(tmp[,3])))
    
    tmp[, clustlab := sketchDensityAssignDT(tmp[,1:3, with=F], sketch, F)]
    
    tp = length(tmp[target == "target" & clustlab == "target",]$clustlab)
    tn = length(tmp[target == "other" & clustlab == "other",]$clustlab)
    fp = length(tmp[target == "other" & clustlab == "target",]$clustlab)
    #precision = tp / (tp + fp)
    prec = tp / (tp + fp)
    #recall = tp / total number of positives (fn) 
    rec = tp / length(tmp[target=="target",]$clustlab)
    f1 = (2 * prec * rec) / (prec + rec)
    #Accuracy = (TP+TN)/(TP+TN+FP+FN)
    acc = (tp + tn) / length(tmp$target)
    
    uniqTisHit = length(unique(tmp[clustlab == "target" & target != "target",]$tissue.cancer))
    
    return(list(F1=f1, prec=prec, recall=rec, uniqTisHit=uniqTisHit, accuracy=acc))
  }
  
  if(length(genes) == 1) # single gene case
  {
    sketch = sketch[sketch$type == "tissue" | sketch$tissue.cancer == cn,]
    sketch = sketch[,c(genes[1], "tissue.cancer"), with =F]
    sketch[, "target" := ifelse(sketch$tissue.cancer == cn, "target", "other")]
    
    tmp = tmp[tmp$type == "tissue" | tmp$tissue.cancer == cn,]
    tmp = tmp[,c(genes[1], "tissue.cancer"), with=F]
    tmp[, "target" := ifelse(tmp$tissue.cancer == cn, "target", "other")]
    tmp[,1] = ifelse(as.numeric(unlist(tmp[,1])) < 0, 0, as.numeric(unlist(tmp[,1])))
    
    #tmp[, clustlab := sketchDensityAssign(tmp[,1, with=F], sketch, T, T)]
    #tmp[, clustlab := sketchDensityAssignStd(tmp[,1, with=F], sketch, T, T)]
    tmp[, clustlab := sketchDensityAssignDT(tmp[,1, with=F], sketch, T)]
    
    # gg1 = ggplot(tmp, aes_string(x=names(tmp)[1], y=1)) + geom_point(aes(color=target), alpha = 0.4) + #facet_wrap(~target) +
    #   scale_color_manual(values=c("#34495e", "#e74c3c", "#bdc3c7")) +
    #   theme_minimal() + theme(legend.position="bottom", plot.title = element_text(size=12))
    # 
    # gg2 = ggplot(tmp, aes_string(x=names(tmp)[1], y=1)) + geom_point(aes(color=clustlab), alpha = 0.4) + #facet_wrap(~target) +
    #   scale_color_manual(values=c("#34495e", "#e74c3c", "#bdc3c7")) +
    #   theme_minimal() + theme(legend.position="bottom", plot.title = element_text(size=12))
    # 
    # plt = grid.arrange(gg1, gg2, ncol = 2)
    
    tp = length(tmp[target == "target" & clustlab == "target",]$clustlab)
    tn = length(tmp[target == "other" & clustlab == "other",]$clustlab)
    fp = length(tmp[target == "other" & clustlab == "target",]$clustlab)
    #precision = tp / (tp + fp)
    prec = tp / (tp + fp)
    #recall = tp / total number of positives (fn) 
    rec = tp / length(tmp[target=="target",]$clustlab)
    f1 = (2 * prec * rec) / (prec + rec)
    #Accuracy = (TP+TN)/(TP+TN+FP+FN)
    acc = (tp + tn) / length(tmp$target)
    
    uniqTisHit = length(unique(tmp[clustlab == "target" & target != "target",]$tissue.cancer))
    
    return(list(F1=f1, prec=prec, recall=rec, uniqTisHit=uniqTisHit, accuracy=acc))
  }
  
}


# return gate type (LH, HH, HL, LL, L, or H) 
# corresponding to where the boundary of the pair is drawn
GetClassType <- function(row, dt)
{
  genes = unlist(strsplit(row["pair"], ":"))
  cn = row["can"]
  sketchL = sketchLAll[can==cn]$dataset
  tmp = dt[dataset %in% sketchL,]
  
  if (length(genes) == 1) # single gene case
  {
    tmp = tmp[tmp$type == "tissue" | tmp$tissue.cancer == cn,]
    tmp = tmp[,c(genes[1], "tissue.cancer"), with=F]
    tmp[, "target" := ifelse(tmp$tissue.cancer == cn, "target", "other")]
    
    tmp[,1] = ifelse(as.numeric(unlist(tmp[,1])) < 0, 0, as.numeric(unlist(tmp[,1])))
    
    cm = colMeans(tmp[target == "other",1, with=F])
    cm = rbind(cm, colMeans(tmp[target == "target",1, with=F]))
    row.names(cm) = NULL
    cm = as.matrix(cm)
    fl = "H"
    if (as.numeric(cm[1,]) > as.numeric(cm[2,]))
    {
      fl = "L"
    }
    #memm = rbind(memm, data.table(pair=paste0(genes[1], "-", can), val=fl))
    return(fl)
  }
  if (length(genes) == 2) # pair of genes
  {
    tmp = tmp[tmp$type == "tissue" | tmp$tissue.cancer == cn,]
    tmp = tmp[,c(genes[1], genes[2], "tissue.cancer"), with=F]
    tmp[, "target" := ifelse(tmp$tissue.cancer == cn, "target", "other")]
    
    tmp[,1] = ifelse(as.numeric(unlist(tmp[,1])) < 0, 0, as.numeric(unlist(tmp[,1])))
    tmp[,2] = ifelse(as.numeric(unlist(tmp[,2])) < 0, 0, as.numeric(unlist(tmp[,2])))
    
    cm = colMeans(tmp[target == "other",1:2, with=F])
    cm = rbind(cm, colMeans(tmp[target == "target",1:2, with=F]))
    row.names(cm) = c("other", "target")
    cm = as.matrix(cm)
    fl = "H"
    sl = "H"
    if (as.numeric(cm[1,1]) > as.numeric(cm[2,1]))
    {
      fl = "L"
    }
    if (as.numeric(cm[1,2]) > as.numeric(cm[2,2]))
    {
      sl = "L"
    }
    # ggplot(tmp, aes_string(x=names(tmp)[1], y=names(tmp)[2])) + geom_point(aes(color=target), alpha = 0.4) + #facet_wrap(~target) +
    #   scale_color_manual(values=c("#34495e", "#e74c3c", "#bdc3c7")) +
    #   geom_point(data = as.data.frame(cm), size = 3, shape=19,  color="purple") +
    #   theme_minimal() + theme(legend.position="bottom", plot.title = element_text(size=12))
    #memm = rbind(memm, data.table(pair=paste0(genes[1], "-", cn), val=fl))
    #memm = rbind(memm, data.table(pair=paste0(genes[2], "-", cn), val=sl))
    return (paste0(fl, sl))
  }
  if (length(genes) > 2) # pair of genes
  {
    if((paste0(genes[1], "-", cn) %in% memm$pair) & (paste0(genes[2], "-", cn) %in% memm$pair) & (paste0(genes[3], "-", cn) %in% memm$pair))
    {
      return (paste0(memm[pair==paste0(genes[1], "-", cn)]$val, memm[pair==paste0(genes[2], "-", cn)]$val, memm[pair==paste0(genes[3], "-", cn)]$val))
    }
    
    tmp = tmp[tmp$type == "tissue" | tmp$tissue.cancer == cn,]
    tmp = tmp[,c(genes[1], genes[2], genes[3], "tissue.cancer"), with=F]
    tmp[, "target" := ifelse(tmp$tissue.cancer == cn, "target", "other")]
    
    tmp[,1] = ifelse(as.numeric(unlist(tmp[,1])) < 0, 0, as.numeric(unlist(tmp[,1])))
    tmp[,2] = ifelse(as.numeric(unlist(tmp[,2])) < 0, 0, as.numeric(unlist(tmp[,2])))
    tmp[,3] = ifelse(as.numeric(unlist(tmp[,3])) < 0, 0, as.numeric(unlist(tmp[,3])))
    
    cm = colMeans(tmp[target == "other",1:3, with=F])
    cm = rbind(cm, colMeans(tmp[target == "target",1:3, with=F]))
    row.names(cm) = c("other", "target")
    cm = as.matrix(cm)
    fl = "H"
    sl = "H"
    en = "H"
    if (as.numeric(cm[1,1]) > as.numeric(cm[2,1]))
    {
      fl = "L"
    }
    if (as.numeric(cm[1,2]) > as.numeric(cm[2,2]))
    {
      sl = "L"
    }
    if (as.numeric(cm[1,3]) > as.numeric(cm[2,3]))
    {
      en = "L"
    }
    #memm = rbind(memm, data.table(pair=paste0(genes[1], "-", cn), val=fl))
    #memm = rbind(memm, data.table(pair=paste0(genes[2], "-", cn), val=sl))
    #memm = rbind(memm, data.table(pair=paste0(genes[3], "-", cn), val=en))
    return (paste0(fl, sl, en))
  }
}

# return gate type (LH, HH, HL, LL, L, or H) 
# corresponding to where the boundary of the pair is drawn
GetClassType2 <- function(genes, can, tmp)
{
  #sketchL = readLines(paste0("scripts/clustering-method/sketches/", gsub(" ", "-", tolower(can)), "-sketch.txt"))
  
  #sketchL = sketchLAll[can==can]$dataset # THIS IS BAD
  #tmp = tmp[dataset %in% sketchL,]
  
  if (length(genes) == 1) # single gene case
  {
    tmp = tmp[tmp$type == "tissue" | tmp$tissue.cancer == can,]
    tmp = tmp[,c(genes[1], "tissue.cancer"), with=F]
    tmp[, "target" := ifelse(tmp$tissue.cancer == can, "target", "other")]
    
    tmp[,1] = ifelse(as.numeric(unlist(tmp[,1])) < 0, 0, as.numeric(unlist(tmp[,1])))
    
    cm = colMeans(tmp[target == "other",1, with=F])
    cm = rbind(cm, colMeans(tmp[target == "target",1, with=F]))
    row.names(cm) = NULL
    cm = as.matrix(cm)
    fl = "H"
    if (as.numeric(cm[1,]) > as.numeric(cm[2,]))
    {
      fl = "L"
    }
    return(fl)
  }
  if (length(genes) == 2) # pair of genes
  {
    tmp = tmp[tmp$type == "tissue" | tmp$tissue.cancer == can,]
    tmp = tmp[,c(genes[1], genes[2], "tissue.cancer"), with=F]
    tmp[, "target" := ifelse(tmp$tissue.cancer == can, "target", "other")]
    
    tmp[,1] = ifelse(as.numeric(unlist(tmp[,1])) < 0, 0, as.numeric(unlist(tmp[,1])))
    tmp[,2] = ifelse(as.numeric(unlist(tmp[,2])) < 0, 0, as.numeric(unlist(tmp[,2])))
    
    cm = colMeans(tmp[target == "other",1:2, with=F])
    cm = rbind(cm, colMeans(tmp[target == "target",1:2, with=F]))
    row.names(cm) = c("other", "target")
    cm = as.matrix(cm)
    fl = "H"
    sl = "H"
    if (as.numeric(cm[1,1]) > as.numeric(cm[2,1]))
    {
      fl = "L"
    }
    if (as.numeric(cm[1,2]) > as.numeric(cm[2,2]))
    {
      sl = "L"
    }
    # ggplot(tmp, aes_string(x=names(tmp)[1], y=names(tmp)[2])) + geom_point(aes(color=target), alpha = 0.4) + #facet_wrap(~target) +
    #   scale_color_manual(values=c("#34495e", "#e74c3c", "#bdc3c7")) +
    #   geom_point(data = as.data.frame(cm), size = 3, shape=19,  color="purple") +
    #   theme_minimal() + theme(legend.position="bottom", plot.title = element_text(size=12))
    return (paste0(fl, sl))
  }
  if (length(genes) > 2)
  {
    tmp = tmp[tmp$type == "tissue" | tmp$tissue.cancer == can,]
    tmp = tmp[,c(genes[1], genes[2], genes[3], "tissue.cancer"), with=F]
    tmp[, "target" := ifelse(tmp$tissue.cancer == can, "target", "other")]
    
    tmp[,1] = ifelse(as.numeric(unlist(tmp[,1])) < 0, 0, as.numeric(unlist(tmp[,1])))
    tmp[,2] = ifelse(as.numeric(unlist(tmp[,2])) < 0, 0, as.numeric(unlist(tmp[,2])))
    tmp[,3] = ifelse(as.numeric(unlist(tmp[,3])) < 0, 0, as.numeric(unlist(tmp[,3])))
    
    cm = colMeans(tmp[target == "other",1:3, with=F])
    cm = rbind(cm, colMeans(tmp[target == "target",1:3, with=F]))
    row.names(cm) = c("other", "target")
    cm = as.matrix(cm)
    fl = "H"
    sl = "H"
    en = "H"
    if (as.numeric(cm[1,1]) > as.numeric(cm[2,1]))
    {
      fl = "L"
    }
    if (as.numeric(cm[1,2]) > as.numeric(cm[2,2]))
    {
      sl = "L"
    }
    if (as.numeric(cm[1,3]) > as.numeric(cm[2,3]))
    {
      en = "L"
    }
    return (paste0(fl, sl, en))
  }
  
}

# optimized function just for 3
GetClassType3 <- function(row, tmp)
{
  genes = unlist(strsplit(row["pair"], ":"))
  
  tmp = tmp[,c(genes[1], genes[2], genes[3], "tissue.cancer", "target"), with=F]
  tmp[,1] = ifelse(as.numeric(unlist(tmp[,1])) < 0, 0, as.numeric(unlist(tmp[,1])))
  tmp[,2] = ifelse(as.numeric(unlist(tmp[,2])) < 0, 0, as.numeric(unlist(tmp[,2])))
  tmp[,3] = ifelse(as.numeric(unlist(tmp[,3])) < 0, 0, as.numeric(unlist(tmp[,3])))
  
  cm = colMeans(tmp[target == "other",1:3, with=F])
  cm = rbind(cm, colMeans(tmp[target == "target",1:3, with=F]))
  row.names(cm) = c("other", "target")
  cm = as.matrix(cm)
  fl = "H"
  sl = "H"
  en = "H"
  if (as.numeric(cm[1,1]) > as.numeric(cm[2,1]))
  {
    fl = "L"
  }
  if (as.numeric(cm[1,2]) > as.numeric(cm[2,2]))
  {
    sl = "L"
  }
  if (as.numeric(cm[1,3]) > as.numeric(cm[2,3]))
  {
    en = "L"
  }

  return (paste0(fl, sl, en))

}


GetClassType4 <- function(row, tmp)
{
  genes = unlist(strsplit(row["pair"], ":"))
  
  tmp = tmp[,c(genes[1], genes[2], "tissue.cancer", "target"), with=F]
  tmp[,1] = ifelse(as.numeric(unlist(tmp[,1])) < 0, 0, as.numeric(unlist(tmp[,1])))
  tmp[,2] = ifelse(as.numeric(unlist(tmp[,2])) < 0, 0, as.numeric(unlist(tmp[,2])))
  
  cm = colMeans(tmp[target == "other",1:2, with=F])
  cm = rbind(cm, colMeans(tmp[target == "target",1:2, with=F]))
  row.names(cm) = c("other", "target")
  cm = as.matrix(cm)
  fl = "H"
  sl = "H"
  if (as.numeric(cm[1,1]) > as.numeric(cm[2,1]))
  {
    fl = "L"
  }
  if (as.numeric(cm[1,2]) > as.numeric(cm[2,2]))
  {
    sl = "L"
  }
  
  return (paste0(fl, sl))
  
}


# get the pair type and calculate the manhattan distance in one go


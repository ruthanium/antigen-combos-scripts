library("ggnet")
library("sna")
library(network)

setwd("/Genomics/ogtr04/rd6/BioGPS/")
setwd("/home/rd6/BioGPS/") # for lumos

tm = fread("data/transmembrane-markers/filteredTransmembraneList-updated-12-19.txt", header=F)
rnaseqmat = readRDS("scripts/clustering-method/rdata/combat-rnaseq-tm-clin-backcompatible.RData")
cans = unique(rnaseqmat[rnaseqmat$type == "tumor",]$tissue.cancer)
clinSmall = fread("data/transmembrane-markers/clin-paper-dec-19.txt", header=F)
clinSmall = clinSmall$V1
clinSmall = clinSmall[clinSmall %in% names(rnaseqmat)]
tm = tm$V1
tm = tm[tm %in% names(rnaseqmat)]

tripleres = readRDS("scripts/clustering-method/rdata/all-rnaseq-combat-sketch-0.2-triples-S-reduction.RData")
singles = readRDS("scripts/clustering-method/rdata/singles-complete.RData")
#doubles = readRDS("scripts/clustering-method/rdata/doubles-topPerClin.RData") # missing doubles
# recreate top doubles for plot - can now just load the RDS file see below
# bigres = readRDS("scripts/clustering-method/rdata/bigres-no-LL-for-sara.RData")
# bigres = bigres[order(-score),]
# bigres[, ttype := ifelse(ttype == "N:C", "C:N", ttype)]
# 
# # build the doubles list again - trying to preserve more of the pair types
# lim = 500
# d2s = rbindlist(mclapply(cans, function(cn)
# {
#   print (cn);
#   
#   res = bigres[type == cn,]
#   r1 = res[ttype == "N:N",]
#   # if (nrow(r1) >= lim)
#   # {
#   #   r1 = r1[1:lim]
#   # }
#   r2 = res[ttype == "C:N",]
#   # if (nrow(r2) >= lim)
#   # {
#   #   r2 = r2[1:lim]
#   # }
#   r3 = res[ttype == "C:C",]
#   # if (nrow(r3) >= lim)
#   # {
#   #   r3 = r3[1:lim]
#   # }
#   
#   res = rbind(r1, r2, r3)
#   rm(r1, r2, r3)
#   
#   # only allow gene pairing 2 times
#   keep = character()
#   for (gt in unique(res$ttype))
#   {
#     one = character()
#     two = character()
#     count = 1;
#     for (p in res[ttype==gt,]$pair)
#     {
#       if (count > lim)
#       {
#         break;
#       }
#       ps = unlist(strsplit(p, ":"))
#       genes = c(ps[1], ps[2])
#       rm(ps)
#       if (!(genes[1] %in% two) & !(genes[2] %in% two)) # if antigen has come up twice print another one for top
#       {
#         if (genes[1] %in% one)
#         {
#           two = append(two, genes[1])
#           one = one[!one %in% c(genes[1])]
#         }
#         else
#         {          
#           one = append(one, genes[1])
#         }
#         if (genes[2] %in% one)
#         {
#           two = append(two, genes[2])
#           one = one[!one %in% c(genes[2])]
#         }
#         else
#         {
#           one = append(one, genes[2])
#         }
#         keep = append(keep, as.character(p))
#         count = count + 1;
#       }
#     }
#   }
#   res[pair %in% keep,]
# }, mc.cores=25))

#saveRDS(d2s, "scripts/clustering-method/rdata/doubles-500-from-bigres.RData")
d2s = readRDS("scripts/clustering-method/rdata/doubles-500-from-bigres.RData")
d2s = d2s[, head(.SD, 50), by=c("type", "ttype")] # limit to max of top 50
names(d2s)[names(d2s) =='type'] <- 'can'

singles = singles[!is.na(db),]
singles[,dbscaled := 1-((log(db) - min(log(db)))/(max(log(db)) - min(log(db))))] # make score between 0 and 1
singles[,distscaled := (log(dist.man) - min(log(dist.man)))/(max(log(dist.man)) - min(log(dist.man)))]
singles[,score := apply(singles[,c("dbscaled", "distscaled"), with=F], 1, min)]
singles = singles[ttype != "L",]
singles = singles[order(-score),]

topsingles = singles[, head(.SD, 50), by=c("can", "glab")]

# the f1 for the doubles is incorrect!!!
d2 = rbindlist(mclapply(cans, function(cann)
{
  tpc = d2s[can==cann,]
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
    clist = append(clist, as.character(row$pair))
    acc = append(acc, acl$accuracy)
    
  }
  data.table(pair=clist, can=tcan, f1, prec, rec, uth, acc)
}, mc.cores=25))


dbls2 = merge(d2s, d2, by=c("pair", "can"))
#doubles2 = merge(doubles[,c("pair", "can", "db", "dist.man", "dbscaled", "distscaled", "score", "gene1", "gene2", "g1lab", "g2lab", "ttype", "ctype", "score2")], d2, by=c("pair", "can"))
#doubles2[,score2 := apply(doubles2[,c("score", "f1"), with=F], 1, mean)]
#saveRDS(dbls2, "scripts/clustering-method/rdata/doubles-topPerClin-updated-040720.RData")
dbls2[,score2 := apply(dbls2[,c("score", "f1"), with=F], 1, mean)]

# f1 for singles is incorrect!!
s2 = rbindlist(mclapply(cans, function(cann)
{
  tpc = topsingles[can==cann,]
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
    acl = assignClusterLabels(cann, as.character(row$gene), rnaseqmat)
    f1 = append(f1, acl$F1)
    prec = append(prec, acl$prec)
    rec = append(rec, acl$recall)
    uth = append(uth, acl$uniqTisHit)
    tcan = append(tcan, cann)
    clist = append(clist, as.character(row$gene))
    acc = append(acc, acl$accuracy)
    
  }
  data.table(pair=clist, can=tcan, f1, prec, rec, uth, acc)
}, mc.cores=33))

singles2 = merge(topsingles[,c("gene", "can", "db", "dist.man", "dbscaled", "distscaled", "score", "glab", "ttype")], s2, by.x=c("gene", "can"), by.y=c("pair", "can"))
saveRDS(singles2, "scripts/clustering-method/rdata/topSingles.RData")
singles2[,score2 := apply(singles2[,c("score", "f1"), with=F], 1, mean)]
singles2 = singles2[order(-score),]


toptriples = fread("scripts/clustering-method/res-tables/topAllTriples-100.txt")

# make the pie chart
dfc = as.data.table(table(toptriples$ctype))
names(dfc) <- c("gate", "count")
dfc[, gate := ifelse(gate == "HLH", "HHL", gate)]
dfc[, gate := ifelse(gate == "LHH", "HHL", gate)]
dfc[, gate := ifelse(gate == "LHL", "HLL", gate)]
dfc[, gate := ifelse(gate == "LLH", "HLL", gate)]
dfc = dfc[, .(count = sum(count)), by = gate]
dfc[, percent := round((count / sum(dfc$count)), 2) * 100]

ggplot(dfc, aes(x="", y=count, fill=gate)) + 
  geom_bar(width = 1, stat = "identity") + coord_polar("y") +
  scale_fill_manual(values =c("#e74c3c", "#8e44ad", "#3498db")) + theme_classic() + xlab("") + ylab("") +
  geom_text(aes(label = paste0(percent, "%")), position = position_stack(vjust = 0.5)) +
  theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())

# make plot showing C, CN, CCN or CNN
# f1 plot, precision plot, recall plot
# top 10 per gate type per cancer
# get top 10 singles
k = 10
singles2 = singles2[order(-score),]
df = singles2[, head(.SD, k), by=c("can", "glab")]
df = singles2
df[, category := "single", ]
df[, score2 := NULL]
names(df)[names(df) =='ttype'] <- 'ctype'
names(df)[names(df) =='glab'] <- 'ttype'
names(df)[names(df) =='gene'] <- 'pair'

dbls2 = dbls2[order(-score),]
df2 = dbls2[, head(.SD, k), by=c("can", "ttype")]
df2[, category := "double"]

toptriples = toptriples[order(-score),]
df3 = toptriples[, head(.SD, k), by=c("can", "ttype")]
df3[, category := "triple"]
df3

cols = c("pair", "can", "score", "f1", "prec", "rec", "db", "dist.man", "ttype", "ctype", "category")
df.all = rbind(df[,cols, with=F], df2[,cols, with=F], df3[,cols, with=F])
df = df.all
rm(df2, df3, df.all)
df$category = factor(df$category, levels=c("single", "double", "triple"))

gg <- ggplot(df, aes(x=category, y=f1)) + geom_boxplot(outlier.shape = NA, notch=T) + 
  theme_bw() + ylim(0,1) + xlab("") + ylab("F1 score") + 
  theme(legend.position="bottom", axis.title = element_text(size=12), 
        axis.text= element_text(size=12), panel.grid.minor=element_blank(),
        plot.background=element_blank(), panel.background=element_blank(), panel.grid.major=element_blank())
gg
ggsave(paste0("scripts/clustering-method/res-figures/single-double-triple-f1.pdf"),
       gg, width=3.5, height=4.5, units="in", device="pdf", useDingbats=FALSE)

wilcox.test(f1 ~ category, data=df[category %in% c("single", "double")])$p.value
wilcox.test(f1 ~ category, data=df[category %in% c("double", "triple")])$p.value


tmp = df[ttype %in% c("C", "C:N", "C:C:N"),]
tmp$ttype = factor(tmp$ttype, levels=c("C", "C:N", "C:C:N"))
gg <- ggplot(tmp, aes(x=ttype, y=f1)) + geom_boxplot(outlier.shape = NA, notch=T) + 
  theme_bw() + ylim(0,1) + xlab("") + ylab("F1 score") + 
  theme(legend.position="bottom", axis.title = element_text(size=12), 
        axis.text= element_text(size=12), panel.grid.minor=element_blank(),
        plot.background=element_blank(), panel.background=element_blank(), panel.grid.major=element_blank())
gg
ggsave(paste0("scripts/clustering-method/res-figures/c-cn-ccn.pdf"),
       gg, width=3.5, height=4.5, units="in", device="pdf", useDingbats=FALSE)

wilcox.test(f1 ~ ttype, data=tmp[ttype %in% c("C", "C:N"),])$p.value
wilcox.test(f1 ~ ttype, data=tmp[ttype %in% c("C:N", "C:C:N"),])$p.value

tmp = df[ttype %in% c("C", "N", "C:C", "C:N", "N:N"),]
tmp$ttype = factor(tmp$ttype, levels=c("C", "N", "C:C", "C:N", "N:N"))
gg <- ggplot(tmp, aes(x=ttype, y=f1)) + geom_boxplot(outlier.shape = NA, notch=T) + 
  theme_bw() + ylim(0,1) + xlab("") + ylab("F1 score") + 
  theme(legend.position="bottom", axis.title = element_text(size=12), 
        axis.text= element_text(size=12), panel.grid.minor=element_blank(),
        plot.background=element_blank(), panel.background=element_blank(), panel.grid.major=element_blank())
gg
ggsave(paste0("scripts/clustering-method/res-figures/c-n-cc-cn-nn.pdf"),
       gg, width=4, height=4.5, units="in", device="pdf", useDingbats=FALSE)

wilcox.test(f1 ~ ttype, data=tmp[ttype %in% c("C", "N"),])$p.value
wilcox.test(f1 ~ ttype, data=tmp[ttype %in% c("N", "C:C"),])$p.value
wilcox.test(f1 ~ ttype, data=tmp[ttype %in% c("C:C", "C:N"),])$p.value
wilcox.test(f1 ~ ttype, data=tmp[ttype %in% c("C:N", "N:N"),])$p.value


# precision vs recall
ggplot(df, aes(x=prec, y=rec, color=category)) + geom_point() + 
  theme_bw() + ylim(0,1) +
  theme(legend.position="bottom", axis.title = element_text(size=12), 
        axis.text= element_text(size=12), panel.grid.minor=element_blank(),
        plot.background=element_blank(), panel.background=element_blank(), panel.grid.major=element_blank())

df2 = melt(df[,c("pair", "can", "category", "prec", "rec")], id.vars=c("pair", "can", "category"))

carrots = ggplot(df2, aes(x=category, y=value, color=variable)) + geom_quasirandom(position ="dodge", alpha=0.8) + theme_bw() + facet_wrap(~variable) +
  scale_color_manual(values=c("#27ae60", "#e67e22")) + 
  #stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 0.5, color="black")
carrots
ggsave(paste0("scripts/clustering-method/res-figures/prec-rec-carrot.pdf"),
       carrots, width=6, height=4, units="in", device="pdf", useDingbats=FALSE)

gg = ggplot(df2, aes(x=category, y=value, color=variable)) + geom_quasirandom(position ="dodge", alpha=0.8) + theme_bw() + facet_wrap(~variable) +
  scale_color_manual(values=c("#999999", "#999999")) + ylab("") + xlab("") + 
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 0.5, color="black")
gg
ggsave(paste0("scripts/clustering-method/res-figures/prec-rec-compare.pdf"),
       gg, width=6, height=4, units="in", device="pdf", useDingbats=FALSE)


# lollipopsss
# just show top and top
df1 = df[category == "single" & ttype=="C",]
df2 = df[category == "double",]
df3 = df[category == "triple",]

maxS = df1[df1[, .I[which.max(f1)], by=can]$V1]
maxD = df2[df2[, .I[which.max(f1)], by=can]$V1]
maxT = df3[df3[, .I[which.max(f1)], by=can]$V1]
rm(df1, df2, df3)
maxPairs = merge(maxS, maxD, by="can")
maxPairs = merge(maxPairs, maxT, by="can")

names(maxPairs)[names(maxPairs) =='pair.x'] <- 'single'
names(maxPairs)[names(maxPairs) =='pair.y'] <- 'double'
names(maxPairs)[names(maxPairs) =='pair'] <- 'triple'
maxPairs[,category.x := NULL]
maxPairs[,category.y := NULL]
maxPairs[,category := NULL]

# simple just show best C pair to best double
maxPairs$can <- factor(maxPairs$can, levels=rev(unique(maxPairs[order(-f1.y), ]$can)))
ggplot(maxPairs) + geom_segment(aes(x=can, xend=can, y=f1.x, yend=f1.y), size=1.5, color="#dadfe3") +
  geom_point(aes(x=can, y=f1.y), color="#16a085", size=2 ) +
  geom_point(aes(x=can, y=f1.x), color="#8e44ad", size=2 ) + xlab("") + ylab("F1") + 
  #geom_point(aes(x=can, y=f1), color="#e74c3c", size=2 ) + 
  scale_x_discrete(label=function(x) abbreviate(x, minlength=15)) + coord_flip() + theme_minimal()


# show the gains
tmp = maxPairs[f1 > f1.y,]
tmp$can <- factor(tmp$can, levels=rev(unique(tmp[order(-f1.y), ]$can)))
lali = ggplot(tmp) + geom_segment(aes(x=can, xend=can, y=f1.x, yend=f1.y), size=1.5, color="#dadfe3") +
  geom_point(aes(x=can, y=f1.y), color="#34495e", size=2 ) +
  geom_point(aes(x=can, y=f1.x), color="#dadfe3", size=2 ) + xlab("") + ylab("F1") + ylim(0, 1.5) +
  geom_segment(aes(x=can, xend=can, y=f1.y, yend=f1), size=1.5, color="#34495e") +
  geom_point(aes(x=can, y=f1, color=ttype), size=2 ) + scale_color_manual(values=c("#e67e22", "#16a085", "#2980b9")) +
  scale_x_discrete(label=function(x) abbreviate(x, minlength=15)) + coord_flip() + 
  geom_text(data=tmp, aes(label=triple, x=can, y=f1 + 0.05), hjust = 0, check_overlap = T) + theme_minimal() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.25), legend.position = "bottom")
lali
ggsave(paste0("scripts/clustering-method/res-figures/increased-lolli.pdf"),
       lali, width=6.5, height=4, units="in", device="pdf", useDingbats=FALSE)

# show all of them
tmp = maxPairs
tmp$can <- factor(tmp$can, levels=rev(unique(tmp[order(-f1.y), ]$can)))
lali = ggplot(tmp) + geom_segment(aes(x=can, xend=can, y=f1.x, yend=f1.y), size=1.5, color="#dadfe3") +
  geom_point(aes(x=can, y=f1.y), color="#34495e", size=2 ) +
  geom_point(aes(x=can, y=f1.x), color="#dadfe3", size=2 ) + xlab("") + ylab("F1") + ylim(0, 1.5) +
  geom_segment(aes(x=can, xend=can, y=f1.y, yend=f1), size=1.5, color="#34495e") +
  geom_point(aes(x=can, y=f1, color=ttype), size=2 ) + scale_color_manual(values=c("#e67e22", "#16a085", "#2980b9")) +
  scale_x_discrete(label=function(x) abbreviate(x, minlength=15)) + coord_flip() + 
  geom_text(data=tmp, aes(label=triple, x=can, y=f1 + 0.05), hjust = 0, check_overlap = T) + theme_minimal() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.25), legend.position = "bottom")
lali
ggsave(paste0("scripts/clustering-method/res-figures/increased-lolli-all.pdf"),
       lali, width=6.5, height=4, units="in", device="pdf", useDingbats=FALSE)


# updated lollipops to show best C to N and doubles without NNs
df1 = df[category == "single",]
df1 = df1[ttype == "C",]
df2 = df[category == "single",]
df2 = df2[ttype == "N",]

maxC = df1[df1[, .I[which.max(f1)], by=can]$V1]
maxN = df2[df2[, .I[which.max(f1)], by=can]$V1]
rm(df1, df2)
maxPairs = merge(maxC, maxN, by="can")

names(maxPairs)[names(maxPairs) =='pair.x'] <- 'clinical'
names(maxPairs)[names(maxPairs) =='pair.y'] <- 'novel'
maxPairs[,category.x := NULL]
maxPairs[,category.y := NULL]

maxPairs$can <- factor(maxPairs$can, levels=rev(unique(maxPairs[order(-f1.y), ]$can)))
gglolli = ggplot(maxPairs) + geom_segment(aes(x=can, xend=can, y=f1.y, yend=f1.x), size=1.5, color="#dadfe3") +
  geom_point(aes(x=can, y=f1.x, color=ttype.x), size=2 ) + scale_color_manual(values=c("#8e44ad", "#e67e22", "#16a085", "#2980b9")) +
  geom_text(aes(label=clinical, x=can, y=f1.x),hjust=1.2, vjust=0.5)  +
  geom_point(aes(x=can, y=f1.y), color="#16a085", size=2) +
  #geom_text(aes(label=novel, x=can, y=f1.y),hjust=-1, vjust=0.5) + xlab("") + ylab("F1") + 
  scale_x_discrete(label=function(x) abbreviate(x, minlength=15)) + coord_flip() +
  geom_text(data=maxPairs, aes(label=novel, x=can, y= f1.y+ .05), hjust = 0, check_overlap = T) + theme_minimal() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.25), legend.position = "bottom")
gglolli
ggsave(paste0("scripts/clustering-method/res-figures/paper/maxF1LollipopCtoN.pdf"), 
       gglolli, width=5, height=6, units="in", device="pdf", useDingbats=FALSE)


# simple just show best C pair to best double (excluding NNs)





# make plots for best potential doubles
doubles2 = doubles2[order(-score2),]
for (cn in cans)
{
  print(cn)
  apply (head(doubles2[can == cn & ctype == "HH"], 15), 1, function(row) { 
    genes = unlist(strsplit(as.character(row["pair"]), ":"))  
    mSP(cn, genes, printToFile = T)
  })
  apply (head(doubles2[can == cn & ctype != "HH"], 15), 1, function(row) { 
    genes = unlist(strsplit(as.character(row["pair"]), ":"))  
    mSP(cn, genes, printToFile = T)
  })
}

# make plots for top triples
tttmp = toptriples[ttype != "N:N:N",]
tttmp = tttmp[order(-f1),]
tttmp = tttmp[order(-score),]
for (cn in cans)
{
  print(cn)
  apply (head(tttmp[can == cn,], 5), 1, function(row) { 
    genes = unlist(strsplit(as.character(row["pair"]), ":"))
    if (sum(genes %in% names(rnaseqmat)) == 3)
    {
      make3DPlots(cn, genes, F, T)
    }
    else
    {
      print(genes)
    }
  })
}
rm(tttmp)

# other top triple plots
makeSplitViolins("Mesothelioma", c("CA9", "FAP", "KISS1R"))
makeSplitViolins.facet("Mesothelioma", c("CA9", "FAP", "KISS1R"))


makeSplitViolins("Mesothelioma", c("CA9", "FAP"))
makeSplitViolins("Mesothelioma", c("CA9", "KISS1R"))

target = "Mesothelioma"
genes =  c("CA9", "FAP", "KISS1R")


makeBeeswarms("Mesothelioma", c("CA9", "FAP", "KISS1R"))

plotDensity("Kidney Renal Clear Cell Carcinoma", "CD70")
plotDensity("Kidney Renal Clear Cell Carcinoma", "AXL")

plotDensity("Brain Lower Grade Glioma", "SYT11")
plotDensity("Brain Lower Grade Glioma", "EPCAM")
# maybe do beeswarm on sketch??
#sketchL = sketchLAll[can==can]$dataset
#sketch = rnaseqmat[dataset %in% sketchL,]

mSP("Mesothelioma", c("FAP", "CA9"), T)
mSP("Mesothelioma", c("FAP", "KISS1R"), T)

mSP.withCutoff("Colon Adenocarcinoma", c("CA9", "GRIN2D"), 7.5, T)
mSP.withCutoff("Liver Hepatocellular Carcinoma", c("ASGR1", "GABRD"), 9, T)
mSP.withCutoff("Skin Cutaneous Melanoma", c("MUC1", "MLANA"), 9, T)
mSP.withCutoff("Glioblastoma Multiforme", c("B4GALNT1", "AOC3"), 8.5, T)
mSP.withCutoff("Testicular Germ Cell Tumors", c("LEPR", "PTPRN"), 8, T)
mSP.withCutoff("Lung Squamous Cell Carcinoma", c("CA9", "ABCA8"), 10, T)
mSP.withCutoff("Pancreatic Adenocarcinoma", c("B4GALNT1", "MSLN"), 8, T)
mSP.withCutoff("Lung Squamous Cell Carcinoma", c("CA9", "B4GALNT1"), 7.5, T)

# doubles for 3D
mSP.nolabels("Glioblastoma Multiforme", c("B4GALNT1", "NTSR2"), T)
mSP.nolabels("Glioblastoma Multiforme", c("B4GALNT1", "GPR19"), T)

mSP.nolabels("Head and Neck Squamous Cell Carcinoma", c("CA9", "KREMEN2"), T)
mSP.nolabels("Head and Neck Squamous Cell Carcinoma", c("CA9", "PROM1"), F)

mSP.nolabels("Head and Neck Squamous Cell Carcinoma", c("KREMEN2", "ULBP2"), T)
mSP.nolabels("Head and Neck Squamous Cell Carcinoma", c("KREMEN2", "PROM1"), T)

mSP.nolabels("Mesothelioma", c("CA9", "CD70"), T)
mSP.nolabels("Mesothelioma", c("CA9", "FAP"), T)

mSP.nolabels("Skin Cutaneous Melanoma", c("MC1R", "EGFR"), T)
mSP.nolabels("Skin Cutaneous Melanoma", c("GPR143", "EGFR"), T)

mSP.nolabels("Pancreatic Adenocarcinoma", c("FAP", "CA9"), T)
mSP.nolabels("Pancreatic Adenocarcinoma", c("CELSR3", "CA9"), T)





make3DPlots("Mesothelioma", c("CA9", "FAP", "KISS1R"))

plotDensity("Lung Adenocarcinoma", "CA9")
plotDensity("Lung Adenocarcinoma", "TREM1")

#KISS1R, PRR7, LIFR, PTPRN, OR2B6
plotDensityAll("KISS1R")
plotDensityAll("PRR7")
plotDensityAll("LIFR")
plotDensityAll("PTPRN")
plotDensityAll("OR2B6")
plotDensityAll("GRIN2D")

plotDensity("Pancreatic Adenocarcinoma", "FAP")
plotDensity("Pancreatic Adenocarcinoma", "MSLN")

plotDensity("Kidney Renal Clear Cell Carcinoma", "CDH6")
plotDensity("Kidney Renal Clear Cell Carcinoma", "AXL")

# look at top connected CNs and print out their subnets
topCN = dbls2[ttype == "C:N",]

gs = topCN[g1lab == "N",]$gene1
gs = append(gs, topCN[g2lab == "N",]$gene2)
Ncounts = data.table(gene=gs)
Ncounts = Ncounts[,.N, by="gene"]
Ncounts = Ncounts[order(-N),]
rm(gs)

col = c("event" = "#2980b9", "actor" = "#e67e22")
Ncounts = head(Ncounts, 10)
totcan = numeric()
for(gene in Ncounts$gene)
{
  print(gene)
  tc = length(unique(topCN[gene1==gene | gene2==gene,.N, by="can"]$can)) # pairs for 23 cans
  totcan = append(totcan, tc)
  
  Cs = topCN[gene1==gene,]$gene2
  Cs = append(Cs, topCN[gene2==gene,]$gene1)
  Cs = data.table(gene=Cs)
  Cs = Cs[,.N, by="gene"]
  Cs = Cs[order(-N),]
  
  net = data.frame(target = Cs$N,
                   dummy = rep(0, length(Cs$N)),
                   row.names = Cs$gene)
  
  net = network(net,
                matrix.type = "bipartite",
                ignore.eval = FALSE,
                names.eval = "weights")
  ggn = ggnet2(net, size = 4, color = "mode", label = TRUE, label.size = 5, edge.label = "weights", edge.size = "weights", palette=col) + ggtitle(gene)
  
  ggsave(paste0("scripts/clustering-method/res-figures/paper/net-", gene,"-CN.pdf"),
         ggn, width=5.5, height=4.5, units="in", device="pdf", useDingbats=FALSE)
  
}
Ncounts[,tc := 1]
Ncounts$tc = totcan


# *******************************
# stats for the paper
# *******************************
tmp = df[ttype %in% c("C", "C:N", "C:C:N"),]
wilcox.test(f1 ~ ttype, data=tmp[ttype %in% c("C", "C:N"),])$p.value
wilcox.test(f1 ~ ttype, data=tmp[ttype %in% c("C:N", "C:C:N"),])$p.value
summary(tmp[ttype=="C",]$f1)
summary(tmp[ttype=="C:N",]$f1)
summary(tmp[ttype=="C:C:N",]$f1)


tmp = df[ttype %in% c("C", "N", "C:C", "C:N", "N:N"),]
wilcox.test(f1 ~ ttype, data=tmp[ttype %in% c("N", "C:N"),])$p.value
wilcox.test(f1 ~ ttype, data=tmp[ttype %in% c("C", "C:C"),])$p.value
summary(tmp[ttype=="N:N",]$f1)
summary(tmp[ttype=="C:N",]$f1)
summary(tmp[ttype=="C:C",]$f1)
summary(tmp[ttype=="C",]$f1)
summary(tmp[ttype=="N",]$f1)

summary(doubles2[can=="Cholangiocarcinoma"]$f1)
summary(toptriples[can=="Cholangiocarcinoma",]$f1)


#topsingles = singles[, head(.SD, 10), by=c("can", "glab")] # not accurate b/c want all clinicals
topsingles = tmp[category == "single",]
topdoubles = doubles2[, head(.SD, 10), by=c("can", "ttype")]
toptriples2 = toptriples[, head(.SD, 10), by=c("can", "ttype")]

summary(topsingles$prec)
summary(topdoubles$prec)
summary(toptriples2$prec)

summary(topsingles$rec)
summary(topdoubles$rec)
summary(toptriples2$rec)

mean(maxPairs$f1.y - maxPairs$f1.x)
mean(maxPairs$prec.y - maxPairs$prec.x)
mean(maxPairs$rec.y - maxPairs$rec.x)

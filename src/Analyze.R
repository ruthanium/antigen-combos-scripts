# **************************************
# Author: Ruth Dannenfelser
#
# Analyze the top combinations, make
# result figures.
# Assumes the functions in Functions.R
# are loaded and the *.total data.tables 
# have been built from Evaluate.R
# **************************************
library(ggplot2)
library(ggrepel)
library(ggbeeswarm)
library(plotly)
library(network)
library(GGally)

# load in the combined clustering + evaluation scores
singles.total = readRDS("results/singles.total.RData")
doubles.total = readRDS("results/doubles.total.RData")
triples.total = readRDS("results/triples.total.RData")

# scale the scores between 0 and 1 for all combinations
st = cleanAndScale(singles.total)
dt = cleanAndScale(doubles.total)
tt = cleanAndScale(triples.total)

# get the top k per gate type per cancer allowing
# a gene to occur in a combination no more than
# two times (this prevents the top lists from being
# dominated by a pair or triplet with one fixed gene)
k = 10 # number of top to take

# get top singles
topSingles = st[, head(.SD, k), by=c("can","lab")]
maxClin = st[st[lab=="C", .I[which.max(f1)], by=can]$V1] # best performing clinical per tumor type

# get top doubles
dt[, lab := ifelse(lab == "N:C", "C:N", lab)] # standardize the labels so tree types are consistent
treeTypes = c("N:N", "C:N", "C:C")
topDoublesNames = rbindlist(lapply(treeTypes, function(tt) {
  rbindlist(mclapply(cans, function(cn,tt)
  {
    print(cn)
    one = character()
    two = character()
    count = 0
    
    res = dt[can==cn,]
    res = res[lab == tt]
    
    tmp.combo = character()
    tmp.can = character()
    
    i = 1;
    
    # genes can only be present a maximum of 2 times
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
        
        tmp.combo = append(tmp.combo, as.character(row$combo))
        tmp.can = append(tmp.can, as.character(row$can))
      }
    }
    data.table(combo=tmp.combo, can=tmp.can)
    
  }, tt, mc.cores=length(unique(dt$can))))
}))
topDoubles = merge(topDoublesNames, dt, by=c("combo", "can"))
rm(topDoublesNames)

# get top triples
# standardize the labels so tree types are consistent
tt[, lab := ifelse(lab == "C:N:C", "C:C:N", lab)] 
tt[, lab := ifelse(lab == "N:C:C", "C:C:N", lab)] 
tt[, lab := ifelse(lab == "N:N:C", "C:N:N", lab)] 
tt[, lab := ifelse(lab == "N:C:N", "C:N:N", lab)] 

treeTypes = c("N:N:N", "C:C:C", "C:C:N", "C:N:N")
topTriplesNames = rbindlist(lapply(treeTypes, function(ttree) {
  rbindlist(mclapply(cans, function(cn,ttree)
  {
    print(cn)
    one = character()
    two = character()
    count = 0
    
    res = tt[can==cn,]
    res = res[lab == ttree]
    
    tmp.combo = character()
    tmp.can = character()
    
    i = 1;
    
    # genes can only be present a maximum of 2 times
    while ((count < k)  & (i < nrow(res)))
    {
      row = res[i,]
      i = i + 1;
      genes = c(row$gene1, row$gene2, row$gene3)
      
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
        
        tmp.combo = append(tmp.combo, as.character(row$combo))
        tmp.can = append(tmp.can, as.character(row$can))
      }
    }
    data.table(combo=tmp.combo, can=tmp.can)
    
  }, ttree, mc.cores=length(unique(tt$can))))
}))
topTriples = merge(topTriplesNames, tt, by=c("combo", "can"))
rm(topTriplesNames)



# -------------------
#       plots
# -------------------
topSingles[, category := "single"]
topDoubles[, category := "double"]
topTriples[, category := "triple"]
df = rbind(topSingles, topDoubles[,.SD, .SDcols = !c('gene1', 'gene2', 'g1lab', 'g2lab')], 
           topTriples[,.SD, .SDcols = !c('gene1', 'gene2', 'gene3', 'g1lab', 'g2lab', 'g3lab')])

# F1 changes from single to double by antigen type 
tmp = df[lab %in% c("C", "N", "C:C", "C:N", "N:N"),]
tmp$lab = factor(tmp$lab, levels=c("C", "N", "C:C", "C:N", "N:N"))
ggplot(tmp, aes(x=lab, y=f1)) + geom_boxplot(outlier.shape = NA, notch=T) + 
  theme_bw() + ylim(0,1) + xlab("") + ylab("F1 score") + 
  theme(legend.position="bottom", axis.title = element_text(size=12), 
        axis.text= element_text(size=12), panel.grid.minor=element_blank(),
        plot.background=element_blank(), panel.background=element_blank(), panel.grid.major=element_blank())

wilcox.test(f1 ~ lab, data=tmp[lab %in% c("C", "N"),])$p.value
wilcox.test(f1 ~ lab, data=tmp[lab %in% c("N", "C:C"),])$p.value
wilcox.test(f1 ~ lab, data=tmp[lab %in% c("C:C", "C:N"),])$p.value
wilcox.test(f1 ~ lab, data=tmp[lab %in% c("C:N", "N:N"),])$p.value

rm(tmp)

# best C to best pair (score)
df1 = df[category == "single" & lab =="C",]
df2 = df[category == "double",]
df3 = df[category == "triple",]

maxS = df1[df1[, .I[which.max(f1)], by=can]$V1]
maxD = df2[df2[, .I[which.max(f1)], by=can]$V1]
maxT = df3[df3[, .I[which.max(f1)], by=can]$V1]
rm(df1, df2, df3)
maxPairs = merge(maxS, maxD, by="can")
maxPairs = merge(maxPairs, maxT, by="can")

names(maxPairs)[names(maxPairs) =='combo.x'] <- 'single'
names(maxPairs)[names(maxPairs) =='combo.y'] <- 'double'
names(maxPairs)[names(maxPairs) =='combo'] <- 'triple'
maxPairs[,category.x := NULL]
maxPairs[,category.y := NULL]
maxPairs[,category := NULL]

# best C pair to best double
maxPairs$can <- factor(maxPairs$can, levels=rev(unique(maxPairs[order(-f1.y), ]$can)))
ggplot(maxPairs) + geom_segment(aes(x=can, xend=can, y=score.x, yend=f1.y), size=1.5, color="#dadfe3") +
  geom_point(aes(x=can, y=f1.y), color="#16a085", size=2 ) +
  geom_point(aes(x=can, y=f1.x), color="#8e44ad", size=2 ) + xlab("") + ylab("F1") + 
  scale_x_discrete(label=function(x) abbreviate(x, minlength=15)) + coord_flip() + theme_minimal()

# 2D plot
makeScatterPlot("Glioblastoma Multiforme", c("OR4F5", "SLC2A5"), rnaseqmat, T)

# 3D plot 
make3DPlots("Glioblastoma Multiforme", c("SLC2A5", "CLCNKA", "EPCAM"), rnaseqmat)

# tradeoff of increasing number of antigens in the combination
df2 = melt(df[,c("combo", "can", "category", "prec", "rec")], id.vars=c("combo", "can", "category"))
df2$category = factor(df2$category, levels=c("single", "double", "triple"))
gg = ggplot(df2, aes(x=category, y=value, color=variable)) + geom_quasirandom(position ="dodge", alpha=0.8) + theme_bw() + facet_wrap(~variable) +
  scale_color_manual(values=c("#999999", "#999999")) + ylab("") + xlab("") + 
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 0.5, color="black")
gg
rm(df2)

# look at top connected CNs and print out their subnets
topCN = topDoubles[lab == "C:N",]

gs = topCN[g1lab == "N",]$gene1
gs = append(gs, topCN[g2lab == "N",]$gene2)
Ncounts = data.table(gene=gs) # counts for each N in the C:N pair
Ncounts = Ncounts[,.N, by="gene"]
Ncounts = Ncounts[order(-N),]
rm(gs)

col = c("event" = "#2980b9", "actor" = "#e67e22")
Ncounts = head(Ncounts, 3) # look at top 3
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
  ggnet2(net, size = 4, color = "mode", label = TRUE, label.size = 5, edge.label = "weights", edge.size = "weights", palette=col) + ggtitle(gene)
  
}
Ncounts[,tc := 1]
Ncounts$tc = totcan



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




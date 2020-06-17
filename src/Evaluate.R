# *************************************
# Author: Ruth Dannenfelser
#
# Using decision trees calculate 
# evaluation metrics for single,
# double, and triple antigen gates.
#
# warning: will take a long time to 
# run to calculate all combinations
#**************************************
library(data.table)
library(rpart)

# load prereqs
# start with data.table of dataset x gene expression values
rnaseqmat = fread("data/test-normalized-matrix.txt")
rnaseqmat = as.data.table(rnaseqmat) # make sure it is data.table and not data.frame
rnaseqmat[,"batch" := NULL] # remove batch - b/c redundant info with type

cans = unique(rnaseqmat[rnaseqmat$type=="cancer",]$tissue.cancer)

# if not already loaded
singles = readRDS("results/singles.RData")
doubles = readRDS("results/doubles.RData")
triples = readRDS("results/triples.RData")

# calculate metrics for singles
singles = singles[ttype == "H",] # filter out Ls (as they are not targetable)
singles.perf = rbindlist(mclapply(cans, function(cann)
{
  f1=numeric()
  prec = numeric()
  rec = numeric()
  tcan = character()
  gene = character()
  for (g in unique(singles[can==cann,]$combo))
  {
    acl = assignClusterLabels(cann, c(g), rnaseqmat)
    f1 = append(f1, acl$F1)
    prec = append(prec, acl$prec)
    rec = append(rec, acl$recall)
    tcan = append(tcan, cann)
    gene = append(gene, g)
    
  }
  data.table(combo=gene, can=tcan, f1, prec, rec)
}, mc.cores=length(cans)))

singles.total = merge(singles, singles.perf, by=c("combo", "can"))
rm(singles.perf)

# calculate metrics for doubles
doubles = doubles[ttype != "LL",] # filter out Ls (as they are not targetable)

doubles.perf = rbindlist(mclapply(cans, function(cann)
{
  f1=numeric()
  prec = numeric()
  rec = numeric()
  tcan = character()
  gene = character()
  for (g in unique(doubles[can==cann,]$combo))
  {
    genes = unlist(strsplit(g, ":"))
    acl = assignClusterLabels(cann, genes, rnaseqmat)
    f1 = append(f1, acl$F1)
    prec = append(prec, acl$prec)
    rec = append(rec, acl$recall)
    tcan = append(tcan, cann)
    gene = append(gene, g)
    
  }
  data.table(combo=gene, can=tcan, f1, prec, rec)
}, mc.cores=length(cans)))

doubles.total = merge(doubles, doubles.perf, by=c("combo", "can"))
rm(doubles.perf)

# calculate metrics for triples
triples = triples[ttype != "LLL",] # filter out Ls (as they are not targetable)

triples.perf = rbindlist(mclapply(cans, function(cann)
{
  f1=numeric()
  prec = numeric()
  rec = numeric()
  tcan = character()
  gene = character()
  for (g in unique(triples[can==cann,]$combo))
  {
    genes = unlist(strsplit(g, ":"))
    acl = assignClusterLabels(cann, genes, rnaseqmat)
    f1 = append(f1, acl$F1)
    prec = append(prec, acl$prec)
    rec = append(rec, acl$recall)
    tcan = append(tcan, cann)
    gene = append(gene, g)
    
  }
  data.table(combo=gene, can=tcan, f1, prec, rec)
}, mc.cores=length(cans)))

triples.total = merge(triples, triples.perf, by=c("combo", "can"))
rm(triples.perf)

saveRDS(singles.total, "results/singles.total.RData")
saveRDS(doubles.total, "results/doubles.total.RData")
saveRDS(triples.total, "results/triples.total.RData")


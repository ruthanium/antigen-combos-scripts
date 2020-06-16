# *************************************
# Author: Ruth Dannenfelser
#
# Using decision trees calculate 
# evaluation metrics for single,
# double, and triple antigen gates.
#**************************************
library(data.table)
library(rpart)


# start with data.table of dataset x gene expression values
rnaseqmat = fread("data/test-normalized-matrix.txt")
rnaseqmat = as.data.table(rnaseqmat) # make sure it is data.table and not data.frame
rnaseqmat[,"batch" := NULL] # remove batch - b/c redundant info with type

cans = unique(rnaseqmat[rnaseqmat$type=="cancer",]$tissue.cancer)

# calculate metrics for singles
singles = singles[ttype == "H",] # filter out Ls (as they are not targetable)
single.perf = rbindlist(mclapply(cans, function(cann)
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



# *************************************
# Author: Ruth Dannenfelser
#
# Calculate clustering scores 
# run using parallelization
#
#**************************************
library(parallel)
library(data.table)

library(clusterSim)
library(ggplot2)
library(cluster)
library(clusterCrit)

# start with data.table of dataset x gene expression values
rnaseqmat = fread("data/test-normalized-matrix.txt")
rnaseqmat = as.data.table(rnaseqmat) # make sure it is data.table and not data.frame
rnaseqmat[,"batch" := NULL] # remove batch - b/c redundant info with type

cans = unique(rnaseqmat[rnaseqmat$type=="cancer",]$tissue.cancer)

clinical = fread("data/clin-antigens.txt", header=F)
clinical = clinical$V1
tm = fread("data/transmembrane-genes.txt", header=F) # already contains the clinicals
tm = tm$V1

# calculate clustering scores for single antigens
singles = (mclapply(cans, calcClustMetricsForSingles, mat=rnaseqmat, mc.cores=2))
setattr(singles, 'names', cans)
singles = rbindlist(singles, use.names=T, fill=T, idcol=T)
names(singles) = c("can", "combo", "db", "dist.man", "ttype")


# label clinicals vs novels and remove Ls
saveRDS(singles, "results/singles.RData")

# load in the sketches
sketchLAll = rbindlist(lapply(cans,
                              function(cn) data.table(can=cn,
                                                      dataset=readLines(paste0("results/sketches/", gsub(" ", "-", tolower(cn)), "-sketch.txt")))))









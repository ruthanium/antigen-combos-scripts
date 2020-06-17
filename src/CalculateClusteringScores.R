# *************************************
# Author: Ruth Dannenfelser
#
# Calculate clustering scores 
# run using parallelization
#
#**************************************
library(parallel)
library(data.table)

tr = F # if true, test run (limit number of doubles and triples calculated) 

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
singles = (mclapply(cans, calcClustMetricsForSingles, mat=rnaseqmat, mc.cores=length(cans)))
setattr(singles, 'names', cans)
singles = rbindlist(singles, use.names=T, fill=T, idcol=T)
names(singles) = c("can", "combo", "db", "dist.man", "ttype")

# label clinicals vs novels
singles[, lab := "N"]
singles[combo %in% clinical,]$lab = "C"

saveRDS(singles, "results/singles.RData")

# calculate clustering scores for doubles
doubles = (mclapply(cans, calcClustMetrics, mat=rnaseqmat, testRun=tr, mc.cores=length(cans)))
setattr(doubles, 'names', cans)
doubles = rbindlist(doubles, use.names=T, fill=T, idcol=T)
names(doubles) = c("can", "combo", "db", "dist.man", "ttype")

# label clinicals vs novels
doubles[, c("gene1", "gene2") := tstrsplit(combo, ":", fixed=TRUE)]
doubles[, g1lab := ifelse(gene1 %in% clinical, "C", "N")]
doubles[, g2lab := ifelse(gene2 %in% clinical, "C", "N")]
doubles[, lab := paste0(g1lab, ":", g2lab)]

saveRDS(doubles, "results/doubles.RData")

# calculate clustering scores for triples
singfilt = singles[db <= 5 & dist.man > 2,] # reduction to top performing singles
triples = (mclapply(cans, calcClustMetricsTriples, mat=rnaseqmat, singles=singfilt, testRun=tr, mc.cores=length(cans)))
setattr(triples, 'names', cans)
triples = rbindlist(triples, use.names=T, fill=T, idcol=T)
names(triples) = c("can", "combo", "db", "dist.man", "ttype")

# label clinicals vs novels
triples[, c("gene1", "gene2", "gene3") := tstrsplit(combo, ":", fixed=TRUE)]
triples[, g1lab := ifelse(gene1 %in% clinical, "C", "N")]
triples[, g2lab := ifelse(gene2 %in% clinical, "C", "N")]
triples[, g3lab := ifelse(gene3 %in% clinical, "C", "N")]
triples[, lab := paste0(g1lab, ":", g2lab, ":", g3lab)]

saveRDS(triples, "results/triples.RData")








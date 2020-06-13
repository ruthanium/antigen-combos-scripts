# *************************************
# Author: Ruth Dannenfelser
#
# Calculate clustering scores and 
# running using parallelization
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

# make sure it is data.table and not data.frame
rnaseqmat = as.data.table(rnaseqmat)






rnaseqmat[rnaseqmat$type=="tumor",]$type = "cancer"
params = unique(rnaseqmat[rnaseqmat$type=="cancer",]$tissue.cancer)
rnaseqmat = within(rnaseqmat, rm(batch)) # remove batch


sketchLAll = rbindlist(lapply(params,
                              function(cn) data.table(can=cn,
                                                      dataset=readLines(paste0("scripts/clustering-method/sketches/", gsub(" ", "-", tolower(cn)), "-sketch.txt")))))

clinical = fread("data/transmembrane-markers/clin-paper-dec-19.txt", header=F)
clinical = clinical$V1

tm = fread("data/transmembrane-markers/filteredTransmembraneList-updated-12-19.txt", header=F)
tm = tm$V1

singles = readRDS("scripts/clustering-method/rdata/singles-complete.RData")
singles = singles[db <= 5 & dist.man > 2,] # reduction to top performing singles

#subset = names(rnaseqmatClean)[1:500]
#subset = append(subset, c("type", "tissue.cancer"))
#sub = rnaseqmatClean[,subset, with=F]
#mclapply(params, calcClustMetricsFaster, am=rnaseqmat, fnpost="-rnaseq-combat-sketch-0.2-triples", mc.cores=20)

# test to make sure things are working correctly 
tmpres = calcClustMetricsSinglesReduction("Lung Adenocarcinoma", rnaseqmatClean, "-triple-test", T)
rm(tmpres)

mclapply(params, calcClustMetricsSinglesReduction, am=rnaseqmatClean, fnpost="-rnaseq-combat-sketch-0.2-triples-S-reduction", mc.cores=33)


arraymat = readRDS("scripts/clustering-method/rdata/arraymat-clean-transmembrane.RData")
params = unique(arraymat[arraymat$type=="cancer",]$tissue.cancer)

bigres = mclapply(params, calcClustMetricsFaster, am=arraymat, fnpost="-array-vs-normals-res-allTM-fixed", mc.cores=(length(params)*1.5))

tm = fread("data/transmembrane-markers/filteredTransmembraneList-updated-12-19.txt", header=F)
tm = tm$V1
tm = tm[tm %in% names(rnaseqmat)]


#rnaseqmat = readRDS("scripts/clustering-method/rnaseq-tm-clin.RData")
#params = unique(rnaseqmat[rnaseqmat$type=="cancer",]$tissue.cancer)

#bigres = mclapply(params, calcClustMetricsFaster, am=rnaseqmat, fnpost="-rnaseq-res-allTM", mc.cores=(length(params)*1.5))


# updated to run on the combat corrected RNAseq in TPMs
rnaseqmat = readRDS("scripts/clustering-method/rdata/combat-rnaseq-tm-clin-zeroadj.RData")

rnaseqmat = fread("scripts/clustering-method/rnaseq-mat-combat-zeroadj.txt", sep="\t")

rnaseqmat[rnaseqmat$type=="tumor",]$type = "cancer"
params = unique(rnaseqmat[rnaseqmat$type=="cancer",]$tissue.cancer)
rnaseqmat = within(rnaseqmat, rm(batch)) # remove batch

subset = names(rnaseqmat)[1:100]
subset = append(subset, c("type", "tissue.cancer"))

sub = rnaseqmat[,subset, with=F]

mclapply(params, calcClustMetricsFaster, am=rnaseqmat, fnpost="-rnaseq-combat-sketch-0.2", mc.cores=(length(params)*1.5))

mclapply(params, calcClustMetricsFasterClin, am=rnaseqmat, fnpost="-rnaseq-combat-sketch-0.2-clin", mc.cores=(length(params)*1.5))
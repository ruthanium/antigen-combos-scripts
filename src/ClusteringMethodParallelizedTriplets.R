# *************************************
# Author: Ruth Dannenfelser
# Date: September 30, 2019
#
# More systematic clustering metric
# script - extended to triplets
#
#**************************************
library(parallel)
library(data.table)
library(clusterSim)
library(ggplot2)
library(cluster)
#library(clusterCrit)

setwd("/Genomics/ogtr04/rd6/BioGPS/")

# updated to run on the combat corrected RNAseq in TPMs
rnaseqmat = readRDS("scripts/clustering-method/rdata/combat-rnaseq-tm-clin-zeroadj.RData")

#rnaseqmat = fread("scripts/clustering-method/rnaseq-mat-combat-zeroadj.txt", sep="\t")

rnaseqmat[rnaseqmat$type=="tumor",]$type = "cancer"
params = unique(rnaseqmat[rnaseqmat$type=="cancer",]$tissue.cancer)
rnaseqmat = within(rnaseqmat, rm(batch)) # remove batch

rnaseqmatClean = copy(rnaseqmat)
rnaseqmatClean[rnaseqmatClean < 0] = 0 # set values < 0 in the rnaseqmat to 0

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



calcClustMetricsSinglesReduction <- function(target, am, fnpost, testRun=F)
{
  dt = am[am$type == "tissue" | am$tissue.cancer == target,]
  sketch = sketchLAll[can==target]$dataset
  dt = dt[dataset %in% sketch,]
  
  sn = singles[can==target,]
  combos = as.data.frame(t(combn(as.character(sn$gene), 3)))
  if(testRun)
  {
    combos = combos[1:50,] 
  }
  ares = apply(combos, 1, calcscore, dtt=dt, target=target)
  res = do.call(rbind.data.frame, ares)
  res = as.data.table(res)
  names(res) = c("pair", "db", "dist.man")
  
  saveRDS(res, paste0("scripts/clustering-method/rdata/",gsub(" ", "_", target), fnpost, ".RData"))
  return (res);
}


calcClustMetricsFasterClin <- function(target, am, fnpost)
{
  dt = am[am$type == "tissue" | am$tissue.cancer == target,]
  sketch = sketchLAll[can==target]$dataset
  dt = dt[dataset %in% sketch,]
  
  tmpc = names(dt)[2:ncol(dt)-2]
  tmpc = tmpc[-1]
  tmpc = tmpc[tmpc %in% clinical]
  combos = as.data.frame(t(combn(tmpc, 3)))
  rm(tmpc)
  
  ares = apply(combos, 1, calcscore, dtt=dt, target=target)
  res = do.call(rbind.data.frame, ares)
  res = as.data.table(res)
  names(res) = c("pair", "db", "dist.man")
  
  saveRDS(res, paste0("scripts/clustering-method/rdata/",gsub(" ", "_", target), fnpost, ".RData"))
  return (res);
}

calcClustMetricsFaster <- function(target, am, fnpost)
{
  dt = am[am$type == "tissue" | am$tissue.cancer == target,]
  sketch = sketchLAll[can==target]$dataset
  dt = dt[dataset %in% sketch,]
  
  tmpc = names(dt)[2:ncol(dt)-2]
  tmpc = tmpc[-1]
  combos = as.data.frame(t(combn(tmpc, 3)))
  rm(tmpc)
  #tmpc = combos[1:50,] # need to change this
  ares = apply(combos, 1, calcscore, dtt=dt, target=target)
  res = do.call(rbind.data.frame, ares)
  res = as.data.table(res)
  names(res) = c("pair", "db", "dist.man")
  
  saveRDS(res, paste0("scripts/clustering-method/rdata/",gsub(" ", "_", target), fnpost, ".RData"))
  return (res);
}


calcscore <- function(row, dtt, target) { # calc scores
  dtt = as.data.table(dtt)
  genes = c(as.character(row["V1"]), as.character(row["V2"]), as.character(row["V3"]))
  tmp = dtt[,c(genes[1], genes[2], genes[3], "tissue.cancer"), with=F]
  tmp$target = ifelse(tmp$tissue.cancer == target, "target", "other")
  clust <- ifelse(tmp$target == "target", 1, 2)
  m = as.matrix(tmp[,1:3, with=FALSE])
  mode(m) = "numeric"
  db = index.DB.fixed(m, as.integer(clust), centrotypes = "centroids", p=1)$DB
  cm = colMeans(tmp[target == "other",1:3, with=F])
  cm = rbind(cm, colMeans(tmp[target == "target",1:3, with=F]))
  row.names(cm) = NULL
  cm = as.matrix(cm)
  cent = dist(cm, method="manhattan")[1]
  return (list(pair=paste(genes[1], genes[2], genes[3], sep=":"), db, cent))
}


## override index.DB - b/c package is broken
index.DB.fixed <- function(x,cl,d=NULL,centrotypes="centroids",p=2,q=2){
  if(sum(c("centroids","medoids")==centrotypes)==0)
    stop("Wrong centrotypes argument")
  if("medoids"==centrotypes && is.null(d))
    stop("For argument centrotypes = 'medoids' d cannot be null")
  if(!is.null(d)){
    if(!is.matrix(d)){
      d<-as.matrix(d)
    }
    row.names(d)<-row.names(x)
  }
  if(is.null(dim(x))){
    dim(x)<-c(length(x),1)
  }
  x<-as.matrix(x)
  n <- length(cl)
  k <- max(cl)
  #print(n)
  dAm<-d
  centers<-matrix(nrow=k,ncol=ncol(x))
  if (centrotypes=="centroids"){
    for(i in 1:k)
    {
      for(j in 1:ncol(x))
      {
        centers[i,j]<-mean(x[cl==i,j])
      }
    }
  }
  else if (centrotypes=="medoids"){
    #print("start")
    #print(dAm)
    for (i in 1:k){
      clAi<-dAm[cl==i,cl==i]
      if (is.null(clAi)){
        centers[i,]<-NULL
      }
      else{
        #print("przed centers")
        #print(x[cl==i,])
        #print(clAi)
        centers[i,]<-.medoid(x[cl==i,],dAm[cl==i,cl==i])
        #print("po centers")
        #print(centers[i])
      }
    }   
    #print("stop")
  }
  else{
    stop("wrong centrotypes argument")
  }
  S<-rep(0,k)
  for(i in 1:k){                             # For every cluster
    ind <- (cl==i)
    if (sum(ind)>1){
      centerI<-centers[i,]
      centerI<-rep(centerI,sum(ind))
      centerI<-matrix(centerI,nrow=sum(ind),ncol=ncol(x),byrow=TRUE)
      #S[i] <- mean(sqrt(apply((x[ind,] - centerI)^2,1,sum))^q)^(1/q)
      S[i] <- mean(apply((x[ind,] - centerI)^2,1,sum)^q)^(1/q)
      
    }
    else
      S[i] <- 0                         
  }
  M<-as.matrix(dist(centers, method="minkowski", p=p))
  R <- array(Inf,c(k,k))
  r = rep(0,k)
  for (i in 1:k){
    for (j in 1:k){
      R[i,j] = (S[i] + S[j])/M[i,j]
    }
    r[i] = max(R[i,][is.finite(R[i,])])
  } 
  DB = mean(r[is.finite(r)])        
  resul<-list(DB=DB,r=r,R=R,d=M,S=S,centers=centers)
  resul
}


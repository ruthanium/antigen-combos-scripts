# *************************************
# Author: Ruth Dannenfelser
# Date: Sept 17, 2019
#
# Rerun clustering metric
# with collapsed zeros for top pairs
#
#**************************************
library(parallel)
library(data.table)
library(clusterSim)
library(ggplot2)
library(cluster)
library(clusterCrit)

setwd("/Genomics/ogtr04/rd6/BioGPS/")

bigres = readRDS("scripts/clustering-method/rdata/all-rnaseq-combat-sketch-0.2.RData")
rnaseqmat = readRDS("scripts/clustering-method/rdata/combat-rnaseq-tm-clin.RData")
cans = unique(rnaseqmat[rnaseqmat$type == "tumor",]$tissue.cancer)

rnaseqmatClean = copy(rnaseqmat)
rnaseqmatClean[rnaseqmatClean$type=="tumor",]$type = "cancer"
params = unique(rnaseqmatClean[rnaseqmatClean$type=="cancer",]$tissue.cancer)
rnaseqmatClean = within(rnaseqmatClean, rm(batch)) # remove batch
rnaseqmatClean[rnaseqmatClean < 0] = 0 # set values < 0 in the rnaseqmat to 0


sketchLAll = rbindlist(lapply(cans,
                              function(can) data.table(can=can,
                                                       dataset=readLines(paste0("scripts/clustering-method/sketches/", gsub(" ", "-", tolower(can)), "-sketch.txt")))))

# get the pairs for each cancer for the already prioriized pairs
cutoff = 0.85
bigres[score >= cutoff,.N,by=.(type)]
pairs = bigres[score >= cutoff,c("pair", "type")]
names(pairs) = c("pair", "can")

trashh = mclapply(params, calcClustMetricsSubset, am=rnaseqmatClean, fnpost="-rnaseq-combat-sketch-0.2-corrected", pars=pairs, mc.cores=33)

# combine with new results in bigres
updated.bigres = rbindlist(lapply(cans, updateRes))
saveRDS(updated.bigres, "scripts/clustering-method/rdata/all-rnaseq-combat-sketch-0.2-top-fixed.RData")


updateRes <- function(can)
{
  br = bigres[type==can]
  
  # grab the new pairs
  updated = readRDS(paste0("scripts/clustering-method/rdata/", gsub(" ", "_", can), "-rnaseq-combat-sketch-0.2-corrected.RData"))
  # replace the db and dist.man columns
  br[pair %in% updated$pair,]$db = as.numeric(updated$db)
  br[pair %in% updated$pair,]$dist.man = as.numeric(updated$dist.man)
  
  br[,dbscaled := 1-((log(db) - min(log(db)))/(max(log(db)) - min(log(db))))] # make score between 0 and 1
  br[,distscaled := (log(dist.man) - min(log(dist.man)))/(max(log(dist.man)) - min(log(dist.man)))]
  br[,score := apply(br[,c("dbscaled", "distscaled"), with=F], 1, min)]
  br = br[order(-score),]
  
  # resort the results and return the updated res
  return(br)
}

calcClustMetricsSubset <- function(target, am, fnpost, pars)
{  
  dt = am[am$type == "tissue" | am$tissue.cancer == target,]
  sketch = sketchLAll[can==target]$dataset
  dt = dt[dataset %in% sketch,]

  ps = pairs[can==target]
  ps[, c("V1", "V2") := tstrsplit(pair, ":", fixed=TRUE)]
  
  ares = apply(ps[,c("V1", "V2")], 1, calcscore, dtt=dt, target=target)
  res = do.call(rbind.data.frame, ares)
  res = as.data.table(res)
  names(res) = c("pair", "db", "dist.man")
  
  saveRDS(res, paste0("scripts/clustering-method/rdata/",gsub(" ", "_", target), fnpost, ".RData"))
  return (res);
}



calcscore <- function(row, dtt, target) { # calc scores
  dtt = as.data.table(dtt)
  genes = c(as.character(row["V1"]), as.character(row["V2"]))
  tmp = dtt[,c(genes[1], genes[2], "tissue.cancer"), with=F]
  tmp$target = ifelse(tmp$tissue.cancer == target, "target", "other")
  clust <- ifelse(tmp$target == "target", 1, 2)
  m = as.matrix(tmp[,1:2, with=FALSE])
  mode(m) = "numeric"
  #db = index.DB.fixed(m, as.integer(clust), centrotypes = "centroids", p=2)$DB
  db = index.DB.fixed(m, as.integer(clust), centrotypes = "centroids", p=1)$DB
  cm = colMeans(tmp[target == "other",1:2, with=F])
  cm = rbind(cm, colMeans(tmp[target == "target",1:2, with=F]))
  row.names(cm) = NULL
  cm = as.matrix(cm)
  cent = dist(cm, method="manhattan")[1]
  return (list(pair=paste(genes[1], genes[2], sep=":"), db, cent))
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


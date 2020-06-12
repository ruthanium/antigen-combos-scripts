# *************************************
# Author: Ruth Dannenfelser
# Date: March 11, 2018
#
# More systematic clustering metric
# script - trying to improve the speed
#
#**************************************
library(parallel)
library(data.table)
library(clusterSim)
library(ggplot2)
library(cluster)
library(clusterCrit)

setwd("/Genomics/ogtr04/rd6/BioGPS/")

## REGENERATE ARRAYMAT USING DATA IN ARRAYDF
#arraydf = fread("data/GEO/db-input-formatted-filtered-array-clean.txt", sep="\t")
#annots = arraydf[,c("type", "dataset", "tissue.cancer")]
#annots = unique.data.frame(annots)
#
## dataset, genes, tissue.cancer, type
#arraymat = dcast.data.table(arraydf, dataset~gene, value.var="UPC")
#arraymat = merge(arraymat, annots, by="dataset")
#
#saveRDS(arraymat,"scripts/clustering-method/rdata/arraymat-clean-transmembrane.RData")

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



# rnaseqmatTG = rnaseqmat[rnaseqmat$type == "tissue",]
# # remove the weird tissue annotations
# rnaseqmatTG = rnaseqmatTG[rnaseqmatTG$tissue.cancer != "V7460",]
# rnaseqmatTG = rnaseqmatTG[rnaseqmatTG$tissue.cancer != "V7057",]
# params = unique(rnaseqmatTG$tissue.cancer)
# 
# bigres = mclapply(params, calcClustMetricsFasterTG, am=rnaseqmatTG, fnpost="-tissueGPS-res-allTM", mc.cores=(length(params)*1.5))



start.time <- Sys.time()
calcClustMetricsFaster("Skin Cutaneous Melanoma", sub, "-toss")
end.time <- Sys.time()
end.time - start.time

calcClustMetricsFasterSingles("Skin Cutaneous Melanoma", sub, "-singles", tm)

singTmp = (mclapply(params, calcClustMetricsFasterSingles, am=rnaseqmatClean, fnpost="-rnaseq-combat-sketch-0.2-singles", tm=tm, mc.cores=33))
setattr(singTmp, 'names', params)
singTmp = rbindlist(singTmp, use.names=T, fill=T, idcol=T)
names(singTmp) = c("can", "pair", "db", "dist.man", "ttype")
saveRDS(singTmp, "scripts/clustering-method/rdata/singles.RData")

calcClustMetricsFaster <- function(target, am, fnpost)
{
  dt = am[am$type == "tissue" | am$tissue.cancer == target,]
  fn =  file(paste0("scripts/clustering-method/sketches/", gsub(" ", "-", tolower(target)) ,"-sketch.txt"))
  sketch = readLines(fn)
  print(length(sketch))
  close(fn)
  
  dt = dt[dataset %in% sketch,]
  tmpc = names(dt)[2:ncol(dt)-2]
  tmpc = tmpc[-1]
  combos = as.data.frame(t(combn(tmpc, 2)))
  rm(tmpc)
  #tmpc = combos[1:50,] # need to change this
  ares = apply(combos, 1, calcscore, dtt=dt, target=target)
  res = do.call(rbind.data.frame, ares)
  res = as.data.table(res)
  names(res) = c("pair", "db", "dist.man")
  
  saveRDS(res, paste0("scripts/clustering-method/rdata/",gsub(" ", "_", target), fnpost, ".RData"))
  return (res);
}

calcClustMetricsFasterSingles <- function(target, am, fnpost, tm)
{
  dt = am[am$type == "tissue" | am$tissue.cancer == target,]
  fn =  file(paste0("scripts/clustering-method/sketches/", gsub(" ", "-", tolower(target)) ,"-sketch.txt"))
  sketch = readLines(fn)
  print(length(sketch))
  close(fn)
  
  dt = dt[dataset %in% sketch,]
  tm2 = data.frame(V1=tm[tm %in% names(dt)])

  ares = apply(tm2, 1, calcscoresingle, dtt=dt, target=target)
  res = do.call(rbind.data.frame, ares)
  res = as.data.table(res)
  names(res) = c("pair", "db", "dist.man", "ttype")
  
  saveRDS(res, paste0("scripts/clustering-method/rdata/",gsub(" ", "_", target), fnpost, ".RData"))
  return (res);
}

calcClustMetricsFasterClin <- function(target, am, fnpost)
{
  dt = am[am$type == "tissue" | am$tissue.cancer == target,]
  fn =  file(paste0("scripts/clustering-method/sketches/", gsub(" ", "-", tolower(target)) ,"-sketch.txt"))
  sketch = readLines(fn)
  print(length(sketch))
  close(fn)
  
  clinical = fread("data/clinical-all.txt", sep="\t")
  
  dt = dt[dataset %in% sketch,]
  tmpc = names(dt)[2:ncol(dt)-2]
  tmpc = tmpc[-1]
  combos = as.data.frame(t(combn(tmpc, 2)))
  rm(tmpc)
  
  combos = combos[combos$V1 %in% clinical$symbol | combos$V2 %in% clinical$symbol,]
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

calcscoresingle <- function(row, dtt, target) { # calc scores
  dtt = as.data.table(dtt)
  g = as.character(row["V1"])
  tmp = dtt[,c(g, "tissue.cancer"), with=F]
  tmp$target = ifelse(tmp$tissue.cancer == target, "target", "other")
  clust <- ifelse(tmp$target == "target", 1, 2)
  m = as.matrix(tmp[,1, with=FALSE])
  mode(m) = "numeric"
  #db = index.DB.fixed(m, as.integer(clust), centrotypes = "centroids", p=2)$DB
  db = index.DB.fixed(m, as.integer(clust), centrotypes = "centroids", p=1)$DB
  cm = colMeans(tmp[target == "other",1, with=F])
  cm = rbind(cm, colMeans(tmp[target == "target",1, with=F]))
  row.names(cm) = NULL
  cm = as.matrix(cm)
  cent = dist(cm, method="manhattan")[1]
  
  fl = "H"
  if (as.numeric(cm[1,]) > as.numeric(cm[2,]))
  {
    fl = "L"
  }
  
  return (list(pair=g, db, cent, fl))
}

# calcscore <- function(row, dtt, target) {
#   genes = c(as.character(row["V1"]), as.character(row["V2"]))
#   tmp = dtt[,c(genes[1], genes[2], "tissue.cancer"), with=F]
#   tmp$target = ifelse(tmp$tissue.cancer == target, "target", "other")
#   clust <- ifelse(tmp$target == "target", 1, 2)
#   m = as.matrix(tmp[,1:2, with=F])
#   mode(m) = "numeric"
#   db = index.DB(m, as.integer(clust), centrotypes = "centroids")$DB
#   cb = index.G1(m, as.integer(clust), centrotypes = "centroids")
#   #sil = silhouette(as.integer(clust), dist(m)) # use silhouette from cluster its faster, still slows down to less than 1 a sec
#   #sil = mean(sil[sil[,"cluster"] == 1,3])
#   
#   return (list(pair=paste(genes[1], genes[2], sep=":"), db, cb))#, sil))
# }



# version for tissueGPS - tissue only comparisons
calcClustMetricsFasterTG <- function(target, am, fnpost)
{
  dt = am
  tmpc = names(dt)[2:ncol(dt)-2]
  tmpc = tmpc[-1]
  combos = as.data.frame(t(combn(tmpc, 2)))
  rm(tmpc)
  #tmpc = combos[1:50,] # need to change this
  ares = apply(combos, 1, calcscore, dtt=dt, target=target)
  res = do.call(rbind.data.frame, ares)
  res = as.data.table(res)
  names(res) = c("pair", "db", "ch")
  
  saveRDS(res, paste0("scripts/clustering-method/rdata/",gsub(" ", "_", target), fnpost, ".RData"))
  return (res);
}




## testing the calcscore script
#system.time({ares = apply(tmpc, 1, calcscore)}) #181.472 for 50 
#res.df = do.call(rbind.data.frame, ares)
#res.df = as.data.table(res.df)

#system.time({blah = apply(tmpc, 1, calcscoreNew)}) #0.7 for 50
#res.blah = do.call(rbind.data.frame, blah)
#res.blah = as.data.table(res.blah)


#calcscore <- function(row) {
#  genes = c(as.character(row["V1"]), as.character(row["V2"]))
#  tmp = dt[,c(genes[1], genes[2], "tissue.cancer"), with=F]
#  tmp$target = ifelse(tmp$tissue.cancer == target, "target", "other")
#  clust <- ifelse(tmp$target == "target", 1, 2)
#  m = as.matrix(tmp[,1:2, with=F])
#  mode(m) = "numeric"
#  ic = intCriteria(m, as.integer(clust), c("sd_dis", "sd_scat", "davies_bouldin", "calinski_harabasz", "silhouette"))
#  
#  return (list(pair=paste(genes[1], genes[2], sep=":"), sd=ic$sd_dis, sd.scat=ic$sd_scat, db=ic$davies_bouldin, ch=ic$calinski_harabasz, sil=ic$silhouette))
#}

## for the new matrix script need to speed up distance matrix computation
#system.time({blah = dist(m)}) #0.611
#system.time({sil = index.S(blah,as.integer(clust))}) # 9.5 seconds
#system.time({sil = silhouette(as.integer(clust), blah)}) # 0.327 seconds
#system.time({sil = silhouette(as.integer(clust), dist(m))}) # 0.7 seconds

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


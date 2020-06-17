# ***********************************
# functions for clustering scores
# ***********************************

# function called by calcClustMetricsForSingles
# does the work of calculating the clustering
# scores for a single gene
calcscoresingle <- function(row, dtt, target) 
{
  dtt = as.data.table(dtt)
  g = as.character(row["V1"]) # get the gene
  
  # get cluster values for a single gene
  tmp = dtt[,c(g, "tissue.cancer"), with=F]
  tmp$target = ifelse(tmp$tissue.cancer == target, "target", "other")
  clust <- ifelse(tmp$target == "target", 1, 2)
  m = as.matrix(tmp[,1, with=FALSE])
  mode(m) = "numeric"
  db = index.DB.fixed(m, as.integer(clust), centrotypes = "centroids", p=1)$DB
  cm = colMeans(tmp[target == "other",1, with=F])
  cm = rbind(cm, colMeans(tmp[target == "target",1, with=F]))
  row.names(cm) = NULL
  cm = as.matrix(cm)
  cent = dist(cm, method="manhattan")[1]
  
  # get H/L type
  fl = "H"
  if (as.numeric(cm[1,]) > as.numeric(cm[2,]))
  {
    fl = "L"
  }
  
  return (list(pair=g, db, cent, fl))
}

# function called by calcClustMetrics
# does the work of calculating the clustering
# scores for an antigen pair
calcscore <- function(row, dtt, target)
{ 
  dtt = as.data.table(dtt)
  genes = c(as.character(row["V1"]), as.character(row["V2"]))
  
  # get cluster values per pair
  tmp = dtt[,c(genes[1], genes[2], "tissue.cancer"), with=F]
  tmp$target = ifelse(tmp$tissue.cancer == target, "target", "other")
  clust <- ifelse(tmp$target == "target", 1, 2)
  m = as.matrix(tmp[,1:2, with=FALSE])
  mode(m) = "numeric"
  db = index.DB.fixed(m, as.integer(clust), centrotypes = "centroids", p=1)$DB
  cm = colMeans(tmp[target == "other",1:2, with=F])
  cm = rbind(cm, colMeans(tmp[target == "target",1:2, with=F]))
  row.names(cm) = NULL
  cm = as.matrix(cm)
  cent = dist(cm, method="manhattan")[1]
  
  # get type (HH, HL / LH, LL)
  fl = "H"
  sl = "H"
  if (as.numeric(cm[1,1]) > as.numeric(cm[2,1]))
  {
    fl = "L"
  }
  if (as.numeric(cm[1,2]) > as.numeric(cm[2,2]))
  {
    sl = "L"
  }
  
  return (list(pair=paste(genes[1], genes[2], sep=":"), db, cent, paste0(fl, sl)))
}

# function called by calcClustMetricsTriples
# does the work of calculating the clustering
# scores for an antigen pair
calcscoretriplets <- function(row, dtt, target) 
{
  dtt = as.data.table(dtt)
  genes = c(as.character(row["V1"]), as.character(row["V2"]), as.character(row["V3"]))
  
  # get cluster values per triple
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
  
  # get type (HHH, LLL, HLL, LHH, etc)
  fl = "H"
  sl = "H"
  en = "H"
  if (as.numeric(cm[1,1]) > as.numeric(cm[2,1]))
  {
    fl = "L"
  }
  if (as.numeric(cm[1,2]) > as.numeric(cm[2,2]))
  {
    sl = "L"
  }
  if (as.numeric(cm[1,3]) > as.numeric(cm[2,3]))
  {
    en = "L"
  }
  
  return (list(pair=paste(genes[1], genes[2], genes[3], sep=":"), db, cent, paste0(fl, sl, en)))
}

# takes a cancer (target) and the entire expression matrix
# returns a data table of clustering scores for every gene's
# power to separate target from all other normal tissue samples
# uses only the samples in the sketch
# db = Davies Bouldin, dist.man = Manhattan distance b/w the 
# cluster centroids, ttype = whether gene is on (H) or off (L) in
# the target
calcClustMetricsForSingles <- function(target, mat)
{
  dt = mat[mat$type == "tissue" | mat$tissue.cancer == target,]
  fn =  file(paste0("results/sketches/", gsub(" ", "-", tolower(target)) ,"-sketch.txt"))
  sketch = readLines(fn)
  close(fn)
  
  dt = dt[dataset %in% sketch,]
  genes = names(dt)
  genes = genes[!(genes %in% c("type", "tissue.cancer", "batch", "dataset"))]
  genes = data.frame(V1=genes)
  
  ares = apply(genes, 1, calcscoresingle, dtt=dt, target=target)
  res = do.call(rbind.data.frame, ares)
  res = as.data.table(res)
  names(res) = c("combo", "db", "dist.man", "ttype")
  
  return (res);
}

# takes a cancer (target) and the entire expression matrix
# returns a data table of clustering scores for every possible antigen pair's
# ability to separate target samples from all other normal tissue samples
# uses only the samples in the sketch
# db = Davies Bouldin, dist.man = Manhattan distance b/w the 
# cluster centroids, ttype = whether gene is on (H) or off (L) in
# the target
calcClustMetrics <- function(target, mat, testRun=F)
{
  dt = mat[mat$type == "tissue" | mat$tissue.cancer == target,]
  fn =  file(paste0("results/sketches/", gsub(" ", "-", tolower(target)) ,"-sketch.txt"))
  sketch = readLines(fn)
  close(fn)
  
  dt = dt[dataset %in% sketch,]
  tmpc = names(dt)[2:ncol(dt)-2]
  tmpc = tmpc[-1]
  combos = as.data.frame(t(combn(tmpc, 2)))
  rm(tmpc)
  
  if(testRun)
  {
    combos = combos[1:50,] 
  }
  
  ares = apply(combos, 1, calcscore, dtt=dt, target=target)
  res = do.call(rbind.data.frame, ares)
  res = as.data.table(res)
  names(res) = c("combo", "db", "dist.man", "ttype")
  
  return (res);
}

# takes a cancer (target) and the entire expression matrix
# returns a data table of clustering scores for an antigen triplet's
# ability to separate target samples from all other normal tissue samples
# uses only the samples in the sketch
# the input singles data narrows down the total number of triples
# by keeping only ones that show some minimum potential per tumor type
# db = Davies Bouldin, dist.man = Manhattan distance b/w the 
# cluster centroids, ttype = whether gene is on (H) or off (L) in
# the target
calcClustMetricsTriples <- function(target, mat, singles, testRun=F)
{
  dt = mat[mat$type == "tissue" | mat$tissue.cancer == target,]
  fn =  file(paste0("results/sketches/", gsub(" ", "-", tolower(target)) ,"-sketch.txt"))
  sketch = readLines(fn)
  close(fn)
  
  dt = dt[dataset %in% sketch,]
  sn = singles[can==target,]
  combos = as.data.frame(t(combn(as.character(sn$combo), 3)))
  
  if(testRun)
  {
    combos = combos[1:50,] 
  }
  
  ares = apply(combos, 1, calcscoretriplets, dtt=dt, target=target)
  res = do.call(rbind.data.frame, ares)
  res = as.data.table(res)
  names(res) = c("combo", "db", "dist.man", "ttype")
  
  return (res);
}


# taken from the ClusterSim R package
# but their davies bouldin implementation had
# an error - copied their function and fixed the bug
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

# *****************************
#       Eval functions
# *****************************

# takes a tumor type (cn), genes, expression data and
# returns evaluation metrics using decision trees
assignClusterLabels <- function(cn, genes, dt)
{
  # get the sketch subset
  fn =  file(paste0("results/sketches/", gsub(" ", "-", tolower(cn)) ,"-sketch.txt"))
  sketch = readLines(fn)
  close(fn)
  
  train = dt[dataset %in% sketch,] # training set
  train = train[train$type == "tissue" | train$tissue.cancer == cn,]
  train[, "target" := ifelse(train$tissue.cancer == cn, "target", "other")]
  
  test = dt[!(dataset %in% sketch),] # test set for eval
  test = test[test$type == "tissue" | test$tissue.cancer == cn,]
  test[, "target" := ifelse(test$tissue.cancer == cn, "target", "other")]
  
  if(length(genes) == 1) # single gene case
  {
    train = train[,c(genes[1], "tissue.cancer", "target"), with =F]
    test = test[,c(genes[1], "tissue.cancer", "target"), with=F]
    
    test[, clustlab := assignDT(test[,1, with=F], train, T)] # run the DTs
  }
  
  if (length(genes) == 2) # double
  {
    train = train[,c(genes[1], genes[2], "tissue.cancer", "target"), with =F]
    test = test[,c(genes[1], genes[2], "tissue.cancer", "target"), with=F]
    
    test[, clustlab := assignDT(test[,1:2, with=F], train, F)] # run the DTs
  }
  
  if (length(genes) == 3) # triple
  {
    train = train[,c(genes[1], genes[2], genes[3], "tissue.cancer", "target"), with =F]
    test = test[,c(genes[1], genes[2], genes[3], "tissue.cancer", "target"), with=F]
    
    test[, clustlab := assignDT(test[,1:3, with=F], train, F)]
  }
  
  # calc performance metrics  
  tp = length(test[target == "target" & clustlab == "target",]$clustlab)
  tn = length(test[target == "other" & clustlab == "other",]$clustlab)
  fp = length(test[target == "other" & clustlab == "target",]$clustlab)
  
  prec = tp / (tp + fp) # precision = tp / (tp + fp)
  rec = tp / length(test[target=="target",]$clustlab) # recall = tp / total number of positives (fn) 
  f1 = (2 * prec * rec) / (prec + rec)
  
  return(list(F1=f1, prec=prec, recall=rec))
}

# use decision tree boundaries on the training data
# to evaluate performance on the test set
assignDT <- function(test, train, single)
{
  ctmp = train
  if(single) # if single gene
  {
    names(ctmp) = c("g1", "tissue.cancer", "target")
    ctmp[, y := ifelse(target == "target", 1, -1)]
    ctmp$y = as.factor(ctmp$y)
    
    names(test) = c("g1")
    
    dt = rpart(y ~ g1, method="class", data=ctmp[,c("g1", "y")], cp=-1, maxdepth=1)
    test$prediction = predictOnSplits(test$g1, dt$splits[4], dt$splits[2])   
  }
  else
  {
    if (length(names(ctmp)) == 4) # double
    {
      names(ctmp) = c("g1", "g2", "tissue.cancer", "target")
      ctmp[, y := ifelse(target == "target", 1, -1)]
      ctmp$y = as.factor(ctmp$y)
      
      names(test) = c("g1", "g2")
      
      dt = rpart(y ~ g1, method="class", data=ctmp[,c("g1", "y")], cp=-1, maxdepth=1)
      test$pg1 = predictOnSplits(test$g1, dt$splits[4], dt$splits[2])
      
      dt = rpart(y ~ g2, method="class", data=ctmp[,c("g2", "y")], cp=-1, maxdepth=1)
      test$pg2 = predictOnSplits(test$g2, dt$splits[4], dt$splits[2])
      
      test[, prediction := pmin(pg1, pg2)]
      
    }
    else # triple
    {
      names(ctmp) = c("g1", "g2", "g3", "tissue.cancer", "target")
      ctmp[, y := ifelse(target == "target", 1, -1)]
      ctmp$y = as.factor(ctmp$y)
      
      names(test) = c("g1", "g2", "g3")
      
      dt = rpart(y ~ g1, method="class", data=ctmp[,c("g1", "y")], cp=-1, maxdepth=1)
      test$pg1 = predictOnSplits(test$g1, dt$splits[4], dt$splits[2])
      
      dt = rpart(y ~ g2, method="class", data=ctmp[,c("g2", "y")], cp=-1, maxdepth=1)
      test$pg2 = predictOnSplits(test$g2, dt$splits[4], dt$splits[2])
      
      dt = rpart(y ~ g3, method="class", data=ctmp[,c("g3", "y")], cp=-1, maxdepth=1)
      test$pg3 = predictOnSplits(test$g3, dt$splits[4], dt$splits[2])
      
      test[, prediction := pmin(pg1, pg2, pg3)]
      
    }
  }
  clusters <- ifelse(test$prediction >= 1, "target", "other")
  return (clusters)
}

# use the cutoff line to predict
# label of tumor vs normal
predictOnSplits <- function(x, cutoff, ncat)
{
  if (is.null(cutoff))
  {
    return(NA)
  }
  pred = ifelse(x < cutoff, 1, -1)
  if (ncat == -1) # greater than
  {
    pred = ifelse(x > cutoff, 1, -1)
  }
  return (pred)
}

# *****************************
#       Analysis functions
# *****************************

# takes a data.table of scores and scales both the DB 
# and distance into a single clustering score that ranges
# from 0 to 1 based on the distribution in the data.table
# returns a data.table now sorted by these values
cleanAndScale <- function(dt)
{
  dt = dt[!is.na(db)] # remove NA values
  dt[,dbscaled := 1-((log(db) - min(log(db)))/(max(log(db)) - min(log(db))))] # make score between 0 and 1
  dt[,distscaled := (log(dist.man) - min(log(dist.man)))/(max(log(dist.man)) - min(log(dist.man)))]
  dt[,score := apply(dt[,c("dbscaled", "distscaled"), with=F], 1, min)]
  
  return (dt[order(-score),]) 
}





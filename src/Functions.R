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
  
  fl = "H"
  if (as.numeric(cm[1,]) > as.numeric(cm[2,]))
  {
    fl = "L"
  }
  
  return (list(pair=g, db, cent, fl))
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


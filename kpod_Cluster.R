#####Edits to kpod method to increase number of random starts and return sum of within group SS
#####Last updated 1 Aug 2017 by Noah Lottig (nrlottig@wisc.edu)
require(cluster)

#' Make test data
#' 
#' \code{makeData} Function for making test data
#' 
#' @param p Number of features (or variables)
#' @param n Number of observations
#' @param k Number of clusters
#' @param sigma Variance
#' @param missing Desired missingness percentage
#' @param seed (Optional) Seed (default seed is 12345)
#' 
#' @export
#' 
#' @examples
#' p <- 2
#' n <- 100
#' k <- 3
#' sigma <- 0.25
#' missing <- 0.05
#' 
#' X <- makeData(p,n,k,sigma,missing)$Orig
#' 
#' @author Jocelyn T. Chi
makeData <- function(p,n,k,sigma,missing,seed=12345){
  if(p <= 0 | n <= 0 | k <= 0){
    return(cat('Please select positive values for p, n, and k.'))
  }
  if (missing < 0 | missing > 1){
    return(cat('Please select a missingness percentage between 0 and 1.'))
  }
  
  # Make complete data
  set.seed(seed)
  M <- matrix(rnorm(k*p),k,p)
  assignment <- sample(1:k,n,replace=TRUE)
  X <- M[assignment,] + sigma*matrix(rnorm(n*p),n,p)
  
  # Make missing data
  X_missing <- X
  missing_ix <- sample(1:(n*p),(n*p*missing),replace=TRUE)
  X_missing[missing_ix] <- NA
  
  return(list(Orig=X,Missing=X_missing,truth=assignment))  
}

#' Function for assigning clusters to rows in a matrix
#' 
#' \code{assign_clustpp} Function for assigning clusters to rows in a matrix
#' 
#' @param X Data matrix containing missing entries whose rows are observations and columns are features
#' @param init_centers Centers for initializing k-means
#' @param kmpp_flag (Optional) Indicator for whether or not to initialize with k-means++
#' @param max_iter (Optional) Maximum number of iterations
#' 
#' @export
#' 
#' @examples
#' p <- 2
#' n <- 100
#' k <- 3
#' sigma <- 0.25
#' missing <- 0.05
#' Data <- makeData(p,n,k,sigma,missing)
#' X <- Data$Missing
#' Orig <- Data$Orig
#' 
#' clusts <- assign_clustpp(Orig, k)
#' 
#' @author Jocelyn T. Chi
assign_clustpp <- function(X,init_centers,kmpp_flag=TRUE,max_iter=10000){
  res <- kmeans(X, init_centers,iter.max = 10000,nstart = 200)
  clusts <- res$cluster
  obj <- res$totss
  # fit <- 1-(sum(res$withinss)/res$totss)
  fit <- sum(res$withinss)
  centers <- res$centers
  m.sil = silhouette(res$cluster, dist(X))
  m.sil = mean(m.sil[,3])
  if (kmpp_flag == TRUE) {
    ## Try to find a better assignment
    for (iter in 1:max_iter) {
      centers_kmpp <- kmpp(X,length(res$size))
      sol <- kmeans(X, centers_kmpp,iter.max = 10000,nstart = 200)
      if (sol$totss < obj) {
        obj <- sol$totss
        clusts <- sol$cluster
        # fit <- 1-(sum(sol$withinss)/sol$toss)
        fit <- sum(sol$withinss)
        centers <- sol$centers
        m.sil = silhouette(res$cluster, dist(X))
        m.sil = mean(m.sil[,3])
        break
      }
    }
  }
  return(list(clusts=clusts,obj=obj,centers=centers,fit=fit,silh=m.sil))
}

#' Function for finding indices of missing data in a matrix
#' 
#' \code{findMissing} Function for finding indices of missing data in a matrix
#' 
#' @param X Data matrix containing missing entries whose rows are observations and columns are features
#' 
#' @return A numeric vector containing indices of the missing entries in X
#' 
#' @export
#' 
#' @examples
#' p <- 2
#' n <- 100
#' k <- 3
#' sigma <- 0.25
#' missing <- 0.05
#' Data <- makeData(p,n,k,sigma,missing)
#' X <- Data$Missing
#' missing <- findMissing(X)
#' 
#' @author Jocelyn T. Chi
findMissing <- function(X){
  missing_all <- which(is.na(X))
  return(missing_all)
}

#' Function for initial imputation for k-means
#' 
#' \code{initialImpute} Initial imputation for k-means
#' 
#' @param X Data matrix containing missing entries whose rows are observations and columns are features
#' 
#' @return A data matrix containing no missing entries
#' 
#' @export
#' 
#' @examples
#' p <- 2
#' n <- 100
#' k <- 3
#' sigma <- 0.25
#' missing <- 0.05
#' Data <- makeData(p,n,k,sigma,missing)
#' X <- Data$Missing
#' X_copy <- initialImpute(X)
#' 
#' @author Jocelyn T. Chi
initialImpute <- function(X){
  avg <- mean(X,na.rm=TRUE)
  X[which(is.na(X))] <- avg
  return(X)
}

#' Function for performing k-POD
#' 
#' \code{kpod} Function for performing k-POD, a method for k-means clustering on partially observed data
#' 
#' @param X Data matrix containing missing entries whose rows are observations and columns are features
#' @param k Number of clusters
#' @param kmpp_flag (Optional) Indicator for whether or not to initialize with k-means++
#' @param maxiter (Optional) Maximum number of iterations
#' 
#' @return cluster: Clustering assignment obtained with k-POD
#' @return cluster_list: List containing clustering assignments obtained in each iteration
#' @return obj_vals: List containing the k-means objective function in each iteration
#' @return fit: Fit of clustering assignment obtained with k-POD (calculated as total withinss))
#' @return fit_list: List containing fit of clustering assignment obtained in each iteration
#' 
#' @export
#' 
#' @import clues
#' 
#' @examples
# p <- 5
# n <- 200
# k <- 3
# sigma <- 0.15
# missing <- 0.20
# Data <- makeData(p,n,k,sigma,missing)
# X <- Data$Missing
# Orig <- Data$Orig
# truth <- Data$truth
# 
# kpod_result <- kpod(X,k)
# kpodclusters <- kpod_result$cluster
#' 
#' @author Jocelyn T. Chi
#' 
kpod <- function(X,k,kmpp_flag=TRUE,maxiter=1000){
  
  n <- nrow(X)
  p <- ncol(X)
  
  cluster_vals <- vector(mode="list",length=maxiter)
  obj_vals <- double(maxiter)
  fit <- double(maxiter)
  m.sil <- double(maxiter)
  sil_val <- vector(mode="list",length=maxiter)
  
  missing <- findMissing(X)
  
  # Run first iteration
  X_copy <- initialImpute(X)
  
  ## Use kmpp to select initial centers
  init_centers <- kmpp(X_copy, k)
  temp <- kmeans(X_copy,init_centers,iter.max = 1000,nstart = 200)
  clusts <- temp$cluster
  centers <- temp$centers
  # fit[1] <- 1-(sum(temp$withinss)/temp$totss)
  fit[1] <- sum(temp$withinss)        
  t.sil = silhouette(temp$cluster, dist(X_copy))
  m.sil[1] = mean(t.sil[,3])
  sil_val[[1]] = t.sil[,3]

  # Make cluster matrix
  clustMat <- centers[clusts,]
  
  # Update with centers
  X_copy[missing] <- clustMat[missing]
  
  #obj_vals[1] <- temp$obj
  obj_vals[1] <- sum((X[-missing]-clustMat[-missing])^2)
  cluster_vals[[1]] <- clusts
  # Run remaining iterations
  for (i in 2:maxiter){
    temp <- assign_clustpp(X_copy,centers,kmpp_flag)
    clusts <- temp$clusts
    centers <- temp$centers
    fit[i] <- temp$fit
    t.sil = silhouette(temp$clusts, dist(X_copy))
    m.sil[i] = mean(t.sil[,3])
    sil_val[[i]] = t.sil[,3]
    
    # Impute clusters
    clustMat <- centers[clusts,]
    X_copy[missing] <- clustMat[missing]
    
    obj_vals[i] <- sum((X[-missing]-clustMat[-missing])**2)
    cluster_vals[[i]] <- clusts
    message(paste('Number of kmeans iterations =',i))
    if (all(cluster_vals[[i]] == cluster_vals[[i-1]])){
      message('Clusters have converged.')
      return(list(cluster=clusts,cluster_list=cluster_vals[1:i],obj_vals=obj_vals[1:i],fit=fit[i],fit_list=fit[1:i],silh=m.sil[i],silh_list=m.sil[1:i],silh_vals=sil_val[i]))
      break
    }
  }
  return(list(cluster=clusts,cluster_list=cluster_vals[1:i],obj_vals=obj_vals[1:i],fit=fit[i],fit_list=fit[1:i],silh=m.sil[i],silh_list=m.sil[1:i],silh_vals=sil_val[i]))
}

#' k-means++
#'
#' \code{kmpp} Computes initial centroids via kmeans++
#'
#' @param X Data matrix whose rows are observations and columns are features
#' @param k Number of clusters.
#' 
#' @return A data matrix whose rows contain initial centroids for the k clusters
#' 
#' @export
#' 
#' @examples
#' n <- 10
#' p <- 2
#' X <- matrix(rnorm(n*p),n,p)
#' k <- 3
#' kmpp(X,k)
#' 
kmpp <- function(X, k) {
  n <- nrow(X)
  p <- ncol(X)
  C <- integer(k)
  C[1] <- sample(1:n, 1)
  for (i in 2:k) {
    S <- matrix(NA,n,i-1)
    for (j in 1:(i-1)) {
      S[,j] <- apply(X -
                       matrix(X[C[j],],n,p,byrow=TRUE),1,FUN=function(x)
                       {norm(as.matrix(x),'f')**2})
    }
    D <- apply(S,1,min)
    pr <- D/sum(D)
    C[i] <- sample(1:n, 1, prob = pr)
  }
  return(X[C,])
}

######Function to estimate random Cluster Within SS for determining optimal number of clusters

quant_no_clusters = function(X,no_clusters=2:7,boots=100){
  #data should be standardized, rows are sites, columns are observations
  require(GMD)
  n.cluster = length(no_clusters)
  rand_fit = rep_len(NA,n.cluster)
  best_fit = rep_len(NA,n.cluster)
  sil_fit = rep_len(NA,n.cluster)
  sil_val <- vector(mode="list",length=n.cluster)
  cluster_vals <- vector(mode="list",length=n.cluster)
  pb <- txtProgressBar(min = 1, max = max(no_clusters), style = 3)
  for (j in 1:n.cluster){
    temp.fit = rep_len(NA, boots)
    for(i in 1:boots){
      set.seed(NULL)
      X_copy = initialImpute(X)
      temp = css(dist(X_copy),sample(x = c(1:no_clusters[j]),replace = TRUE,size = nrow(X_copy)))
      
          # tryCatch({
        # seg=NA
        # (temp <- kmeans(X_copy,no_clusters[j],iter.max = 1,nstart = 1))
      # },error=function(e){})
      
      temp.fit[i] = temp$totwss
    }
    rand_fit[j] = mean(temp.fit)
    kpod_fit = kpod(X,no_clusters[j],kmpp_flag = TRUE)
    best_fit[j] = kpod_fit$fit
    sil_fit[j] = kpod_fit$silh
    sil_val[[j]] = kpod_fit$silh_vals
    cluster_vals[[j]] = kpod_fit$cluster
    setTxtProgressBar(pb, no_clusters[j])
  }
  close(pb)
  
  for(j in 1:n.cluster){
  sil_sum1 = data.frame(cluster_size=rep(no_clusters[j],nrow(X)),Sil_values=sil_val[[j]][[1]],cluster_ID=cluster_vals[[j]])
  if(j==1) sil_sum = sil_sum1
  if(j>1) sil_sum = rbind(sil_sum,sil_sum1)
  }
  
  sum.sil =aggregate(sil_sum$Sil_values,by=list(sil_sum$cluster_size,sil_sum$cluster_ID),FUN=mean)
  sum.silmedian = aggregate(sil_sum$Sil_values,by=list(sil_sum$cluster_size),FUN=median)
  
  withinss_range = range(rand_fit,best_fit)
  par(mfrow=c(4,1),mar=c(2.2,4,0,0),oma=c(3,0,.25,0.25))
  plot(no_clusters,best_fit,type="b",col="blue",ylab="Total Within SS",xlab="")
  plot(no_clusters,rand_fit,ylim=withinss_range,type="b",col="red",ylab="Total Within SS",xlab="")
  lines(no_clusters,best_fit,type="b",col="blue")
  legend('topright',legend=c("Best Fit","Random Fit"),lty=1,col=c("blue","red"))
  plot(no_clusters,(rand_fit-best_fit),type="b",ylab="Difference between Random & Best",xlab="")
  plot(sil_sum[,1:2],col="grey",ylab = "Silhouette Coefficient",xlab ="")
  points(no_clusters,sil_fit,col="green",pch=16,cex=2,type="b")
  points(sum.silmedian$Group.1,sum.silmedian$x,col="red",pch=16,cex=2,type="b")
  points(sum.sil$Group.1,sum.sil$x,col="black",pch=16,cex=1)
  mtext(side=1,line=2.5,text = "Number of Clusters",cex=0.75)
  # plot(Cluster_Vec,(best_fit-rand_fit),type="l”)
  
  return(list(rand_SS=rand_fit,kpod_SS=best_fit,msil_val=sil_fit,sil_vals=sil_val))
}

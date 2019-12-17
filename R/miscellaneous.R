##############
# miscellaneous functions for DESeq2 implementation

# Function nbinomLogLike: calculate sum of log likelihood of counts data, given dispersion parameter and mean
nbinomLogLike = function(countsdata,mu,dispersions){
  
  n = nrow(countsdata) # number of genes
  m = ncol(countsdata) # number of samples
  ll = 0
  
  for (i in 1:n){
    alpha = dispersions[i]
    for (j in 1:m){
      ll = ll + dnbinom(x = countsdata[i,j], size = alpha, mu = mu[i,j], log=T)
    }
  }
  return(ll)
}

# Function getBaseMeanAndVariances

getBaseMeansAndVariances <- function(object) {
  
  baseMean = unname(rowMeans(counts(object,normalized=TRUE)))
  baseVar = unname(apply(counts(object,normalized=TRUE),1,var))
  allZero = unname(rowSums(counts(object)) == 0)
  
  object[["baseMean"]] = baseMean
  object[["baseVar"]] = baseVar
  object[["allZero"]] = allZero
  return(object) 
}



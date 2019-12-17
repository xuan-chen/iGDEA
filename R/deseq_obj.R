##############
# supplementary functions for initializing, querying NB-GLM objects, which is a list
# source("main.R")

# Function deseq_data: initialize a list object for metadata and data storage
#' @param data counts data, n by m, n: number of genes, m: number of samples
#' @param designMatrix the matrix containing the treatment condition for each sample, m by k, m: number of samples, k: covariates
#' @param nf user-supplied normalization factor
#' @return a list, the NBGLM object, containing countdata, designMatrix information and normalization factors

deseq_data = function(data,designMatrix,nf=NULL){
  object = list(counts = data, design = designMatrix, nf = nf, dispersions = NULL, dispGeneEst = NULL, dispGeneIter = NULL)
  return(object)
}

# Function counts: take out the counts data for the deseq_data object and check if it is matrix
#' @param object the NBGLM object
#' @param normalized boolean, whether one wants the normalized count, normalization is done using the definition of s_{ij} in p.14 of the references
#' @param idx row index, which can be used to subset the count data by row
#' @return count data

counts = function(object,normalized=FALSE,idx = NULL){
  countdata = object$counts
  if (!is.matrix(countdata) & is.data.frame(countdata)){
     countdata = data.matrix(countdata)
  }
  
  if (normalized){
    nf = getNormFactors(object, matformat = T)
    countdata = countdata/nf 
  }
  
  if(is.null(idx)){
    return(countdata)
  }else{
    return(countdata[idx,])
  }
}


# Functin getmodelMatrix: query the model matrix (the design matrix containing the treatment information of each sample)
#' @param object the NBGLM object
#' @return the design matrix of the object

getmodelMatrix = function(object){
  return(object$design)
}

# Function getNormFactors: get the normalization factor of each sampel
#' @param object the deseq_data object
#' @param matformat boolean, whether one needs an n by m normalization factor matrix instea of a numeric vector of length m
#' @return normalization factor vector/matrix

getNormFactors = function(object,matformat=T){
  countdata = counts(object,normalized = F)
  genes = nrow(countdata)
  
  if (!is.null(object$nf)){
    nf = object$nf
    newnf = nf
    if (length(dim(nf))==2 & matformat==F){
      newnf = as.vector(nf[1,]) 
      #object$nf <<- newnf
    }else if(is.vector(nf) & matformat==T){
      newnf = matrix(rep(nf,each=genes),genes,length(nf))
      #object$nf <<- newnf
    }
    return(newnf)
  }
  
  nf = estimateSizeFactors(object)
  if(matformat==T){
    newnf = matrix(rep(nf,each=genes),genes,length(nf))
  }else{
    newnf = nf
  }
  #object$nf <<- newnf
  return(newnf)
}


# Function dispersions: get the dispersion parameter of the deseq_data object
#' @param object the NBGLM object
#' @return the dispersion parameter vector for all genes

dispersions = function(object){
  if(!is.null(object$dispersions)){
    return(object$dispersions)
  }else{
    stop("dispersions of DESeq2 data object is null") # need more coding here!
  }
}

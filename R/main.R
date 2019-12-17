####################
# Project for Biostats 615: A re-implementation of package "DESeq2"
# for gene differential expression analysis (DEA)
####################

source('./R/disp_sup.R')
source('./R/fromDESeq2.R')
source("./R/deseq_obj.R")
source("./R/miscellaneous.R")
source("./R/wald_sup.R")

source("./R/Create_DE_object.R")
source("./R/Filter_Data.R")
source("./R/Normalize_Data.R")
source("./R/DE_test.R")
source("./R/Get_DEG.R")



library(locfit)

# functions for doing DEA:
# estimateSizeFactors --- estimate the size factor for each sample j
# estimateDispersionGeneEst --- estimate the dipersion parameter for each gene i, which is the gene-wise dispersion parameter
# estimateDispersionFit --- A trend fitting function to fit all gene-wise dispersion parameter and find their prior mean
# estimateDispersionMAP --- Finding the maximum a posteriori estimate of each gene's dispersion parameter
# nbinomWaldTest --- use NB-GLM model to estimate the beta for each treatment condition, given dispersion parameter for each gene, and test the
#                    significance of each beta using Wald Test
# DEA_deseq --- A function that integrates all of the functions above into an organized workflow, which is the function to call

#' @references
#'
#' DESeq2 reference:
#'
#' Love, M.I., Huber, W., Anders, S. (2014) Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology, 15:550. \url{http://dx.doi.org/10.1186/s13059-014-0550-8}
#'

# function for estimating the size factor of each sample j
#' @param object our own deseq2 data object
estimateSizeFactors = function(object){
  countsdata= counts(object,normalized = FALSE) # get the counts data
  nozero = which(apply(countsdata,1,function(x) ifelse(sum(x==0)>0,1,0))==0) # check if a row contains zero
  loggeoMeans = rowMeans((log(countsdata))) # the logarithm of the geometric mean of each gene, the K_i^R
  res = apply(countsdata,2,function(cnts) {exp(stats::median(log(cnts[nozero])-loggeoMeans[nozero]))}) # the size factor for each sample
  #object$nf <<- res # load the normalization factors into the object
  return(res)
}


# Function estimateDispersionGeneEst: estimate each gene's dispersion parameter, namely the alpha_i^{GW} on p.15 of DESeq2 paper
#' @param object: the deseq_data object
#' @param minDisp: the minimum dispersion parameter allowed
#' @param kappa_0: the initial step size for backtracking search
#' @param dispTol: the tolerance level of convergence for the backtracking search
#' @param maxit: the maximum iteration of backtracking search for each gene's dispersion estimation
#' @param useCR: boolean, indicator of whether using Cox-Reid adjustment or not
#' @param modelMatrix: the design matrix, m by k
#' @param niter: number of iterations where each iteration consist of estimating the mean based on dispersion and estimate dispersion based on fitted mean
#' @param minmu: minimum fitted mean allowed
#' @param linearMu: boolean, whether one uses linear model to estimate the mean of counts data or negative binomial GLM
#' @return object: contains the gene-wise dispersion estimate

estimateDispersionGeneEst = function(object, minDisp=1e-8, kappa_0=1, dispTol=1e-6, maxit=100, useCR=TRUE,
                                     modelMatrix=NULL, niter=1,minmu=0.5, linearMu = FALSE){
  if (log(minDisp/10) <= -30) {
    stop("for computational stability, log(minDisp/10) should be above -30") # this criterion is copied from DESeq2 source code
  }

  if (is.null(modelMatrix)){
    modelMatrix = getmodelMatrix(object)
  }
  countsdata = counts(object,normalized = FALSE)

  # now to kick off the iteration of NB-GLM fitting, we need to give a starting point for the dispersion parameter
  # the following block of codes are taken from the official package of DESeq2, with only minor modification
  roughDisp <- roughDispEstimate(y = counts(object,normalized=TRUE),x = modelMatrix)
  momentsDisp <- momentsDispEstimate(object)
  alpha_hat <- pmin(roughDisp, momentsDisp) # initial guess of alphas
  maxDisp <- max(10, ncol(countsdata))
  alpha_hat <- alpha_hat_new <- alpha_init <- pmin(pmax(minDisp, alpha_hat), maxDisp)
  stopifnot(length(niter) == 1 & niter > 0)

  fitidx <- rep(TRUE,nrow(countsdata))
  mu <- matrix(0, nrow=nrow(countsdata), ncol=ncol(countsdata))
  dispIter <- numeric(nrow(countsdata))

  normalizationfactors = getNormFactors(object,matformat = T)

  for(iter in seq_len(niter)){
    newobject = deseq_data(data = countsdata[fitidx,], designMatrix = modelMatrix, nf = normalizationfactors[fitidx,])
    if(linearMu){
      fitMu = linearModelMuNormalized(newobject,x = modelMatrix)
    }else{
      fit = fitNBinomGLM(newobject,x = modelMatrix,alpha_hat = alpha_hat[fitidx],useQR = T)
      fitMu = fit$mu # get the fitted mean, ready for estimating dispersion
    }


    fitMu[fitMu<minmu] = minmu
    mu[fitidx,] = fitMu

    # with the fitted mean value, we estimate the dispersion parameter
    dispRes = fitDisp(countsdata[fitidx,], modelMatrix, fitMu, log(alpha_hat)[fitidx], log_alpha_mean = log(alpha_hat)[fitidx],
                      log_alpha_s2 = 1, min_log_alpha = log(minDisp/10),
                      kappa_0 = kappa_0, tol = dispTol, maxit = maxit, use_prior = F, useCR = useCR)

    dispIter[fitidx] = dispRes$iter
    alpha_hat_new[fitidx] = pmin(exp(dispRes$log_alpha),maxDisp)
    fitidx = abs(log(alpha_hat_new) - log(alpha_hat)) > 0.05
    alpha_hat = alpha_hat_new
    if(sum(fitidx)==0) break
  }

  # this block is taken from DESeq2, used for banning any dispersion estimates whose new estimate does not improve the
  # fitted log likelihood to increase by more than 1-millionth
  dispGeneEst <- alpha_hat
  if (niter == 1) {
    noIncrease <- dispRes$last_lp < dispRes$initial_lp + abs(dispRes$initial_lp)/1e6
    dispGeneEst[which(noIncrease)] <- alpha_init[which(noIncrease)]
  }

  # check the convergence status for each dispersion estimate
  dispGeneEstConv <- dispIter < maxit & !(dispIter == 1)
  refitDisp <- !dispGeneEstConv & dispGeneEst > minDisp*10 # refit these genes' dispersion parameter
  if (sum(refitDisp)>0){
    minLogAlpha = log(1e-8)
    maxLogAlpha = log(max(10, ncol(countsdata)))
    thegrid <- seq(from=minLogAlpha, to=maxLogAlpha, length=20)

    dispGrid = fitDispGrid(countsdata[refitDisp,], modelMatrix, mu[refitDisp,], disp_grid = thegrid, log_alpha_mean = rep(0,sum(refitDisp)),
                           log_alpha_s2 = 1, use_prior = F, useCR = useCR)
    dispGeneEst[refitDisp] = exp(dispGrid$log_alpha)
  }

  dispGeneEst <- pmin(pmax(dispGeneEst, minDisp), maxDisp)
  object$dispGeneEst = dispGeneEst
  object$dispGeneIter = dispIter
  object$mu = mu
  return(object)
}


# Function estimateDispersionFit: fit a trend line between the gene-wise dispersion parameter and the mean of normalized counts
#' @param object: the deseq_data object
#' @param minDisp: the minimum dispersion estimate allowed
#' @param fitType: categorical, the method used for fitting the trend function
#' @return object: contains the fitted trend function and each gene's fitted prior mean

estimateDispersionFit = function(object, minDisp=1e-8, fitType=c("parametric","local")){
  fitType = match.arg(fitType,choices=c("parametric","local"))
  if (is.null(object$allZero)) {
    object = getBaseMeansAndVariances(object)
  }

  # the following two lines are copied from DESeq2 source code, used for checking the plausibility of curve fitting
  useForFit = object$dispGeneEst > 100*minDisp
  if (sum(useForFit) == 0) {
    stop("all gene-wise dispersion estimates are within 2 orders of magnitude
  from the minimum value, and so the standard curve fitting techniques will not work.
  One can instead use the gene-wise estimates as final estimates:
  dds <- estimateDispersionsGeneEst(dds)
  dispersions(dds) <- mcols(dds)$dispGeneEst
  ...then continue with testing using nbinomWaldTest or nbinomLRT")
  }

  if (fitType == "parametric") {
    trial <- try(dispFunction <- parametricDispersionFit(object$baseMean[useForFit],
                                                         object$dispGeneEst[useForFit]),
                 silent=TRUE)
    if (inherits(trial,"try-error")) {
      message("-- note: fitType='parametric', but the dispersion trend was not well captured by the
   function: y = a/x + b, and a local regression fit was automatically substituted.
   specify fitType='local' or 'mean' to avoid this message next time.")
      fitType <- "local"
    }
  }
  if (fitType == "local") {
    dispFunction <- localDispersionFit(means = object$baseMean[useForFit],
                                       disps = object$dispGeneEst[useForFit],
                                       minDisp = minDisp)
  }
  if (!(fitType %in% c("parametric","local"))) {
    stop("unknown fitType")
  }else{
    # estimate the standard deviation of residuals of dispersion fitting
    dispFit = dispFunction(object$baseMean)
    resid = log(object$dispGeneEst) - log(dispFit)
    varLogDispEsts = mad(resid)
    attr(dispFunction,"varLogDispEsts") = varLogDispEsts^2
    object[["dispFunction"]] = dispFunction
    object[["dispFit"]] = dispFit
  }
  # the function will return a fitted trend function as equation 6 in DESeq2 paper
  return(object)
}

# Function estimateDispersionMAP: maximum a posteriori (MAP) estimate of the dispersion parameter
#' @param object: the deseq_data object
#' @param outlierSD: the number of standard deviations used for determining if a gene's gene-wise dispersion estimate is an outlier of the trend function
#' @param dispPriorVar: the prior distribution's variance of all gene's log dispersion
#' @param minDisp: the minimum dispersion allowed
#' @param kappa_0: the initial step size of backtracking search
#' @param dispTol: the tolerance level of convergence of the backtracking search algorithm
#' @param maxit: the maximum iteration of backtracking search for each gene's dispersion estimation
#' @param modelMatrix: the design matrix

estimateDispersionMAP = function(object, outlierSD=2, dispPriorVar, minDisp=1e-8, kappa_0=1,
                                 dispTol=1e-6, maxit=100, modelMatrix=NULL){
  countsdata = counts(object)
  if (is.null(object$allZero)) {
    object <- getBaseMeansAndVariances(object)
  }

  if (is.null(modelMatrix)){
    modelMatrix = object$design
  }

  if(missing(dispPriorVar)){
    # we are going to estimate the prior variance of dispersion parameter
    dispPriorVar <- estimateDispersionsPriorVar(object, modelMatrix=modelMatrix)
  }
  stopifnot(length(dispPriorVar)==1)

  # set prior variance for fitting dispersion
  log_alpha_prior_sigmasq <- dispPriorVar

  # set the NB mean value as the mean value fitted from GLM
  mu = object$mu

  # before using backtrack search for MAP, we need to specify the initial value of dispersion parameter
  dispInit = ifelse(object$dispGeneEst>0.1*object$dispFit,
                    object$dispGeneEst,object$dispFit)
  dispInit[is.na(dispInit)] = object$dispFit[is.na(dispInit)]

  dispResMAP = fitDisp(countsdata,modelMatrix,mu_hat = mu,log_alpha = log(dispInit),log_alpha_mean = log(object$dispFit),
                       log_alpha_s2 = log_alpha_prior_sigmasq, min_log_alpha = log(minDisp/10),kappa_0 = kappa_0,
                       tol = dispTol,maxit = maxit,use_prior = T,useCR =T)

  dispMAP = exp(dispResMAP$log_alpha) # the MAP estimate of the dispersion
  dispConv = dispResMAP$iter < maxit
  refitMAP = !dispConv

  if (sum(refitMAP)>0){
    minLogAlpha = log(1e-8)
    maxLogAlpha = log(max(10, ncol(countsdata)))
    thegrid <- seq(from=minLogAlpha, to=maxLogAlpha, length=20)

    dispGrid = fitDispGrid(countsdata[refitMAP,],x=modelMatrix,mu_hat = mu[refitMAP,],disp_grid = thegrid,
                           log_alpha_mean = log(object$dispFit)[refitMAP],log_alpha_s2 = log_alpha_prior_sigmasq,use_prior = T,useCR = T)
    dispMAP[refitMAP] = exp(dispGrid$log_alpha)
  }

  maxDisp = max(10,ncol(countsdata))
  dispMAP = pmin(pmax(dispMAP,minDisp),maxDisp)
  dispersionFinal = dispMAP

  # check for outliers
  varLogDispEsts = attr(object$dispFunction,"varLogDispEsts")
  dispOutlier = object$dispGeneEst > (object$dispFit + outlierSD * sqrt(varLogDispEsts))
  dispOutlier[is.na(dispOutlier)] = FALSE
  dispersionFinal[dispOutlier] = object$dispGeneEst[dispOutlier]

  resultsList <- list(dispersion = dispersionFinal,
                      dispIter = dispResMAP$iter,
                      dispOutlier = dispOutlier,
                      dispMAP = dispMAP)

  object[['MAPest']] = resultsList
  object[['dispersions']] = dispersionFinal
  return(object)
}


# Function nbinomWaldTest: estimate the log fold change betas, and do wald test to give conclusions on differential gene expression
#' @param object: the deseq_data object
#' @param modelMatrix: the design matrix
#' @param betaTol: the convergence threshold of finding the optimal betas for each gene in the IRLS
#' @param maxit: the maximum iteration of backtracking search for each gene's dispersion estimation
#' @param useT: use Student's t test or just normal distribution for Wald testing each gene's beta
#' @param useQR: use QR decomposition or not for finding the betas for each gene
#' @param minmu: the minimum fitted mean allowed
#' @return object: contains the wald test statistics and p-values


nbinomWaldTest = function(object, modelMatrix=NULL, betaTol=1e-8,
                          maxit=100, useT=TRUE, useQR=TRUE, minmu=0.5){
  if (is.null(object$allZero)){
    object = getBaseMeansAndVariances(object)
  }
  if (is.null(modelMatrix)){
    modelMatrix = getmodelMatrix(object)
  }
  countsdata = counts(object)

  fit <- fitNBinomGLM(object,x=modelMatrix,
                       betaTol=betaTol, maxit=maxit,
                        useQR=useQR,
                       minmu=minmu)
  H <- fit$hat_diagonal
  mu <- fit$mu
  modelMatrix <- fit$modelMatrix
  modelMatrixNames <- fit$modelMatrixNames
  # record the wide prior variance which was used in fitting
  betaPriorVar <- rep(1e6, ncol(fit$modelMatrix))



  #priorFitList = fitGLMsWithPrior(object, betaTol = betaTol, maxit = maxit, useQR = useQR, modelMatrix = modelMatrix, minmu = minmu)
  #fit <- priorFitList$fit
  #H <- priorFitList$H
  #mu <- priorFitList$mu
  #betaPriorVar <- priorFitList$betaPriorVar
  #modelMatrix <- priorFitList$modelMatrix
  #mleBetaMatrix <- priorFitList$mleBetaMatrix

  object[['mu']] = mu
  object[['H']] = H
  object[['betaPriorVar']] = betaPriorVar

  # now we do the wald test
  modelMatrixNames = colnames(modelMatrix)
  betaMatrix = fit$betaMatrix
  colnames(betaMatrix) = modelMatrixNames
  betaSE = fit$betaSE
  colnames(betaSE) = paste0("SE_",modelMatrixNames)
  WaldStatistic = betaMatrix/betaSE
  colnames(WaldStatistic) = paste0("WaldStatistic_",modelMatrixNames)

  # calculate the p-value for the wald statistics
  if(useT){
    num.samps <- rep(ncol(countsdata), nrow(countsdata))
    degree_f = num.samps - ncol(modelMatrix)
    stopifnot(length(degree_f)==nrow(WaldStatistic))
    WaldPvalue <- 2*pt(abs(WaldStatistic),df=degree_f,lower.tail=FALSE)
  }else{
    WaldPvalue <- 2*pnorm(abs(WaldStatistic),lower.tail=FALSE)
  }
  colnames(WaldPvalue) <- paste0("WaldPvalue_",modelMatrixNames)

  # check convergence
  betaConv <- fit$betaConv
  if (any(!betaConv)) {
   message(paste(sum(!betaConv),"rows did not converge in beta, labelled in object$betaConv. Use larger maxit argument with nbinomWaldTest"))
  }

  object[['Wald_Statistics']] = WaldStatistic
  object[['Wald_Pvalue']] = WaldPvalue
  return(object)
}


# Function DEA_deseq: the main function integrating all the functions above for running the NB-GLM model of differential expression analysis
#' @param data: gene's counts data, n by m
#' @param modelMatrix: sample design matrix, m by k
#' @return object: a deseq_data object, which is essentially a list, containing all the intermediate steps's estimations and the resulting wald test
#' p-values

DEA_deseq = function(data,modelMatrix){
  message("Setting up data object")
  object = deseq_data(data =data, designMatrix = modelMatrix)
  message("estimate size factors")
  nf = estimateSizeFactors(object)
  object$nf = nf
  message("estimate gene-wise dispersion parameters")
  object = estimateDispersionGeneEst(object,linearMu = T)
  message("trend function fitting")
  object = estimateDispersionFit(object)
  message("maximum a posteriori estimation of dispersion parameter")
  object = estimateDispersionMAP(object)
  message("wald test for differential expression")
  object = nbinomWaldTest(object,maxit = 100)
  message("analysis completed")
  return(object)
}









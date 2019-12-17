###################
# supplementary function of negative binomial GLM's LFC Wald Test

matchUpperQuantileForVariance <- function(x, upperQuantile=.05) {
  sdEst <- quantile(abs(x), 1 - upperQuantile) / qnorm(1 - upperQuantile/2)
  unname(sdEst)^2
}


estimateBetaPriorVar = function(object,modelMatrix=NULL){
  if(is.null(modelMatrix)){
    modelMatrix = getmodelMatrix(object)
  }
  
  betaMatrix = object$mleBetaMatrix
  betaPriorVar = apply(betaMatrix,2,function(x){
    useFinite = abs(x) < 10
    if(sum(useFinite)==0){
      return(1e6)
    }else{
      # here we use the quantile, unweighted method
      return(matchUpperQuantileForVariance(x))
    }
  })
  
  names(betaPriorVar) = colnames(betaMatrix)
  if ("Intercept" %in% names(betaPriorVar)) {
    betaPriorVar[which(names(betaPriorVar) == "Intercept")] <- 1e6
  }
  
  betaPriorVar
}


fitGLMsWithPrior = function(object, betaTol, maxit, useQR, modelMatrix=NULL, minmu=0.5){
  if(is.null(modelMatrix)){
    modelMatrix = getmodelMatrix(object)
  }
  
  # fit a naive GLM and get a rough estimate of LFC
  fit = fitNBinomGLM(object = object,x = modelMatrix, betaTol = betaTol, maxit = maxit, useQR = useQR, minmu = minmu)
  H = fit$hat_diagonal
  betaMatrix <- fit$betaMatrix
  mu <- fit$mu
  
  modelMatrixNames = colnames(modelMatrix)
  modelMatrixNames[modelMatrixNames == "(Intercept)"] = "Intercept"
  modelMatrixNames = make.names(modelMatrixNames)
  colnames(betaMatrix) = modelMatrixNames
  
  object[['mleBetaMatrix']] = betaMatrix # the estimates of LFC based on MLE of NG-GLM
  betaPriorVar = estimateBetaPriorVar(object,modelMatrix)
  
  if (any(betaPriorVar == 0)) {
    stop("beta prior variances are equal to zero for some variables")
  }
  
  lambda = 1/betaPriorVar
  fit = fitNBinomGLM(object,x=modelMatrix,lambda = lambda,betaTol = betaTol,maxit = maxit,useQR = useQR,minmu = minmu)
  
  res = list(fit=fit, H=H, betaPriorVar=betaPriorVar, mu=mu,
              modelMatrix=modelMatrix, mleBetaMatrix=betaMatrix)
  
  return(res)
}

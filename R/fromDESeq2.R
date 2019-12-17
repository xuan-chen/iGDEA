########### this file contains some functions copied from DESeq2 source code, which are some miscellaneous functions so we do not write our own version

library(locfit)

# Function linearModelMu: fit a linear model of y on x, used for giving an initial estimate of the mean of the negative binomial model of counts data
#' @param y: response variable
#' @param x: design matrix
#' @return the fitted mean, namely hat(y)

linearModelMu <- function(y, x) {
  qrx <- qr(x)    
  Q <- qr.Q(qrx)  
  Rinv <- solve(qr.R(qrx))
  # old code:
  # hatmatrix <- x %*% Rinv %*% t(Q)
  # t(hatmatrix %*% t(y))
  # Wolfgang Huber's rewrite is up to 2 orders of magnitude faster (Sept 2018):
  (y %*% Q) %*% t(x %*% Rinv)
}

# Function: roughtDispEstimate: give an estimation of the dispersion parameter, used for initializing dispersion parameter before backtracking search, 
# method similar to estimating the residual variance of a linear model
#' @param y: response variable
#' @param x: design matrix
#' @return the estimated dispersion parameter

roughDispEstimate <- function(y, x) {
  
  # must be positive
  mu <- linearModelMu(y, x)
  mu <- matrix(pmax(1, mu), ncol=ncol(mu))
  
  m <- nrow(x)
  p <- ncol(x)
  
  # an alternate rough estimator with higher mean squared or absolute error
  # (rowSums( (y - mu)^2/(mu * (m - p)) ) - 1)/rowMeans(mu)
  
  # rough disp ests will be adjusted up to minDisp later
  est <- rowSums( ((y - mu)^2 - mu) / mu^2 ) / (m - p)
  pmax(est, 0)
}

# Function: momentsDispEstimate: use method of moments to estimate dispersion parameter, used for giving an initial guess of the dispersion parameter
# before backtracking search, using method of moments
#' @param object: the deseq_data object
#' @return a rough estimate of dispersion parameter for all genes

momentsDispEstimate <- function(object) {
  xim <- mean(1/colMeans(getNormFactors(object,matformat = T)))
  if(is.null(object$baseVar)){
    object = getBaseMeansAndVariances(object)
  }
  bv <- object$baseVar
  bm <- object$baseMean
  return((bv - xim*bm)/bm^2)
}

# Function: parametericDipsersionFit: use gamma-family GLM to fit the gene-wise dispersion estimator on each gene's normalized counts
#' @param means: each gene's mean normalized counts
#' @param disps: gene-wise dispersion parameter
#' @return fitted trend function

parametricDispersionFit <- function(means,disps) {
  coefs <- c(.1, 1)
  iter <- 0
  while(TRUE) {
    residuals <- disps / ( coefs[1] + coefs[2] / means )
    good <- which( (residuals > 1e-4) & (residuals < 15) )
    # check for glm convergence below to exit while-loop
    suppressWarnings({fit <- glm( disps[good] ~ I(1/means[good]),
                                  family=Gamma(link="identity"), start=coefs )})
    oldcoefs <- coefs
    coefs <- coefficients(fit)
    if ( !all( coefs > 0 ) )
      stop(simpleError("parametric dispersion fit failed"))
    if ( ( sum( log( coefs / oldcoefs )^2 ) < 1e-6 )  & fit$converged )
      break
    iter <- iter + 1
    if ( iter > 10 ) 
      stop(simpleError("dispersion fit did not converge"))
  }
  names( coefs ) <- c( "asymptDisp", "extraPois" )
  ans <- function(q) coefs[1] + coefs[2] / q
  attr( ans, "coefficients" ) <- coefs
  ans

}

# Function: localDispersionFit: if parametric fitting failed, use a localized regression technique to find the trend function linking gene-wise
# dispersion parameter and each gene's normalized counts
#' @param means: each gene's mean normalized counts
#' @param disps: gene-wise dispersion parameter
#' @param minDisp: a lower bound for all gene-wise dispersion parameter, if a certain gene's dispersion is lower than minDisp*10, it is not used 
#' for fitting the localized regression
#' @return fitted trend function


localDispersionFit <- function(means,disps,minDisp ) {
  if (all(disps < minDisp*10)) {
    return(rep(minDisp,length(disps)))
  }
  d <- data.frame(logDisps = log(disps), logMeans = log(means))
  fit <- locfit(logDisps ~ logMeans, data=d[disps >= minDisp*10,,drop=FALSE],
                weights = means[disps >= minDisp*10])
  dispFunction = function(means) exp(predict(fit, data.frame(logMeans=log(means))))
  return(dispFunction)
}

# Function: linearModelMuNormalized:fit linear regression on the normalized counts data
#' @param object: the deseq_data object
#' @param x: the design matrix
#' @return the estimated mean of each count data

linearModelMuNormalized <- function(object, x) {
  cts <- counts(object)
  norm.cts <- counts(object, normalized=TRUE)
  muhat <- linearModelMu(norm.cts, x)
  nf <- getNormFactors(object)
  muhat * nf
}

# Function: fitNbnomGLMsOptim: use the built-in optim function to refit the negative binomial GLM for certain genes that does not converge under
# iterated re-weighted least squres (IRLS)
#' @param object: the deseq_data object
#' @param modelMatrix: the design matrix
#' @param lambda: the tuning parameter of ridge regression, which can also be used to specify the prior distribution of the betas of NB-GLM
#' @param rowsForOptim: genes selected to be put into the Optim procedure
#' @param rowStable: a boolean vector checking if a certain gene's beta estimator is stable
#' @param normalizationFactors: the normalization factor of each sample
#' @param alpha_hat: the dispersion parameter for each gene i, atomic vector of length n
#' @param betaMatrix: the previously fitted betas for each gene
#' @param betaSE: the previously fitted betas' variance
#' @param betaConv: the convergence status of each gene's betas from the fitNBinomGLM function
#' @param beta_mat: an initial guess of the beta, sometimes can replace using the fitted betas from the NB-GLM
#' @param mu the mean of the negative binomial distribution of each gene i's count data based on previously fitted NG-GLM, atomic vector of length m
#' @param logLike: loglikelihood of the counts data under the previously fitted NB-GLM
#' @param minmu: the minimum fitted mean allowed
#' @return list containing the betas of each gene, with some of the betas being replaced by the optim estimate in the function and the rest remaining
#' unchanged from the previous NB-GLM fitting.

fitNbinomGLMsOptim <- function(object,modelMatrix,lambda,
                               rowsForOptim,rowStable,
                               normalizationFactors,alpha_hat,
                               betaMatrix,betaSE,betaConv,
                               beta_mat,
                               mu,logLike,minmu=0.5) {
  x <- modelMatrix
  lambdaNatLogScale <- lambda / log(2)^2
  large <- 30
  for (row in rowsForOptim) {
    betaRow <- if (rowStable[row] & all(abs(betaMatrix[row,]) < large)) {
      betaMatrix[row,]
    } else {
      beta_mat[row,]
    }
    nf <- normalizationFactors[row,]
    k <- counts(object)[row,]
    alpha <- alpha_hat[row]
    objectiveFn <- function(p) {
      mu_row <- as.numeric(nf * 2^(x %*% p))
      logLikeVector <- dnbinom(k,mu=mu_row,size=1/alpha,log=TRUE)
      logLike = sum(logLikeVector)
      logPrior <- sum(dnorm(p,0,sqrt(1/lambda),log=TRUE))
      negLogPost <- -1 * (logLike + logPrior)
      if (is.finite(negLogPost)) negLogPost else 10^300
    }
    o <- optim(betaRow, objectiveFn, method="L-BFGS-B",lower=-large, upper=large)
    ridge <- if (length(lambdaNatLogScale) > 1) {
      diag(lambdaNatLogScale)
    } else {
      as.matrix(lambdaNatLogScale,ncol=1)
    }
    # if we converged, change betaConv to TRUE
    if (o$convergence == 0) {
      betaConv[row] <- TRUE
    }
    # with or without convergence, store the estimate from optim
    betaMatrix[row,] <- o$par
    # calculate the standard errors
    mu_row <- as.numeric(nf * 2^(x %*% o$par))
    # store the new mu vector
    mu[row,] <- mu_row
    mu_row[mu_row < minmu] <- minmu
    w <- diag(weights[row,] * (mu_row^-1 + alpha)^-1)
    xtwx <- t(x) %*% w %*% x
    xtwxRidgeInv <- solve(xtwx + ridge)
    sigma <- xtwxRidgeInv %*% xtwx %*% xtwxRidgeInv
    # warn below regarding these rows with negative variance
    betaSE[row,] <- log2(exp(1)) * sqrt(pmax(diag(sigma),0))
    logLikeVector <- dnbinom(k,mu=mu_row,size=1/alpha,log=TRUE)
    logLike[row] <- sum(logLikeVector)
  }
  return(list(betaMatrix=betaMatrix,betaSE=betaSE,
              betaConv=betaConv,mu=mu,logLike=logLike))
}




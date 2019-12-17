##############
# supplementary functions for estimating the dispersion parameter
# these functions come from the original Rcpp code written in C++ in DESeq2 package
# here we reimplement these functions in R

library(psych)
#source("deseq_obj.R")
#source("miscellaneous.R")

# Function log_posterior: calculate log posterior likelihood of counts data of gene i
#' @param y the count data for gene i, atomic vector of length m
#' @param x the design matrix, m by k (m: number of samples, k: covariates)
#' @param mu the mean of the negative binomial distribution of each gene i's count data, atomic vector of length m
#' @param log_alpha the natural log of dispersion parameter for gene i
#' @param log_alpha_mean the prior mean of the natural log of dispersion parameter for gene i
#' @param log_alpha_s2 the prior variance of the natural log of dispersion parameter for gene i
#' @param use_prior boolean, whether use the prior or not, if not, the function returns a pure log-likelihood of the count data of gene i
#' @param useCR boolean, whether to use Cox-Reid adjustment or not

log_posterior = function(y,x,mu,log_alpha,log_alpha_mean,log_alpha_s2,use_prior=T,useCR=T){
  alpha = exp(log_alpha)
  if (useCR){
    w_diag = diag(1/(alpha+1/mu)) # the weight matrix in cox-reid adjustment
    b = t(x)%*%w_diag%*%x
    cr_term = -0.5 * log(det(b)) # the cox-reid adjustment term
  }else{
    cr_term = 0.0
  }
  
  alpha_neg1 = 1/alpha # 1/alpha
  
  # log-likelihood for all gene i's samples, under its negative binomial distribution specified by mu and alpha
  ll_part = sum(lgamma(y+alpha_neg1) - lgamma(alpha_neg1) - y*log(mu+alpha_neg1) - alpha_neg1 * log(1.0+mu*alpha))
  
  if (use_prior){
    # the prior is a log-normal distribution
    prior_part = -0.5 * (log_alpha-log_alpha_mean)^2/log_alpha_s2
  }else{
    prior_part = 0.0
  }
  
  res = ll_part + prior_part + cr_term
  return(res)
}

# Function dlog_posterior: calculate the derivative of the posterior density of count data of gene i w.r.t log(alpha_i)
#' @param y the counts data for gene i, atomic vector of length m
#' @param x the design matrix, m by k (m: number of samples, k: covariates)
#' @param mu the mean of the negative binomial distribution of each gene i's count data, atomic vector of length m
#' @param log_alpha the natural log of dispersion parameter for gene i
#' @param log_alpha_mean the prior mean of the natural log of dispersion parameter for gene i
#' @param log_alpha_s2 the prior variance of the natural log of dispersion parameter for gene i
#' @param use_prior boolean, whether use the prior or not, if not, the function returns a pure log-likelihood of the count data of gene i
#' @param useCR boolean, whether to use Cox-Reid adjustment or not

dlog_posterior = function(y,x,mu,log_alpha,log_alpha_mean,log_alpha_s2,use_prior=T,useCR=T){
  alpha = exp(log_alpha)

  if (useCR){
    w_diag = diag(1/(alpha+1/mu))
    dw_diag = diag(-1.0*(alpha+1/mu)^(-2))
    b = t(x)%*%w_diag%*%x
    db = t(x)%*%dw_diag%*%x
    ddetb = (det(b)*tr(solve(b)%*%db))
    cr_term = -0.5 * ddetb / det(b)
  }else{
    cr_term = 0.0
  }
  
  alpha_neg1 = 1/alpha
  alpha_neg2 = 1/(alpha)^2
  ll_part = alpha_neg2 * sum(digamma(alpha_neg1) + log(1+mu*alpha) - (mu*alpha)/(1.0+mu*alpha) - 
    digamma(y + alpha_neg1) + y/(mu + alpha_neg1))
  
  if (use_prior){
    prior_part = -1.0 * (log_alpha - log_alpha_mean)/log_alpha_s2
  }else{
    prior_part = 0.0
  }
  res = (ll_part + cr_term) * alpha + prior_part
  return(res)
}


# Function d2log_posterior: calculate the second order derivative of the posterior density of count data of gene i w.r.t log(alpha_i)
#' @param y the counts data for gene i, atomic vector of length m
#' @param x the design matrix, m by k (m: number of samples, k: covariates)
#' @param mu the mean of the negative binomial distribution of each gene i's count data, atomic vector of length m
#' @param log_alpha the natural log of dispersion parameter for gene i
#' @param log_alpha_mean the prior mean of the natural log of dispersion parameter for gene i
#' @param log_alpha_s2 the prior variance of the natural log of dispersion parameter for gene i
#' @param use_prior boolean, whether use the prior or not, if not, the function returns a pure log-likelihood of the count data of gene i
#' @param useCR boolean, whether to use Cox-Reid adjustment or not

d2log_posterior = function(y,x,mu,log_alpha,log_alpha_mean,log_alpha_s2,use_prior=T,useCR=T){
  alpha = exp(log_alpha)
  x_orig = x
  
  if (useCR){
    w_diag = diag(1/(alpha + 1/mu))
    dw_diag  = diag(-1*(alpha + 1/mu)^(-2))
    d2w_diag = diag(2*(alpha + 1/mu)^(-3))
    b = t(x) %*% w_diag %*% x
    b_i = solve(b)
    db = t(x) %*% dw_diag %*% x
    d2b = t(x) %*% d2w_diag %*% x
    ddetb = det(b)*tr(b_i%*%db)
    d2detb = det(b)*(tr(b_i%*%db)^2 - tr(b_i%*%db%*%b_i%*%db) + tr(b_i%*%d2b))
    cr_term = 0.5 * (ddetb/det(b))^2 -0.5*d2detb/det(b)
  }else{
    cr_term = 0.0
  }
  alpha_neg1 = 1/alpha
  alpha_neg2 = alpha^(-2)
  ll_part = -2 * alpha^(-3) *sum(digamma(alpha_neg1) + log(1 + mu*alpha) - mu*alpha/(1+mu*alpha) -
                                   digamma(y + alpha_neg1) + y/(mu+alpha_neg1)) + alpha_neg2 * sum(-1 * alpha_neg2 * trigamma(alpha_neg1)+
                                  mu^2*alpha*(1+mu*alpha)^(-2) + alpha_neg2 * trigamma(y + alpha_neg1) + alpha_neg2 * y * (mu+alpha_neg1)^(-2))
  
  if (use_prior){
    prior_part = -1.0/log_alpha_s2; 
  }else{
    prior_part = 0.0
  }
  
  res = ((ll_part + cr_term) * alpha^2 + dlog_posterior(y,x,mu,log_alpha,log_alpha_mean,log_alpha_s2,use_prior,useCR)) + prior_part
  return(res)
}



# Function fitDisp: fit each gene's dispersion parameter based on MLE of likelihood of each gene's counts data, with or without Cox-Reid Adjustment
#' @param y counts data, n by m, n: number of genes, m: number of samples
#' @param x design matrix, m by k, m: number of samples, k: number of covariates
#' @param mu_hat negative-binomial(NB) mean matrix: n by m, the fitted mean for each gene i, sample j's count
#' @param log_alpha initial gusess for each gene's dispersion, atomic vector of length n
#' @param log_alpha_mean prior mean of each gene's dispersion, atomic vector of length n
#' @param log_alpha_s2 prior variance of each gene's dispersion, a single value, shared by all genes, see DESeq2 paper's p.14 equation (5)
#' @param min_log_alpha minimum log dispersion allowed, set for lower-bounding the dispersion estimate for the sake of numeric stability
#' @param kappa_0 parameter used for backtracking search for finding the MLE dispersion
#' @param tol tolerance for convergence of estimates
#' @param maxit maximum iteration for estimating each gene's disperision
#' @param use_prior boolean, whether prior of dispersion shall be used in calculating the posterior
#' @param useCR boolean, whether Cox-Reid bias adjustment is used for MLE estimation
#' @return list of all important quantities about the dispersion estimate, including the estimate, convergence status, iteration steps and more


fitDisp = function(y, x, mu_hat, log_alpha, log_alpha_mean, log_alpha_s2, min_log_alpha, 
                   kappa_0, tol, maxit, use_prior, useCR){
  
  y_n = nrow(y) # number of genes
  epsilon = 1.0e-4
  
  # book keeping
  initial_lp <- initial_dlp <- last_lp <- last_dlp <- last_d2lp <- last_change <- iter <- iter_accept <- rep(0,y_n)
  
  for (i in 1:y_n){
    yrow = y[i,] # atomic vector of length m
    mu_hat_row = mu_hat[i,] # atomic vector of length m
    a = log_alpha[i] # the initial guess of dispersion
    
    
    lp = log_posterior(yrow,x,mu_hat_row,a,log_alpha_mean[i],log_alpha_s2,use_prior,useCR)
    dlp = dlog_posterior(yrow,x,mu_hat_row,a,log_alpha_mean[i],log_alpha_s2,use_prior,useCR)
    kappa1 = kappa_0
    initial_lp[i] = lp
    initial_dlp[i] = dlp
    change = -1.0
    last_change[i] = -1.0
    
    # use backtracking search to find the MLE for the log_posterior
    for (t in 1:maxit){
      iter[i] = iter[i] + 1
      a_propose = a + kappa1 * dlp # propose log_alpha to be dlp from the current log_alpha
      
      # limit the proposed log_alpha to be no smaller than -30 and no greater than 10 for numeric stability
      if (a_propose < -30.0) {
        kappa1 = (-30.0 - a)/dlp
      }
      
      if (a_propose > 10.0) {
        kappa1 = (10.0 - a)/dlp
      }
      
      theta_kappa1 = -1.0*log_posterior(yrow,x,mu_hat_row,a+kappa1*dlp,log_alpha_mean[i],log_alpha_s2,use_prior,useCR)
      theta_hat_kappa1 = -1.0*lp - kappa1*epsilon*(dlp^2)
      
      # each accepted proposed alpha should satisfy Armijo condition
      if (theta_kappa1<=theta_hat_kappa1){
        iter_accept[i] = iter_accept[i] + 1
        a = a + kappa1 * dlp
        lpnew = log_posterior(yrow,x,mu_hat_row,a,log_alpha_mean[i],log_alpha_s2,use_prior,useCR)
        change = lpnew - lp
        
        # check convergence
        if (change < tol) {
          lp = lpnew
          break
        }

        # check if the accepted log_alpha is too small
        if (a<min_log_alpha){
          break
        }
        
        lp = lpnew
        dlp = dlog_posterior(yrow,x,mu_hat_row,a,log_alpha_mean[i],log_alpha_s2,use_prior,useCR)
        
        # after the accept, we enlarge our step size along the derivative's direction by 10%
        kappa1 = min(kappa1 * 1.1, kappa_0)
        
        # for every 5 accepts, shrink the kappa1 into half
        if (iter_accept[i] %% 5 == 0) {
          kappa1 = kappa1 / 2.0
        }
      }else{
        kappa1 = kappa1/2.0 # this is the key part of "backtracking"
      }
    }
    
    last_lp[i] = lp
    last_dlp[i] = dlp
    last_d2lp[i] = d2log_posterior(yrow,x,mu_hat_row,a,log_alpha_mean[i],log_alpha_s2,use_prior,useCR)
    log_alpha[i] = a # the backtrach search MLE estimate of the dispersion parameter of gene i
    last_change[i] = change
    
  }
  return(list(log_alpha = log_alpha, iter = iter, iter_accept = iter_accept, last_change = last_change, initial_lp = initial_lp,
              initial_dlp = initial_dlp, last_lp = last_lp, last_dlp = last_dlp, last_d2lp = last_d2lp))
}


# Function fitBeta: fit negative binomial GLM for counts data with dispersion for each gene i given
#' @param y counts data, n by m, n: number of genes, m: number of samples
#' @param x design matrix, m by k, m: number of samples, k: number of covariates
#' @param nf normalization factor matrix, n by m
#' @param alpha_hat dispersion estimates for each gene i, atomic vector of length n
#' @param beta_mat initial guess for the beta parameter matrix, n by k
#' @param lambda L2-penalized GLM regression tuning parameter, atomic vector of length k
#' @param tol tolerance for convergence in estimates
#' @param maxit maximum number of IRLS iteration
#' @param useQR boolean, whether use QR decomposition for solving the NB-GLM
#' @param minmu minimum of NB mean value allowed
#' @return list of all important quantities about the beta estimators, including the beta point estimator, variance of the beta, projection matrix
#' and deviance of the fitted GLM model


fitBeta = function(y, x, nf, alpha_hat, beta_mat, lambda, tol, maxit, useQR=F,minmu){
  y_n = nrow(y) # number of genes
  y_m = ncol(y) # number of samples
  x_p = ncol(x) # number of covariates

  beta_var_mat = matrix(0,nrow(beta_mat),ncol(beta_mat))
  hat_diagonals = matrix(0,y_n,y_m)
  
  large = 30.0
  iter <- deviance <- rep(0.0,y_n) # each gene's IRLS iteration times and final deviance of the fitted GLM model
  
  for (i in 1:y_n){
    nfrow = matrix(nf[i,],y_m,1)
    yrow = matrix(y[i,],y_m,1)
    beta_hat = matrix(beta_mat[i,],x_p,1)
    mu_hat = nfrow * exp(x%*%beta_hat)
    # lower bounding the mean value
    for (j in 1:y_m){
      mu_hat[j] = max(mu_hat[j], minmu)
    }
    
    ridge = diag(lambda) # ridge regression matrix lambda*I
    dev <- dev_old <- 0.0 # deviance, used for convergence check after each iteration of IRLS
    
    # start the iteration of IRLS
    if (useQR){
      for (t in 1:maxit){
        iter[i] = iter[i] + 1
        w_vec = diag(as.numeric(mu_hat/(1.0 + alpha_hat[i] * mu_hat)))
        w_sqrt_vec = sqrt(w_vec) # weight matrix of IRLS
        #weighted_x_ridge = rbind(w_sqrt_vec%*%x, sqrt(ridge))
        weighted_x_ridge = t(x)%*%w_vec%*%x + ridge
        qrobj = qr(weighted_x_ridge)
        q = qr.Q(qrobj)
        r = qr.R(qrobj)
        #big_w_diag = matrix(1,y_m+x_p,1)
        #big_w_diag[1:y_m,] = diag(w_vec)
        z = log(mu_hat/nfrow) + (yrow-mu_hat)/mu_hat
        #z_sqrt_w = w_sqrt_vec%*%z
        #big_z_sqrt_w = matrix(0, y_m+x_p, 1)
        #big_z_sqrt_w[1:y_m,] = diag(z_sqrt_w)
        gamma_hat = t(q)%*%t(x)%*%w_vec%*%z
        beta_hat = matrix(solve(r,gamma_hat),x_p,1) # the beta fitted for gene i after 1 IRLS iteration
        
        # any beta too large, break the loop and move on to the next gene
        if (sum(abs(beta_hat) > large) > 0){
          iter[i] = maxit
          break
        }
        
        # update the fitted mean value, lower-bounding it by minmu
        mu_hat = nfrow * exp(x%*%beta_hat)
        for (j in 1:y_m){
          mu_hat[j] = max(mu_hat[j], minmu)
        }
        
        # check convergence condition using deviance
        dev = 0.0
        for (j in 1:y_m){
          dev = dev + -2.0 * dnbinom(x = yrow[j], size = 1.0/alpha_hat[i], mu = mu_hat[j], log = 1) # use the alternative parameterization
        }
        conv_test = abs(dev - dev_old)/(abs(dev) + 0.1)
        
        # if the deviance test statistics cannot be calculated, we do not iterate anymore
        if(is.nan(conv_test)){
          iter[i] = maxit
          break
        }
        
        # we need at least 2 iterations, and we stop iterating more if the fitted model's deviance does not change anymore
        if ((t>1) & (conv_test<tol)){
          break
        }
        dev_old = dev
      }
    }else{
      # if one does not use QR for solving the IRLS, we solve for the NB-GLM solution directly
      for (t in 1:maxit){
        iter[i] = iter[i]+1
        w_vec = diag(as.numeric(mu_hat/(1.0 + alpha_hat[i] * mu_hat)))
        w_sqrt_vec = sqrt(w_vec) # weight matrix of IRLS
        z = matrix(log(mu_hat/nfrow) + (yrow-mu_hat)/mu_hat,y_m,1)
        beta_hat = solve(t(x)%*%w_vec%*%x+ridge,t(x)%*%(z*diag(w_vec))) # IRLS iteration step
      }
      
      # any beta too large, break the loop and move on to the next gene
      if (sum(abs(beta_hat) > large) > 0){
        iter[i] = maxit
        break
      }
      
      # update the fitted mean value, lower-bounding it by minmu
      mu_hat = nfrow * exp(x%*%beta_hat)
      for (j in 1:y_m){
        mu_hat[j] = max(mu_hat[j], minmu)
      }
      
      # check convergence condition using deviance
      dev = 0.0
      for (j in 1:y_m){
        dev = dev + -2.0 * dnbinom(x = yrow[j], size = 1.0/alpha_hat[i], mu = mu_hat[j], log = 1) # use the alternative parameterization
      }
      conv_test = abs(dev - dev_old)/(abs(dev) + 0.1)
      
      # if the deviance test statistics cannot be calculated, we do not iterate anymore
      if(is.nan(conv_test)){
        iter[i] = maxit
        break
      }
      
      # we need at least 2 iterations, and we stop iterating more if the fitted model's deviance does not change anymore
      if ((t>1) & (conv_test<tol)){
        break
      }
      dev_old = dev
      
    }
    
    deviance[i] = dev # the fitted model's deviance
    beta_mat[i,] = t(beta_hat) # the fitted beta_mat
    
    #
    w_vec = diag(as.numeric(mu_hat/(1.0 + alpha_hat[i] * mu_hat)))
    w_sqrt_vec = sqrt(w_vec) # weight matrix of IRLS
    
    
    hat_matrix_diag = matrix(0,nrow(x),1)
    xw = w_sqrt_vec%*%x
    xtwxr_inv = solve(t(x)%*%w_vec%*%x+ridge)
    for (jp in 1:y_m){
      for (idx1 in 1:x_p){
        for (idx2 in 1:x_p){
          hat_matrix_diag[jp] = hat_matrix_diag[jp] + xw[jp,idx1] * (xw[jp,idx2]*xtwxr_inv[idx2,idx1])
        }
      }
    }
    hat_diagonals[i,] = t(hat_matrix_diag)
    
    # sigma is the covariance matrix for the betas, using sandwiched expression
    sigma = solve((t(x)%*%w_vec%*%x+ridge))%*%t(x)%*%(w_vec%*%x)%*%solve((t(x)%*%w_vec%*%x+ridge))
    beta_var_mat[i,] = diag(sigma)
    
  }
  
  return(list(beta_mat = beta_mat, beta_var_mat = beta_var_mat, iter = iter, hat_diagonals = hat_diagonals,
              deviance = deviance))
}

# Function fitDispGrid: when the backtracking search method does not work well for certain gene's dispersion parameter estimation,
# a grid search is applied instead. 

#' @param y counts data, n by m, n: number of genes, m: number of samples
#' @param x design matrix, m by k, m: number of samples, k: number of covariates
#' @param mu_hat negative-binomial(NB) mean matrix: n by m, the fitted mean for each gene i, sample j's count
#' @param disp_grid the grid over which to estimate, atomic vector of length N0
#' @param log_alpha_mean prior mean of each gene's dispersion, atomic vector of length n
#' @param log_alpha_s2 prior variance of each gene's dispersion, a single value, shared by all genes, see DESeq2 paper's p.14 equation (5)
#' @param use_prior boolean, whether prior of dispersion shall be used in calculating the posterior
#' @param useCR boolean, whether Cox-Reid bias adjustment is used for MLE estimation
#' @return list containing the estimated dispersion parameter

fitDispGrid = function(y,x,mu_hat,disp_grid,log_alpha_mean,log_alpha_s2,use_prior=T,useCR=T){
  if(is.vector(y)){
    y = matrix(y,1)
  }
  if(is.vector(mu_hat)){
    mu_hat = matrix(mu_hat,1)
  }
  y_n = nrow(y)
  grid_n = length(disp_grid)
  delta = disp_grid[2] - disp_grid[1]
  
  logpostvec = rep(0.0,grid_n)
  log_alpha = rep(0.0,y_n)
  
  for (i in 1:y_n){
    yrow = y[i,]
    mu_hat_row = mu_hat[i,]
    
    for (t in 1:grid_n){
      a = disp_grid[t]
      logpostvec[t] = log_posterior(yrow,x,mu_hat_row,a,log_alpha_mean[i],log_alpha_s2,use_prior,useCR)
    }
    idx = which(logpostvec==max(logpostvec)) # find the index of the alpha that maximize the log-posterior
    a_hat = disp_grid[idx] # the best log_alpha candidate
    
    # now with the candidate a_hat, we will search around a +/- delta neighborhood of the candidate a_hat to search for 
    # an even better log_alpha
    disp_grid_fine = seq(from = a_hat-delta, to = a_hat+delta, length.out = grid_n)
    
    for (t in 1:grid_n){
      a = disp_grid_fine[t]
      logpostvec[t] = log_posterior(yrow,x,mu_hat_row,a,log_alpha_mean[i],log_alpha_s2,use_prior,useCR)
    }
    idx = which(logpostvec==max(logpostvec)) # find the index of the alpha that maximize the log-posterior
    log_alpha[i] = disp_grid_fine[idx] # the best log_alpha candidate
  }
  return(list(log_alpha = log_alpha))
}


# Function fitNBinomGLM: fit negative binomial GLM
#' @param object a deseq_data object
#' @param x the design matrix, which should also be stored in object$design
#' @param alpha_hat the dispersion parameter for each gene i, atomic vector of length n
#' @param lambda the tuning parameter for ridge regression, penalizing the GLM betas
#' @param betaTol the tolerance level of the IRLS iteration (exactly the same as tol in fitBeta)
#' @param maxit maximum iteration for estimating NB-GLM beta
#' @param useQR boolean, whether using QR for NB-GLM IRLS
#' @param minmu the minimum fitted mean value allowed
#' @return list of all important quantities about the NB-GLM's betas, including beta point estimate, variance and iteration steps of IRLS and
#' convergence status

fitNBinomGLM = function(object,x,alpha_hat,lambda,betaTol=1e-8,maxit=100,useQR=T,minmu=0.5, useOptim=TRUE){
  if(missing(x)){
    x = getmodelMatrix(object)
  }
  
  stopifnot(all(colSums(abs(x)) > 0))
  modelMatrixNames <- colnames(x)
  modelMatrixNames[modelMatrixNames == "(Intercept)"] <- "Intercept"
  modelMatrixNames <- make.names(modelMatrixNames) # clean up the names
  colnames(x) <- modelMatrixNames
  nf <- getNormFactors(object,matformat=T) # get the normalization factors matrix, used for initialized betas
  countsdata = counts(object,normalized = F)
  
  if(missing(alpha_hat)){
    alpha_hat = dispersions(object)
  }
  
  if (length(alpha_hat) != nrow(counts(object))) {
    stop("alpha_hat needs to be the same length as number of genes!")
  }
  
  # set a wide prior for all coefficients
  if (missing(lambda)) {
    lambda <- rep(1e-6, ncol(x))
  }
  
  # initialize the betas (fold change) for each gene
  qrx <- qr(x)
  # if full rank, estimate initial betas for IRLS below
  if (qrx$rank == ncol(x)) {
    Q <- qr.Q(qrx)
    R <- qr.R(qrx)
    y <- t(log(counts(object,normalized=TRUE) + .1)) # m by n matrix
    beta_mat <- t(solve(R, t(Q) %*% y)) # use linear model's estimate as initial guess of NB-GLM estimate
  } else {
    if ("Intercept" %in% modelMatrixNames) {
      beta_mat <- matrix(0, ncol=ncol(x), nrow=nrow(counts(object)))
      # use the natural log as fitBeta occurs in the natural log scale
      logBaseMean <- log(rowMeans(counts(object,normalized=TRUE)))
      beta_mat[,which(modelMatrixNames == "Intercept")] <- logBaseMean
    } else {
      beta_mat <- matrix(1, ncol=ncol(x), nrow=nrow(counts(object)))
    }
  }
  
  
  lambdaNatLogScale <- lambda / log(2)^2
  betaRes <- fitBeta(countsdata, x, nf, alpha_hat, beta_mat, lambdaNatLogScale, betaTol, maxit, useQR, minmu)
  mu <- nf * t(exp(x %*% t(betaRes$beta_mat))) # the fitted mean
  betaMatrix <- log2(exp(1))*betaRes$beta_mat
  colnames(betaMatrix) <- modelMatrixNames
  logLike = nbinomLogLike(countsdata,mu,dispersions = alpha_hat) # calculate log-likelihood
  betaConv = betaRes$iter < maxit # check convergence status
  
  # the standard error for each beta estimator
  betaSE <- log2(exp(1))*sqrt(pmax(betaRes$beta_var_mat,0))
  colnames(betaSE) <- paste0("SE_",modelMatrixNames)
  stopifnot(!any(is.na(betaSE)))
  
  # check for stability
  rowStable = apply(betaRes$beta_mat,1,function(row) sum(is.na(row))) == 0
  
  # check for positive variances
  rowVarPositive = apply(betaRes$beta_var_mat,1,function(row) sum(row <= 0)) == 0
  
  # some of the GLM needs to be fit again
  rowsForOptim <- if (useOptim) {
    which(!betaConv | !rowStable | !rowVarPositive)
  } else {
    which(!rowStable | !rowVarPositive)
  }
  
  if (length(rowsForOptim) > 0) {
    # we use optim if didn't reach convergence with the IRLS code
    resOptim <- fitNbinomGLMsOptim(object,modelMatrix = x,lambda,
                                   rowsForOptim = rowsForOptim,rowStable = rowStable,
                                   normalizationFactors = nf,alpha_hat,
                                   betaMatrix = betaMatrix,betaSE = betaSE,betaConv = betaConv,
                                   beta_mat = beta_mat,
                                   mu = mu,logLike = logLike,minmu = minmu)
    betaMatrix <- resOptim$betaMatrix
    betaSE <- resOptim$betaSE
    betaConv <- resOptim$betaConv
    mu <- resOptim$mu
    logLike <- resOptim$logLike
  }
  
  return(list(logLike = logLike, betaConv = betaConv, betaMatrix = betaMatrix,
       betaSE = betaSE, mu = mu, betaIter = betaRes$iter, modelMatrix=x, 
       nterms=ncol(x), hat_diagonal=betaRes$hat_diagonals))
}


# Function estimateDispersionPriorVar: find the variance of the prior distribution of log-alpha, namely the sigma2_d
# see details on p.16 of DESeq2 paper, here we only consider cases where m>(p+3)
#' @param object: the deseq_data object
#' @param modelMatrix: the design matrix
#' @param minDisp: the minimum dispersion allowed
#' @return the estiumated prior variance of log of dispersion

estimateDispersionsPriorVar = function(object,modelMatrix=NULL,minDisp=1e-8){
  aboveMinDisp = object$dispGeneEst>=minDisp*100
  if(is.null(modelMatrix)){
    modelMatrix = object$design
  }
  
  dispResidual = log(object$dispGeneEst) - log(object$dispFit)
  
  
  if (sum(aboveMinDisp,na.rm=TRUE) == 0) {
    stop("no data found which is greater than minDisp")
  }
  
  varLogDispEsts = mad(dispResidual[aboveMinDisp])^2
  attr(object$dispFunction,"varLogDispEsts") = varLogDispEsts# this is the s_lr
  m = nrow(modelMatrix)
  p = ncol(modelMatrix)
  
  disp_sample_error = trigamma((m-p)/2)
  disp_priorvar = pmax(varLogDispEsts-disp_sample_error,0.25)
  return(disp_priorvar)
}





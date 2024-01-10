#############################################################
#	scMANOVA_H1 function
#	Authors: E. Sabbioni, C. Agostinelli and A. Farcomeni
#	Maintainer e-mail: elena.sabbioni@polito.it
#	Date: 05 October 2023
#	Version: 0.1
#	Copyright (C) 2023 E. Sabbioni, C. Agostinelli and A. Farcomeni
#############################################################

scMANOVA_H1 <- function(x, n, lambda=NULL, lambda0=NULL, lambda.step=0.1, ident=FALSE, tol=1e-8, penalty=function(n, p) log(n), rm.vars=NA) {
  N <- sum(n)
  K <- length(n)
  if (is.null(lambda)) {
    lambda.min <- 0
    lambda.max <- 100
  } else if (length(lambda)==2) {
    lambda.min <- lambda[1]
    lambda.max <- lambda[2]
  } else if (length(lambda)!=1) {
    stop("'lambda' must be 'NULL', a scalar or a vector of length 2")
  }
  if (is.null(lambda0))
    stop("'lambda0' must be a scalar value")

## Estimation
  res <- scMANOVAestimation(x, n, lambda=0, lambda0=lambda0, ident=ident, posdef.check=FALSE, rm.vars=rm.vars)
  mu <- res$mu
  mu0 <- res$mu0
  sigma <- res$sigma
  pp <-  nrow(sigma)
  sigma0 <- res$sigma0
  mu[is.nan(mu)] <- mu0[is.nan(mu)] 
  cor0 <- cov2cor(sigma0)
  dna <- is.na(diag(sigma))
  if (any(dna))
    diag(sigma)[dna] <- diag(sigma0)[dna] 
  ona <- which(is.na(sigma), arr.ind=TRUE)
  if (NROW(ona)) {
    for (i in 1:NROW(ona)) {
      sigma[ona[i,1],ona[i,2]] <- cor0[ona[i,1],ona[i,2]]*sqrt(sigma[ona[i,1],ona[i,1]]*sigma[ona[i,2],ona[i,2]])
    }
  }
  sigma <- (sigma+t(sigma))/2
  logLikPi <- res$logLikPi
  logLik0 <- res$logLik0  
  if (length(res$removed.vars))
    x <-  x[, -res$removed.vars, drop=FALSE]
  x[x==0] <- NA  
  y <- !is.na(x)
  gr <- rep(1:K, times=n)
  xlist <- list()
  for (k in 1:K) {
    pos <- which(gr==k)
    xlist[[k]] <- x[pos,,drop=FALSE]
  }

    #### H_1
    if (pp!=1) {
      if (length(lambda)!=1) {
        s <- if (ident) {
          diag(pp)
        } else {
          diag(diag(sigma), pp)
        }
        def_pos <- is.positive.definite(sigma + lambda.min * s)
        lambda.last <- lambda.min
        while (!def_pos) {
          lambda.last <- lambda.min
          lambda.min <- lambda.min + lambda.step
          def_pos <- is.positive.definite(sigma + lambda.min * s)
          lambda.step <- 2*lambda.step
        }
        d <- lambda.min - lambda.last
        while (d > tol) {
          temp <- (lambda.min+lambda.last)/2
          def_pos <- is.positive.definite(sigma + temp*s)
          if (def_pos) {
            lambda.min <- temp
          } else {
            lambda.last <- temp
          }
          d <- lambda.min - lambda.last
        }

        if (lambda.min > lambda.max)
          stop("'lambda[2]' too small to make the variance-covariance positive definite")

        f1 <- function(lambda) { 
          sigma_ridge <- sigma+lambda*s
          logLik <- logLikPi
          for (k in 1:K) {    
            logLik <- logLik+sum(my.dmvnorm(xlist[[k]], mu[k,], sigma_ridge))
          }
           -2*logLik+penalty(N, pp)*df(y=y, sigma=sigma_ridge)
        }
        lambda <- optimize(f1, lower=lambda.min, upper=lambda.max)$minimum
      }
      sigmaRidge <- sigma+lambda*s
    } else {
      sigmaRidge <- sigma
      if (sigma==0 | is.na(sigma))
        stop("only one variable remained with 0 or NA sample variance under H_1")
    }
  
  #####################################################
  # Likelihood
  #####################################################
    logLik <- logLikPi
    for (k in 1:K) {    
      logLik <- logLik + sum(my.dmvnorm(xlist[[k]], mu[k,], sigmaRidge))
    }
  #####################################################
  # Test statistic
  #####################################################  
    res$statistic <- -2 * (logLik0 - logLik)
    res$sigmaRidge <- sigmaRidge
    res$logLik <- logLik
    res$lambda <- lambda
    res$df <- df(y=y, sigma=sigmaRidge)
    res$aic <- -2*logLik + penalty(N, pp)*res$df
  return(res)
}

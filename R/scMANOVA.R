#############################################################
#	main function
#	Authors: E. Sabbioni, C. Agostinelli and A. Farcomeni
#	Maintainer e-mail: elena.sabbioni@polito.it
#	Date: 17 April 2023
#	Version: 0.2-2
#	Copyright (C) 2023 E. Sabbioni, C. Agostinelli and A. Farcomeni
#############################################################

scMANOVA <- function(x, n, lambda=NULL, lambda0=NULL, lambda.step=0.1, ident=FALSE, tol=1e-8, penalty=function(n, p) log(n), B=500, p.value.perm=FALSE, fixed.lambda=FALSE, ...) {
  N <- sum(n)
  K <- length(n)
  if (is.null(lambda)) {
    ll.min <- lambda.min <- 0
    lambda.max <- 100
  } else if (length(lambda)==2) {
    ll.min <- lambda.min <- lambda[1]
    lambda.max <- lambda[2]
  } else if (length(lambda)!=1) {
    stop("'lambda' must be 'NULL', a scalar or a vector of length 2")
  }
  if (is.null(lambda0)) {
    lambda0.min <- 0
    lambda0.max <- 100
  } else if (length(lambda0)==2) {
    lambda0.min <- lambda0[1]
    lambda0.max <- lambda0[2]
  } else if (length(lambda0)!=1) {
    stop("'lambda0' must be 'NULL', a scalar or a vector of length 2")
  }

## Estimation
##  if (p.value.perm)
##    xx <- x
  res <- scMANOVAestimation(x, n, lambda=ifelse(length(lambda)==1, lambda, 0), lambda0=ifelse(length(lambda0)==1, lambda0, 0), ident=ident, posdef.check=FALSE, rm.vars=NA)
  mu <- res$mu
  mu0 <- res$mu0
  sigma <- res$sigma
  sigma0 <- res$sigma0
  pp <-  nrow(sigma)
  lambda0.step <- lambda.step
  logLikPi <- res$logLikPi
  logLikPi0 <- res$logLikPi0
  rm.vars <- res$removed.vars

  if (p.value.perm) {
    xx <- x
  }
  
  if (length(rm.vars))
    x <- x[, -rm.vars, drop=FALSE]
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
  
  #### H_0
    if (length(lambda0)!=1) {
      s0 <- if (ident) {
        diag(pp)
      } else {
        diag(diag(sigma0), pp)
      }  
      def_pos <- is.positive.definite(sigma0+lambda0.min*s0)
      lambda0.last <- lambda0.min
      while (!def_pos) {
        lambda0.last <- lambda0.min
        lambda0.min <- lambda0.min + lambda0.step
        def_pos <- is.positive.definite(sigma0+lambda0.min*s0)
        lambda0.step <- 2*lambda0.step
      }
      d <- lambda0.min - lambda0.last
      while (d > tol) {
        temp <- (lambda0.min+lambda0.last)/2
        def_pos <- is.positive.definite(sigma0 + temp*s0)
        if (def_pos) {
          lambda0.min <- temp
        } else {
          lambda0.last <- temp
        }
        d <- lambda0.min - lambda0.last
      }

      if (lambda0.min > lambda0.max)
        stop("'lambda0[2]' too small to make the variance-covariance positive definite")

      f0 <- function(lambda) { 
        sigma0_ridge <- sigma0+lambda*s0
        logLik0 <- logLikPi0 + sum(my.dmvnorm(x, mu0, sigma0_ridge))
        -2*logLik0 + penalty(N, pp)*df(y=y, sigma=sigma0_ridge)
      }
      lambda0 <- optimize(f0, lower=lambda0.min, upper=lambda0.max)$minimum
    }
    sigma0Ridge <- sigma0+lambda0*s0
  } else {
    sigmaRidge <- sigma
    lambda <- 0
    sigma0Ridge <- sigma0
    lambda0 <- 0
    if (sigma==0 | is.na(sigma))
      stop("only one variable remained with 0 or NA sample variance under H_1")
    if (sigma0==0 | is.na(sigma0))
      stop("only one variable remained with 0 or NA sample variance under H_0")
  }
  
  #####################################################
  # Likelihood
  #####################################################
  logLik <- logLikPi
  for (k in 1:K) {    
    logLik <- logLik + sum(my.dmvnorm(xlist[[k]], mu[k,], sigmaRidge))
  }
  logLik0 <- logLikPi0 + sum(my.dmvnorm(x, mu0, sigma0Ridge))

  #####################################################
  # Test statistic
  #####################################################  
  statistic <- -2 * (logLik0 - logLik)
  asymptotic <- 1 - pchisq(statistic, 2*pp*(K-1))

  #####################################################
  # Permutation test
  #####################################################  
  if (p.value.perm) {
    if (fixed.lambda) {
      ll <- lambda
    } else {
      ll <- c(ll.min, lambda.max)
    }
    permutation <- scMANOVApermTest(x=xx, n=n, lambda=ll, lambda0=lambda0, ident=ident, B=B, tol=tol, penalty=penalty, only.pvalue=FALSE, rm.vars=rm.vars, ...)
    p.value <- permutation$p.value
    names(p.value) <- "permutation"
  } else {
    permutation <- NULL
    p.value <- asymptotic
    names(p.value) <- "asymptotic"
  }
  
  res$sigmaRidge <- sigmaRidge
  res$sigma0Ridge <- sigma0Ridge
  res$logLik <- logLik
  res$logLik0 <- logLik0
  res$lambda <- lambda
  res$lambda0 <- lambda0
  res$df <- df(y=y, sigma=sigmaRidge)
  res$df0 <- df(y=y, sigma=sigma0Ridge) 
  res$aic <- -2*logLik + penalty(N, pp)*res$df
  res$aic0 <- -2*logLik0 + penalty(N, pp)*res$df0
  res$statistic <- statistic
  res$p.value <- p.value
  res$permutation <- permutation
  class(res) <- c("scMANOVA", class(res))
  return(res)
}

print.scMANOVA <- function(x, prop=FALSE, mean=FALSE, sigma=FALSE, sigma.ridge=FALSE, digits = getOption("digits"), ...) {
  cat("semicontMANOVA estimates\n\n")
  cat("Number of groups: ", nrow(x$mu), "\n")
  if (!is.null(x$rm.vars)) {
    cat("Number of removed variables: ", length(x$rm.vars), "\n\n")
  } else {
    cat("\n")
  }
  cat("Estimated regularized parameter\n")
  cat("(H1)", format(x$lambda, digits = digits), "df =", format(x$df, digits = digits), "AIC =", format(x$aic, digits = digits), "\n")
  cat("(H0)", format(x$lambda0, digits = digits), "df =", format(x$df0, digits = digits), "AIC =", format(x$aic0, digits = digits), "\n\n")
  cat("log Likelihood \n")
  cat("(H1)", format(x$logLik, digits = digits), "\n")
  cat("(H0)", format(x$logLik0, digits = digits), "\n")
  cat("Wilks test", format(x$statistic, digits = digits), "\n")
  cat("p-value (", names(x$p.value), ")", format(x$p.value, digits = digits), "\n\n")    
  if (prop) {
    dimnames(x$pi) <- list(paste0("Group", 1:nrow(x$pi)), paste0("Var", 0:(ncol(x$pi)-1)))
    names(x$pi0) <- paste0("Var", 0:(length(x$pi0)-1))
    cat("Proportions\n")
    cat("(H1)\n")
    print(x$pi, digits = digits)
    cat("(H0)\n")
    print(x$pi0, digits = digits)
    cat("\n")
  }
  if (mean) {
    dimnames(x$mu) <- list(paste0("Group", 1:nrow(x$mu)), paste0("Var", 1:ncol(x$mu)))
    names(x$mu0) <- paste0("Var", 1:length(x$mu0))
    cat("Mean estimates\n")
    cat("(H1)\n")
    print(x$mu, digits = digits)
    cat("(H0)\n")
    print(x$mu0, digits = digits)
    cat("\n")
  }
  if (sigma) {
    dimnames(x$sigma) <- list(paste0("Var", 1:nrow(x$sigma)), paste0("Var", 1:ncol(x$sigma)))
    dimnames(x$sigma0) <- list(paste0("Var", 1:nrow(x$sigma0)), paste0("Var", 1:ncol(x$sigma0)))
    cat("Sigma estimates\n")
    cat("(H1)\n")
    print(x$sigma, digits = digits)
    cat("(H0)\n")
    print(x$sigma0, digits = digits)
    cat("\n")
  }
  if (sigma.ridge) {
    dimnames(x$sigmaRidge) <- list(paste0("Var", 1:nrow(x$sigmaRidge)), paste0("Var", 1:ncol(x$sigmaRidge)))
    dimnames(x$sigma0Ridge) <- list(paste0("Var", 1:nrow(x$sigma0Ridge)), paste0("Var", 1:ncol(x$sigma0Ridge)))
    cat("Sigma Ridge estimates\n")
    cat("(H1)\n")
    print(x$sigmaRidge, digits = digits)
    cat("(H0)\n")
    print(x$sigma0Ridge, digits = digits)
    cat("\n")
  }
  invisible(x)
}

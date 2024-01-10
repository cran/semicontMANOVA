#############################################################
#	estimation functions
#	Authors: E. Sabbioni, C. Agostinelli and A. Farcomeni
#	Maintainer e-mail: elena.sabbioni@polito.it
#	Date: 25 September 2023
#	Version: 0.1-4
#	Copyright (C) 2023 E. Sabbioni, C. Agostinelli and A. Farcomeni
#############################################################
    
scMANOVAestimation <- function(x, n, lambda=NULL, lambda0=NULL, ident=TRUE, posdef.check=TRUE, rm.vars=NA) {
  ## rm.vars:
  ##   NA    : discard variables
  ##   NULL or a vector of 0 length : keep all the variables
  ##   vector: position of the variables to remove, no further variables are removed
  if (is.null(lambda))
    stop("'lambda' must be a scalar value")
  if (is.null(lambda0))
    stop("'lambda0' must be a scalar value")
  K <- length(n)
  p <- ncol(x)
  N <- sum(n)

  x[x==0] <- NA
  y <- is.na(x) * -1 + 1
  gr <- rep(1:K, times=n)
  xlist <- ylist <- list()
  for (k in 1:K) {
    pos <- which(gr==k)
    xlist[[k]] <- x[pos,,drop=FALSE]
    ylist[[k]] <- y[pos,,drop=FALSE]
  }
  
  ######################################
  # Pi estimation without restrictions
  ######################################
  pi <- rep(NA,K*(p+1))
  ng <- rep(NA,K*(p+1))
  logLikPi <- 0
  for (k in 1:K) {
    jnk <- apply(xlist[[k]], 1, function(z) sum(!is.na(z)))
    jnk <- table(factor(jnk, levels=0:p))
    pi[(p + 1) * (k - 1) + 0:p + 1] <- tmp <- jnk/(choose(p, 0:p) * n[k])
    ng[(p + 1) * (k - 1) + 0:p + 1] <- jnk
    logLikPi <- logLikPi + sum((jnk*log(tmp))[jnk > 0])
  }

  ######################################
  # Pi estimation under null hypothesis
  ######################################
  ng0 <- apply(x, 1, function(z) sum(!is.na(z)))
  ng0 <- table(factor(ng0, levels=0:p))
  pi0 <- ng0/(choose(p, 0:p) * N)
  tmp <- ng0 > 0
  logLikPi0 <- sum((ng0*log(pi0))[tmp])
  
  ######################################
  # Mu estimation without restrictions
  ######################################
  # be careful, mu is a pxK matrix!
  mu <- sapply(xlist, function(z) colMeans(z, na.rm=TRUE))
  ##mu[is.na(mu)] <- 0 # Should this be included?
  
  ######################################
  # Mu estimation under the null hypothesis
  ######################################
  # be careful, mu0 is a vector!
  mu0 <- colMeans(x, na.rm = TRUE)
  ##mu0[is.na(mu0)] <- 0 # Should this be included?

  #####################################################
  # Sigma estimation without restrictions
  #####################################################
  sigma <- matrix(0, p, p)
  for (k in 1:K) {
    stry <- (t(xlist[[k]]) - mu[,k])
    stry[is.na(stry)] <- 0
    sigma <- sigma + stry%*%t(stry)
  }
  yy <- t(y)%*%y
  sigma <- sigma/yy

  if (length(rm.vars) > 0 || (length(rm.vars)==1 && is.na(rm.vars))) {
    rm.vars <- which(is.nan(diag(sigma)) | diag(sigma)==0)
    if (length(rm.vars) > 0) {
      sss <- sigma[-rm.vars,-rm.vars]
      keep <- (1:ncol(sigma))[-rm.vars]  
    } else {
      sss <- sigma
      keep <- 1:ncol(sigma)
    }
    while (any(is.nan(sss))) {
      a <- which.max(apply(sss, 1, function(x) sum(is.nan(x))))
      rm.vars <- c(rm.vars, keep[a])
      keep <- keep[-a]
      sss <- sss[-a, -a]  
    }
    rm.vars <- c(rm.vars, unique(unlist(apply(mu, 2, function(z) which(is.na(z))))))
    rm.vars <- unique(as.vector(rm.vars))
    if (length(rm.vars)==0)
      rm.vars <- NULL
  }
  if (length(rm.vars)==p)
    stop("All variables are removed")
  if (length(rm.vars) > 0) {
    p <- p - length(rm.vars)
    x <- x[, -rm.vars, drop=FALSE]
    y <- y[, -rm.vars, drop=FALSE]
#    temp <- apply(y, 1, function(z) sum(is.na(z))==p)
#    x <- x[!temp,,drop=FALSE]
#    y <- y[!temp,,drop=FALSE] 
    mu <- mu[-rm.vars,, drop=FALSE]
    sigma <- sigma[-rm.vars, -rm.vars, drop=FALSE]
    for (k in 1:K) {
      xlist[[k]] <- xlist[[k]][, -rm.vars, drop=FALSE]
#      temp <- apply(xlist[[k]], 1, function(z) sum(is.na(z))==p)
#      xlist[[k]] <- xlist[[k]][!temp,,drop=FALSE]
    }
#    n <- sapply(xlist, nrow) 
  }
  if (any(n==0)) {
    stop('After removing observations one or more groups are empty') 
  }
  #####################################################
  # Sigma estimation under restrictions
  #####################################################
  if (length(rm.vars) > 0) {
    mu0 <- mu0[-rm.vars]
  }
  sigma0 <- t(x) - mu0
  sigma0[is.na(sigma0)] <- 0 #are the entries with NA on x
  sigma0 <- (sigma0%*%t(sigma0))/(t(y)%*%y)

  #####################################################
  # Ridge (with lambda, lambda0 given by the user)
  #####################################################
  s <- if (ident) {
    diag(p)
  } else {
    diag(diag(sigma), p)
  }  
  sigmaRidge <- sigma+lambda*s
  if (posdef.check) {
    def_pos <- is.positive.definite(sigmaRidge)
    if (!def_pos)
      stop("the ridge estimated variance matrix is not positive definite")
  }
  s0 <- if (ident) {
    diag(p)
  } else {
    diag(diag(sigma0), p)
  }  
  sigma0Ridge <- sigma0+lambda0*s0
  if (posdef.check) {
    def_pos <- is.positive.definite(sigma0Ridge)
    if (!def_pos)
      stop("the ridge estimated variance matrix under the null hypothesis is not positive definite")
  }

  #####################################################
  # Likelihood
  #####################################################
  logLik <- logLikPi
  for (k in 1:K) {    
    logLik <- logLik + sum(my.dmvnorm(xlist[[k]], mu[,k], sigmaRidge))
  }
  logLik0 <- logLikPi0 + sum(my.dmvnorm(x, mu0, sigma0Ridge))

  #####################################################
  # Test statistic
  #####################################################  
  statistic <- -2 * (logLik0 - logLik)
    
  result <- list(pi=matrix(pi, nrow=K, byrow=TRUE), mu=t(mu), sigmaRidge=sigmaRidge, sigma=sigma, pi0=pi0, mu0=mu0, sigma0Ridge=sigma0Ridge, sigma0=sigma0, removed.vars=rm.vars, logLikPi=logLikPi, logLik=logLik, logLikPi0=logLikPi0, logLik0=logLik0, statistic=statistic)
  class(result) <- "scMANOVAestimation"
  return(result)
}

print.scMANOVAestimation <- function(x, prop=FALSE, mean=FALSE, sigma=FALSE, sigma.ridge=FALSE, digits = getOption("digits"), ...) {
  cat("semicontMANOVA estimates\n\n")
  cat("Number of groups: ", nrow(x$mu), "\n")
  if (!is.null(x$removed.vars)) {
    cat("Number of removed variables: ", length(x$removed.vars), "\n\n")
  } else {
    cat("\n")
  }
  cat("log Likelihood \n")
  cat("(H1)", format(x$logLik, digits = digits), "\n")
  cat("(H0)", format(x$logLik0, digits = digits), "\n")
  cat("Wilks test", format(x$statistic, digits = digits), "\n\n")
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

scMANOVAridge <- function(object, lambda=NULL, lambda0=NULL, lambda.step=0.01, ident=TRUE, tol=1e-8) {
  if (!inherits(object, "scMANOVAestimation"))
    stop("'object' must be of class 'scMANOVAestimation'")
  sigma <- x$sigma
  sigma0 <- x$sigma0  
  p <- nrow(sigma)
  lambda0.step <- lambda.step
  ## H1
  s <- if (ident) {
    diag(p)
  } else {
    diag(diag(sigma), p)
  }
  if (is.null(lambda)) {
    lambda <- 0  
    def_pos <- is.positive.definite(sigma+lambda*s)
    while (!def_pos) {
      lambda.last <- lambda
      lambda <- lambda + lambda.step
      def_pos <- is.positive.definite(sigma+lambda*s)
      lambda.step <- 2*lambda.step
    }
    d <- lambda - lambda.last
    while (d > tol) {
      temp <- (lambda+lambda.last)/2
      def_pos <- is.positive.definite(sigma+temp*s)
      if (def_pos) {
        lambda <- temp
      } else {
        lambda.last <- temp
      }
      d <- lambda - lambda.last
    }
  }
  sigmaRidge <- sigma+lambda*s

  ## H_0
  s0 <- if (ident) {
    diag(p)
  } else {
    diag(diag(sigma0), p)
  }
  if (is.null(lambda0)) {
    def_pos <- is.positive.definite(sigma0+lambda0*s0)
    while (!def_pos) {
      lambda0.last <- lambda0
      lambda0 <- lambda0 + lambda0.step
      def_pos <- is.positive.definite(sigma0+lambda0*s0)
      lambda0.step <- 2*lambda0.step
    }
    d <- lambda0 - lambda0.last
    while (d > tol) {
      temp <- (lambda0+lambda0.last)/2
      def_pos <- is.positive.definite(sigma0+temp*s0)
      if (def_pos) {
        lambda0 <- temp
      } else {
        lambda0.last <- temp
      }
      d <- lambda0 - lambda0.last
    }
  }
  sigma0Ridge <- sigma0+lambda0*s0
  x$sigmaRidge <- sigmaRidge
  x$sigma0Ridge <- sigma0Ridge
  class(x) <- c("scMANOVAridge", class(x))
  return(x)  
}

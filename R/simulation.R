#############################################################
#	scMANOVAsimulation function
#	Authors: E. Sabbioni, C. Agostinelli and A. Farcomeni
#	Maintainer e-mail: elena.sabbioni@polito.it
#	Date: 03 January 2024
#	Version: 0.2
#	Copyright (C) 2024 E. Sabbioni, C. Agostinelli and A. Farcomeni
#############################################################
  
scMANOVAsimulation <- function(n, p, pmiss=0, rho=0, mu=NULL, sigma=NULL, only.data=TRUE) {
  K <- length(n) # number of groups
  orig <- matrix(NA, nrow=sum(n), ncol=p) # original dataset, without missing values
  if (is.null(mu)) {
    mu <- matrix(0, nrow=K, ncol=p)
  } else if (is.vector(mu)) {
    mu <- matrix(mu, nrow=K, ncol=p, byrow=TRUE)
  } else if (is.matrix(mu)) {
    if (ncol(mu)!=p)
      stop("cols of 'mu' must be equals to 'p', the number of variables")
    if (nrow(mu)!=K)
      stop("rows of 'mu' must be equals to the length of 'n', the number of groups")
  } else {
    stop("'mu' must be a vector, a matrix or 'NULL'")
  }
  if (is.null(sigma)) {  
    sigma <- matrix(rho, nrow=p, ncol=p)
    diag(sigma) <- 1
  } else if (is.matrix(sigma)) {
    if (nrow(sigma)!=p)
      stop("rows of 'sigma' must be equals to 'p'")
    if (ncol(sigma)!=nrow(sigma))
      stop("rows and cols of 'sigma' must be equals")
  } else {
    stop("'sigma' must be a matrix or 'NULL'")
  }
  if (length(pmiss)!=K)
    pmiss <- rep(pmiss, length.out=K)
  
  # generate data without missing values from a multivariate normal distribution
  index <- 1
  for (k in 1:K) {
    t <- rmvnorm(n[k], mu[k,], sigma)
    orig[index:sum(n[1:k]),] <- t
    index <- index + n[k]
  }
  x <- orig

  # introduce missingness
  y <- matrix(1L, nrow=sum(n), ncol=p)
  index <- 1
  for (k in 1:K) {
    y[index:sum(n[1:k]),] <- matrix(1L-rbinom(n[k]*p, 1, pmiss[k]), n[k], p)
    index <- index + n[k]
  }
  x[which(y==0L)] <- 0 # in x we have missing (0) when y is 0

  # remove components with all missing values
  sumComp <- apply(y, 2, sum)
  x <- x[,which(sumComp!=0)]
  if (only.data)
    return(x)
  else {
    res <- list(x=x, y=y, original=orig, mu=mu, sigma=sigma, n=n, p=p, pmiss=pmiss)
    return(res)
  }
}

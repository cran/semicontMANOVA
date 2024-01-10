scMANOVAloglik <- function(lambda, lambda0, x, n, object, ident=TRUE) {
  if (!inherits(object, "scMANOVAestimation"))
    stop("'object' must be of class 'scMANOVAestimation'")
  K <- length(n)
  p <- ncol(x)
  N <- sum(n)
  x[x==0] <- NA
  y <- is.na(x) * -1 + 1
  rm.vars <- object$removed.vars
  if (length(rm.vars) > 0) {
    p <- p - length(rm.vars)
    x <- x[, -rm.vars, drop=FALSE]
    y <- y[, -rm.vars, drop=FALSE]
  }
  gr <- rep(1:K, times=n)
  xlist <- ylist <- list()
  for (k in 1:K) {
    pos <- which(gr==k)
    xlist[[k]] <- x[pos,,drop=FALSE]
    ylist[[k]] <- y[pos,,drop=FALSE]
  }
  
  mu <- object$mu
  sigma <- object$sigma
  mu0 <- object$mu0
  sigma0 <- object$sigma0
  logLikPi <- object$logLikPi
  logLikPi0 <- object$logLikPi0
  s <- if (ident) {
    diag(p)
  } else {
    diag(diag(sigma), p)
  }  
  s0 <- if (ident) {
    diag(p)
  } else {
    diag(diag(sigma0), p)
  }  

  logTemp <- rep(NA, length(lambda))
  for (i in 1:length(lambda)) {
    logTemp[i] <- logLikPi
    for (k in 1:K) {    
      logTemp[i] <- logTemp[i] + sum(my.dmvnorm(xlist[[k]], mu[,k], sigma+lambda[i]*s))
    }
  }
  logLik0 <- logLikPi0 + sum(my.dmvnorm(x, mu0, sigma0+lambda0*s0))
  result <- list(logLik=logTemp, lambda=lambda, logLik0=logLik0, lambda0=lambda0)
  return(result)
}

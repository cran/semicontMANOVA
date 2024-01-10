#############################################################
#	permutation test functions
#	Authors: E. Sabbioni, C. Agostinelli and A. Farcomeni
#	Maintainer e-mail: elena.sabbioni@polito.it
#	Date: 03 January 2024
#	Version: 0.2-3
#	Copyright (C) 2023 E. Sabbioni, C. Agostinelli and A. Farcomeni
#############################################################

scMANOVApermTest <- function(x, n, lambda=NULL, lambda0=NULL, lambda.step=0.1, ident=FALSE, tol=1e-8, penalty=function(n, p) log(n), B=500, parallel = c("no", "multicore", "snow"), ncpus=1L, cl=NULL, only.pvalue=TRUE, rm.vars=NA, ...) {
  dname <- paste(deparse1(substitute(x)))
  if (is.null(lambda))
    stop("'lambda' must be a scalar value or a vector of length 2")
  if (is.null(lambda0))
    stop("'lambda0' must be a scalar value")
  N <- sum(n)
  if (missing(parallel)) 
    parallel <- "no"
  parallel <- match.arg(parallel)
  have_mc <- have_snow <- FALSE
  if (parallel != "no" && ncpus > 1L) {
    if (parallel == "multicore") 
      have_mc <- .Platform$OS.type != "windows"
    else if (parallel == "snow") 
      have_snow <- TRUE
    if (!have_mc && !have_snow) 
      ncpus <- 1L
    loadNamespace("parallel")
  }
  if (length(lambda)==1) {
    res0 <- scMANOVAestimation(x, n, lambda=lambda, lambda0=lambda0, ident=ident, posdef.check=TRUE, rm.vars=rm.vars)
  } else {
    res0 <- scMANOVA_H1(x, n, lambda=lambda, lambda0=lambda0, lambda.step=lambda.step, ident=ident, tol=tol, penalty=penalty, rm.vars=rm.vars)
  }
  
  if (ncpus > 1L && (have_mc || have_snow)) {
    if (length(lambda)==1) {
      fn <- function(i) {
        xn <- x[sample(N),]
        scMANOVAestimation(xn, n, lambda=lambda, lambda0=lambda0, ident=ident, posdef.check=FALSE, rm.vars=res0$removed.vars)$statistic
      }
    } else {
      fn <- function(i) {
        xn <- x[sample(N),]
        scMANOVA_H1(xn, n, lambda=lambda, lambda0=lambda0, lambda.step=lambda.step, ident=ident, tol=tol, penalty=penalty, rm.vars=res0$removed.vars)$statistic
      }
    }
    if (have_mc) {
      ts <- unlist(parallel::mclapply(X=1:B, FUN=fn, mc.cores = ncpus, ...))
    } else if (have_snow) {
      if (is.null(cl)) {
        cl <- parallel::makePSOCKcluster(rep("localhost", ncpus))
        ts <- parallel::clusterMap(cl=cl, fun=fn, 1:B, RECYCLE=FALSE, SIMPLIFY=TRUE)
        parallel::stopCluster(cl)
      } else
        ts <- parallel::clusterMap(cl=cl, fun=fn, 1:B, RECYCLE=FALSE, SIMPLIFY=TRUE)
    }
  } else {
    ts <- rep(NA, B)
    for (i in 1:B) {
      xn <- x[sample(N),]
      if (length(lambda)==1) {
        ts[i] <- scMANOVAestimation(xn, n, lambda=lambda, lambda0=lambda0, ident=ident, posdef.check=FALSE, rm.vars=res0$removed.vars)$statistic
      } else {
        ts[i] <- scMANOVA_H1(xn, n, lambda=lambda, lambda0=lambda0, lambda.step=lambda.step, ident=ident, tol=tol, penalty=penalty, rm.vars=res0$removed.vars)$statistic
      }
    }
  }
  ts <- c(res0$statistic, ts)
  p.value <- mean(ts>=ts[1], na.rm=TRUE)
  if (only.pvalue) {
    return(p.value)
  } else {
    names(res0$statistic) <- "Wilks"
    dimnames(res0$mu) <- list(paste0("Group", 1:nrow(res0$mu)), paste0("Var", 1:ncol(res0$mu)))
    result <- list(statistic=res0$statistic, p.value=p.value, estimate=res0$mu, null.value=res0$mu0, permutation=ts, estimation=res0, method="semicontMANOVA permutation test", data.name=dname, B=sum(!is.na(ts)))
    class(result) <- "htest"
    return(result)
  }
}

  
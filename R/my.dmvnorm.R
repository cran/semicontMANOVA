my.dmvnorm <- function(x, mean, sigma, log=TRUE) {
  p <- ncol(x)
  n <- nrow(x)
  p <- length(mean)
  if (p==1) {
    res <- drop(dnorm(x, mean, sigma))
    res <- res[!is.na(res)]
  } else {
    y <- !is.na(x)
    res <- rep(0, n)
    for (i in 1:n) {
      tmp <- y[i,]
      if (any(tmp==TRUE)) {
        res[i] <- dmvnorm(x[i,tmp], mean=mean[tmp], sigma=sigma[tmp,tmp,drop=FALSE], log=log)
      }
    }
  }
  res
}

df <- function(y=NULL, n=NULL, sigma) {
  pp <- nrow(sigma)
  if (pp!=1) {
    if (is.null(y)) {
      if (is.null(n)) {
        stop("you have to provide either 'y' or 'n'")
      } else {
        df <- sum(diag(solve(sigma)))*n
        return(df)
      }
    } else {
      df <- 0
      for (i in 1:nrow(y)) {
        if (sum(!y[i,]) < pp) {
          df <- df+sum(diag(solve(sigma[y[i,],y[i,]])))
        }
      }
      return(df)
    }
  } else {
    if (is.null(y)) {
      if (is.null(n)) {
        stop("you have to provide either 'y' or 'n'")
      } else {
        df <- n/sigma
        return(df)
      }
    } else {
      df <- nrow(y)/sigma
      return(df)
    }
  }
}

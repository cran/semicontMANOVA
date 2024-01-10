\name{scMANOVAestimation}
\alias{scMANOVAestimation}
\title{
  Multivariate ANalysis Of VAriance Maximum Likelihood Estimation with Ridge Regularization for Semicontinuous High-Dimensional Data
}
\description{
  \code{scMANOVAestimation} computes the regularized Multivariate ANalysis Of
  VAriance (MANOVA) maximum likelihood estimates for semicontinuous
  high-dimensional data. The estimation can be performed also for
  low-dimensional data. The regularization parameters are provided as input
  and the user can decide to perform the regularization adding the identity
  matrix to the raw estimated covariance matrix (default, \code{ident=TRUE})
  or adding the diagonal values of the raw estimated covariance matrix
  (\code{ident=FALSE}). 
}
\usage{
  scMANOVAestimation(x, n, lambda = NULL, lambda0 = NULL,
    ident = TRUE, posdef.check = TRUE, rm.vars = NA)
}

\arguments{
  \item{x}{\code{data.frame} or \code{matrix} of data with units on the rows and variables on the columns 
}
  \item{n}{\code{vector}. The length corresponds to the number of groups, the elements to the number of observations in each group  
}
  \item{lambda}{scalar. Ridge regularization parameter
}
  \item{lambda0}{scalar. Ridge regularization parameter under null hypothesis
}
  \item{ident}{\code{logical}. If \code{TRUE}, \code{lambda} times the identity matrix is added to the raw estimated covariance matrix, if \code{FALSE} the diagonal values of the raw estimated covariance matrix are used instead
}
  \item{posdef.check}{\code{logical}. Check if the estimated covariance matrix is positive definite
}
  \item{rm.vars}{\code{NA}, \code{NULL} or \code{vector}. If \code{NA} variables are removed from the analysis when they do not have enough observations to compute covariances. If \code{NULL} or a zero length \code{vector} all the variables are retained. If it is a \code{vector}, it indicates the position of the variables to remove, no further variables are removed  
}
}
  
\value{
An object of \code{class} scMANOVAestimation which is a list with the following components  
 \item{pi}{\code{matrix}. Estimated proportion of missing values for each group}
 \item{mu}{\code{matrix}. Estimated mean vector for each group}
 \item{sigmaRidge}{\code{matrix}. Estimated covariance matrix with ridge regularization}
 \item{sigma}{\code{matrix}. Estimated covariance matrix by maximum likelihood}
 \item{pi0}{\code{vector}. Estimated proportion of missing values under the null hypothesis}
 \item{mu0}{\code{vector}. Estimated mean vector under the null hypothesis}
 \item{sigma0Ridge}{\code{matrix}. Estimated covariance matrix with ridge regularization under null hypothesis}
 \item{sigma0}{\code{matrix}. Estimated covariance matrix by maximum likelihood under null hypothesis}
 \item{removed.vars}{\code{vector} or \code{NULL}. columns removed in the continuous part of the log-likelihood dues to insufficient number of observations in each group}
 \item{logLikPi}{scalar. Log-likelihood for the discrete part of the model}
 \item{logLik}{scalar. Log-likelihood}
 \item{logLikPi0}{scalar. Log-likelihood for the discrete part of the model under the null hypothesis}
 \item{logLik0}{scalar. Log-likelihood under null hypothesis}
}
  
\references{
  Elena Sabbioni, Claudio Agostinelli and Alessio Farcomeni (2024)
  A regularized MANOVA test for semicontinuous high-dimensional data.
  arXiv: \url{http://arxiv.org/abs/2401.04036}
}
  
\author{
  Elena Sabbioni, Claudio Agostinelli and Alessio Farcomeni
}

\seealso{\code{\link{scMANOVA}} and \code{\link{scMANOVApermTest}}}


  
\examples{
  set.seed(1234)
  n <- c(5,5)
  p <- 20
  pmiss <- 0.1
  x <- scMANOVAsimulation(n=n, p=p, pmiss=pmiss)
  res <- scMANOVAestimation(x=x, n=n, lambda=3.59, lambda0=3.13)
  res
}

\keyword{multivariate}
\keyword{regression}  
\concept{manova}
\concept{regularization}
\concept{ridge}

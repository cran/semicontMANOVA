\name{scMANOVApermTest}
\alias{scMANOVApermTest}
\title{
   Multivariate ANalysis Of VAriance log-likelihood Test with Ridge Regularization for Semicontinuous High-Dimensional Data
}
\description{
  \code{scMANOVApermTest} uses a permutation procedure to perform a test
  based on a Multivariate ANalysis Of VAriance(MANOVA) Likelihood Ratio test statistic with a ridge
  regularization. The statistic is developed for semicontinuous and
  high-dimensional data, but can be used also in low-dimensional scenarios.
}
\usage{
scMANOVApermTest(x, n, lambda = NULL, lambda0 = NULL, lambda.step = 0.1,
  ident = FALSE, tol = 1e-08, penalty = function(n, p) log(n), B = 500,
  parallel = c("no", "multicore", "snow"), ncpus = 1L, cl = NULL,
  only.pvalue = TRUE, rm.vars = NA, ...)
}  
\arguments{
 \item{x}{\code{data.frame} or \code{matrix} of data with units on the rows and variables on the columns 
}
  \item{n}{\code{vector}. The length corresponds to the number of groups, the elements to the number of observations in each group  
}
  \item{lambda}{scalar or a \code{vector} of length 2. Ridge regularization parameter. The optimal value of \code{lambda} is searched in the specified interval when it is a vector of length 2, otherwise it is used as the optimal value 
}
  \item{lambda0}{\code{NULL}, a scalar or a \code{vector} of length 2. Ridge regularization parameter under null hypothesis. The optimal value of \code{lambda0} is searched in the specified interval when it is a vector of length 2, otherwise it is used as the optimal value
}
  \item{lambda.step}{scalar. Step size used in the optimization procedure to find the smallest value of \code{lambda} (and \code{lambda0}) that makes the covariance matrices,  under the alternative and under the null hypothesis, non singular
}
  \item{ident}{\code{logical}. If \code{TRUE}, \code{lambda} times the identity matrix is added to the raw estimated covariance matrix, if \code{FALSE} the diagonal values of the raw estimated covariance matrix are used instead
}
  \item{tol}{scalar. Used in the optimization procedure to find the smallest value of \code{lambda} (and \code{lambda0}) that makes the covariance matrices, under the alternative and under the null hypothesis, non singular 
}
  \item{penalty}{\code{function} with two arguments: sample size (\code{n}) and number of variables (\code{p}) used as penalty function in the definition of the Information Criterion to select the optimal values for \code{lambda} and \code{lambda0}
}
  \item{B}{scalar. Number of permutations to run in the permutation test 
}
  \item{parallel}{The type of parallel operation to be used (if any)
}
  \item{ncpus}{\code{integer}. Number of processes to be used in parallel operation: typically one would chose this to the number of
          available CPUs.
}
  \item{cl}{An optional \code{parallel} or \code{snow} cluster to use if \code{parallel = "snow"}. If not supplied, a cluster on the local
     machine is created for the duration of the call
}
  \item{only.pvalue}{\code{logical}. If \code{TRUE} only the p-value is returned
}
  \item{rm.vars}{\code{vector}. It indicates the position of the variables to remove
}
  \item{\dots}{Further parameters passed to \code{parallel::mclapply} in case of \code{parallel="multicore"}
}
}
\value{
  If \code{only.pvalue=TRUE} (default) a scalar which is the p-value of the Wilks statistic obtain by a permutation procedure, otherwise an object of \code{class} htest  
}
  
\references{
  Elena Sabbioni, Claudio Agostinelli and Alessio Farcomeni (2024)
  A regularized MANOVA test for semicontinuous high-dimensional data.
  arXiv: \url{http://arxiv.org/abs/2401.04036}
}
  
\author{
  Elena Sabbioni, Claudio Agostinelli and Alessio Farcomeni
}

\seealso{\code{\link{scMANOVA}} and \code{\link{scMANOVAestimation}}}
  
\examples{
  set.seed(1234)
  n <- c(5,5)
  p <- 20
  pmiss <- 0.1
  x <- scMANOVAsimulation(n=n, p=p, pmiss=pmiss)
  res <- scMANOVApermTest(x=x, n=n, lambda=3.59, lambda0=3.13,
    only.pvalue=FALSE)
  res
}
  
\keyword{htest}
\keyword{multivariate}
\keyword{regression}  
\concept{manova}
\concept{regularization}
\concept{ridge}
\concept{permutation}
  
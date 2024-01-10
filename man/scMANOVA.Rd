\name{scMANOVA}
\alias{scMANOVA}
\title{
  Multivariate ANalysis Of VAriance Inference and Test with Ridge Regularization for Semicontinuous High-Dimensional Data
}
\description{
  \code{scMANOVA} performs Multivariate ANalysis Of VAriance (MANOVA) inference and test with ridge regularization in presence of 
  semicontinuous high-dimensional data. The test is based on a Likelihood Ratio Test statistic 
  and the p-value can be computed using either asymptotic distribution (\code{p.value.perm = FALSE}) 
  or via permutation procedure (\code{p.value.perm = TRUE}). There is the possibility to provide 
  as input the regularization parameters or to choose them through an optimization procedure. 
}
\usage{
scMANOVA(x, n, lambda = NULL, lambda0 = NULL, lambda.step = 0.1,
  ident = FALSE, tol = 1e-08, penalty = function(n, p) log(n),
  B = 500, p.value.perm = FALSE, fixed.lambda = FALSE, ...)
}
\arguments{
  \item{x}{\code{data.frame} or \code{matrix} of data with units on the rows and variables on the columns 
}
  \item{n}{\code{vector}. The length corresponds to the number of groups, the elements to the number of observations in each group  
}
  \item{lambda}{\code{NULL}, a scalar or a \code{vector} of length 2. Ridge regularization parameter. The optimal value of \code{lambda} is searched in the interval [0,100] if \code{NULL}, and in the specified interval when it is a vector of length 2, otherwise it is used as the optimal value 
}
  \item{lambda0}{\code{NULL}, a scalar or a \code{vector} of length 2. Ridge regularization parameter under null hypothesis. The optimal value of \code{lambda0} is searched in the interval [0,100] if \code{NULL}, and in the specified interval when it is a vector of length 2, otherwise it is used as the optimal value
}
  \item{lambda.step}{scalar. Step size used in the optimization procedure to find the smallest value of \code{lambda} (and \code{lambda0}) that makes the covariance matrices, under the alternative and under the null hypotheses, non singular
}
  \item{ident}{\code{logical}. If \code{TRUE}, \code{lambda} times the identity matrix is added to the raw estimated covariance matrix, if \code{FALSE} the diagonal values of the raw estimated covariance matrix are used instead
}
  \item{tol}{scalar. Used in the optimization procedure to find the smallest value of \code{lambda} (and \code{lambda0}) that makes the covariance matrices, under the alternative and under the null, non singular 
}
  \item{penalty}{\code{function} with two arguments: sample size (\code{n}) and number of variables (\code{p}) used as penalty function in the definition of the Information Criterion to select the optimal values for \code{lambda} and \code{lambda0}
}
  \item{B}{scalar. Number of permutations to run in the permutation test 
}
  \item{p.value.perm}{\code{logical}. If \code{TRUE} a p-value from a permutation test is evaluated, otherwise an asymptotic value is reported
}
  \item{fixed.lambda}{\code{logical}. If \code{TRUE} the optimal values for \code{lambda} and \code{lambda0} are evaluated just once for the observed dataset and kept fixed during the permutation test, otherwise, optimal values are evaluated for each permuted datsets
}
  \item{\dots}{further parameters passed to function \code{\link{scMANOVApermTest}}
}
}
\value{
An object of \code{class} scMANOVA which is a list with the following components  
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
 \item{statistic}{scalar. Wilks statistics}
 \item{lambda}{scalar. Regularization parameter}
 \item{lambda0}{scalar. Regularization parameter under null hypothesis}
 \item{df}{scalar. Model degree of freedom}
 \item{df0}{scalar. Model degree of freedom under null hypothesis}
 \item{aic}{scalar. Information criteria}
 \item{aic0}{scalar. Information criteria under null hypothesis}
 \item{p.value}{scalar. p-value of the Wilks statistic}
}
  
\references{
  Elena Sabbioni, Claudio Agostinelli and Alessio Farcomeni (2024)
  A regularized MANOVA test for semicontinuous high-dimensional data.
  arXiv: \url{http://arxiv.org/abs/2401.04036}
}
  
\author{
  Elena Sabbioni, Claudio Agostinelli and Alessio Farcomeni
}

\seealso{\code{\link{scMANOVAestimation}} and \code{\link{scMANOVApermTest}}}
  
\examples{
  set.seed(1234)
  n <- c(5,5)
  p <- 20
  pmiss <- 0.1
  x <- scMANOVAsimulation(n=n, p=p, pmiss=pmiss)
  res.asy <- scMANOVA(x=x, n=n) # Asymptotic p.value
  res.asy
  \donttest{
    res.perm <- scMANOVA(x=x, n=n, p.value.perm=TRUE) # p-value by permutation test 
    res.perm
  }
}
  
\keyword{htest}
\keyword{multivariate}
\keyword{regression}  
\concept{manova}
\concept{regularization}
\concept{ridge}
\concept{permutation}

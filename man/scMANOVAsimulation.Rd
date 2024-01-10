\name{scMANOVAsimulation}
\alias{scMANOVAsimulation}
\title{
  Simulation of datasets for a semicontinuous scenarios
}
\description{
  Simulation of dataset of semicontinuous data coming from different groups, 
  with specific marginal probabilities of a missing value, specific mean vectors
  and common covariance matrix.
}
\usage{
scMANOVAsimulation(n, p, pmiss = 0, rho = 0, mu = NULL,
  sigma = NULL, only.data = TRUE)
}
\arguments{
  \item{n}{\code{vector}. The length corresponds to the number of groups, the elements to the number of observations in each group  
}
  \item{p}{scalar. Number of variables (columns)
}
  \item{pmiss}{scalar or \code{vector}. Proportion of missingness in each group. If it is a scalar the same proportion is used in each group
}
  \item{rho}{scalar. If \code{sigma=NULL} then \code{sigma} is set as a covariance matrix with covariance \code{rho} equal in every entries that is not on the main diagonal of \code{sigma}, and variance equal to 1
}
  \item{mu}{\code{NULL} or \code{vector} or \code{matrix}. If \code{NULL} the mean of each group is set zero for all the variables, if \code{vector} the different groups have the same mean. If \code{matrix} each row corresponds to the mean vector of the corresponding group
}
  \item{sigma}{\code{NULL} or \code{matrix}. If \code{matrix} it is a covariance matrix. If \code{NULL} the value of \code{rho} is used to build the covariance matrix
}
  \item{only.data}{\code{logical}. If \code{TRUE} only the simulated data are reported
}
}
\value{
  If \code{only.data=TRUE} an object of \code{class} \code{matrix} is reported otherwise a \code{list} with the following components
  \item{x}{\code{matrix}. The simulated dataset}
  \item{y}{\code{matrix}. A matrix with zero when the corresponding entry in \code{x} is zero and one otherwise}
  \item{original}{\code{matrix}. The simulated dataset without missing values}
  \item{mu}{\code{matrix}. Mean vectors, on for each group}
  \item{sigma}{\code{matrix}. Covariance matric}
  \item{n}{As in input}
  \item{p}{As in input}
  \item{pmiss}{\code{vector}. Proportion of missingness in each group}
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
}
  
\keyword{multivariate}
\keyword{datasets}  
\concept{manova}

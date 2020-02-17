\name{mcgk.dist.test}
\alias{mcgk.dist.test}
\title{
Test of independence among several random vectors
}
\description{
A consistent test for independence among several random variables using copula MMD as a dependence characterizer. It utilizes the inter-point distance method.
}
\usage{
mcgk.dist.test(Xlist = NULL , Dlist = NULL , alpha = 0.05, 
parameters = NULL, B.parameters = 1000, B = 1000,
test.preset = NULL)
}
\arguments{
  \item{Xlist}{
A list matrices where each matrix corresponds to each random vector. Each row of a particular matrix is an observation on the corresponding random vector.
}
  \item{Dlist}{
A list of "dist" objects where each distance object corresponds to each random vector. It is required if \code{Xlist} is \code{NULL}.
}
  \item{alpha}{
A numeric value between 0 and 1 specifying the level of the test.
}
  \item{parameters}{
If \code{NULL}, then \code{mcgk.dist.test} automatically calculates the parameters. Otherwise it uses the given set of parameters.
}
  \item{B.parameters}{
An integer value specifying the number of Monte-Carlo iterations to be used for automatic parameter generation procedure.
}
  \item{B}{
An integer value specifying the number of Monte-Carlo iterations made in the permutation test.
}
  \item{test.preset}{
If the list of numeric vectors is precalculated by the \code{mcgk.test.preset} function, supply that list here.
}
}
\details{
See reference.
}
\value{
A list of following components:
  \item{Tdist.sum.stat}{
The value of the test statistic Tdist.sum.
}
  \item{Tdist.max.stat}{
The value of the test statistic Tdist.max.
}
  \item{Tdist.sum.cutoff}{
The critical value for the test statistic Tdist.sum.
}
  \item{Tdist.max.cutoff}{
The critical value for the test statistic Tdist.max.
}
  \item{Tdist.sum.pvalue}{
P-value of the hypothesis test based on Tdist.sum.
}
  \item{Tdist.max.pvalue}{
P-value of the hypothesis test based on Tdist.max.
}
}
\references{
Roy et al. (2019)
}
\author{
Angshuman Roy
}
\examples{
set.seed(1)
x <- rnorm(10)
y <- rnorm(10)
z <- matrix(x+y+rnorm(10,sd=0.1),ncol=1)
xlist <- list(x,y,z)
set.seed(1)
mcgk.dist.test(xlist,alpha=0.05)

##
set.seed(1)
x <- rnorm(10)
y <- rnorm(10)
z <- matrix(x+y+rnorm(10,sd=0.1),ncol=1)
xlist <- list(x,y,z)
dlist <- lapply(xlist,dist)
set.seed(1)
mcgk.dist.test(Dlist=dlist,alpha=0.05)
}
\keyword{nonparametric}
\keyword{htest}

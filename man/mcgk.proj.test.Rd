\name{mcgk.proj.test}
\alias{mcgk.proj.test}
\title{
Test of independence among several random vectors
}
\description{
A consistent test for independence among several random variables using copula MMD as a dependence characterizer. It utilizes the linear projection method.
}
\usage{
mcgk.proj.test(Xlist = NULL , Ilist = NULL , alpha = 0.05,
parameters = NULL, B.parameters = 1000, B = 1000,
test.preset = NULL)
}
\arguments{
  \item{Xlist}{
A list matrices where each matrix corresponds to each random vector. Each row of a particular matrix is an observation on the corresponding random vector.
}
  \item{Ilist}{
A list of inner product gram matrices where each matrix corresponds to each random vector. It is required if \code{Xlist} is \code{NULL}.
}
  \item{alpha}{
A numeric value between 0 and 1 specifying the level of the test.
}
  \item{parameters}{
If \code{NULL}, then \code{mcgk.proj.test} automatically calculates the parameters. Otherwise it uses the given set of parameters.
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
  \item{Tproj.sum.stat}{
The value of the test statistic Tproj.sum.
}
  \item{Tproj.max.stat}{
The value of the test statistic Tproj.max.
}
  \item{Tproj.sum.cutoff}{
The critical value for the test statistic Tproj.sum.
}
  \item{Tproj.max.cutoff}{
The critical value for the test statistic Tproj.max.
}
  \item{Tproj.sum.pvalue}{
P-value of the hypothesis test based on Tproj.sum.
}
  \item{Tproj.max.pvalue}{
P-value of the hypothesis test based on Tproj.max.
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
mcgk.proj.test(xlist,alpha=0.05)

##
set.seed(1)
x <- rnorm(10)
y <- rnorm(10)
z <- matrix(x+y+rnorm(10,sd=0.1),ncol=1)
xlist <- list(x,y,z)
xlist <- lapply(xlist,as.matrix)
ilist <- lapply(xlist,function(x) crossprod(t(x)))
set.seed(1)
mcgk.proj.test(Ilist=ilist,alpha=0.05)
}
\keyword{nonparametric}
\keyword{htest}

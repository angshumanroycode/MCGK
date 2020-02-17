\name{mcgk.test.preset}
\alias{mcgk.test.preset}
\title{
A list of numeric vectors pre-calculated for \code{mcgk.dist.test} or \code{mcgk.proj.test} function
}
\description{
When \code{mcgk.dist.test} or \code{mcgk.proj.test} function has to be calculated repeatedly, we can save a lot of computation by pre-calculating some numeric vectors that \code{mcgk.dist.test} or \code{mcgk.proj.test} function uses by using \code{mcgk.test.preset} function.
}
\usage{
mcgk.test.preset(n, d, parameters = NULL, B.parameters = 1000)
}
\arguments{
  \item{n}{
Sample size.
}
  \item{d}{
Number of random vectors.
}
  \item{parameters}{
If \code{NULL}, then \code{mcgk.test.preset} automatically calculates the parameters. Otherwise it calculates preset using the given set of parameters.
}
  \item{B.parameters}{
An integer value specifying the number of Monte-Carlo iterations to be used for automatic parameter generation procedure.
}
}
\details{
See reference.
}
\value{
A list of numeric vectors.
}
\references{
Roy et al. (2019)
}
\author{
Angshuman Roy
}
\examples{
#automated parameters
set.seed(0)
mcgk.test.preset(10,2,B.parameters=1000)

#manual parameters
mcgk.test.preset(10,2,parameters=c(0.3,0.5,1))
}
\keyword{nonparametric}

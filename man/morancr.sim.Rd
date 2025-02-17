% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/morancr.sim.R
\name{morancr.sim}
\alias{morancr.sim}
\title{Constant-risk Moran's I statistic}
\usage{
morancr.sim(nsim = 1, cases, w, ex)
}
\arguments{
\item{nsim}{The number of simulations from which to
compute the p-value.}

\item{cases}{The number of cases observed in each region.}

\item{w}{A binary spatial adjacency matrix for the
regions.}

\item{ex}{The expected number of cases for each region.
The default is calculated under the constant risk
hypothesis.}
}
\value{
Returns a numeric value.
}
\description{
\code{morancr.stat} computes the constant-risk version of the Moran's I
statistic proposed by Walter (1992).
}
\examples{
data(nydf)
data(nyw)
ex = sum(nydf$cases)/sum(nydf$pop)*nydf$pop
morancr.sim(nsim = 10, cases = nydf$cases, w = nyw, ex = ex)
}
\references{
Walter, S. D. (1992). The analysis of regional patterns in health
  data: I. Distributional considerations. American Journal of Epidemiology,
  136(6), 730-741.
}
\seealso{
\code{\link{morancr.test}}
}
\author{
Joshua French
}

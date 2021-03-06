% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/product_gev.R
\name{dGEVxGEV}
\alias{dGEVxGEV}
\alias{pGEVxGEV}
\alias{qGEVxGEV}
\alias{rGEVxGEV}
\title{Two-component generalized extreme value distribution (GEV)}
\usage{
dGEVxGEV(x, param1, param2)

pGEVxGEV(q, param1, param2)

qGEVxGEV(p, param1, param2)

rGEVxGEV(n, param1, param2)
}
\arguments{
\item{x}{vector of quantiles.}

\item{param1}{three-dimensional vector (loc, scale, shape)' of a GEV from season 1.}

\item{param2}{three-dimensional vector (loc, scale, shape)' of a GEV from season 2.}

\item{q}{vector of quantiles.}

\item{p}{vector of probabilities.}

\item{n}{number of observations.}
}
\description{
Density, distribution function, quantile function and
random generation for a two-component GEV distribution (product of two GEVs).
}
\details{
These functions use the parametrization of the \link[evd]{gev}-functions from the package 'evd'.
The distribution \eqn{F} of a two-component GEV is:
\eqn{F=F_1 \cdot F_2}{F=F_1 * F_2}, where \eqn{F_1} and \eqn{F_2} are two
distribution functions of a GEV.
}
\examples{
# density and distribution function of a two-component GEV:
par(mfrow=c(3,1))
curve(dGEVxGEV(x, c(2,1,0.2), c(3,2,0.4)), from=0, to=20, ylab="Density", n=1000)
curve(pGEVxGEV(x, c(2,1,0.2), c(3,2,0.4)), from=0, to=20, ylab="Probability", n=1000)

# quantiles of a two-component GEV:
qGEVxGEV(p=c(0.9, 0.99), c(2,1,0.2), c(3,2,0.4))

# random numbers of a two-component GEV:
set.seed(23764)
rn <- rGEVxGEV(1000, c(2,1,0.1), c(3,2,0))
hist(rn, breaks=40, freq=FALSE, main="")
curve(dGEVxGEV(x, c(2,1,0.1), c(3,2,0)), from=0, to=20,
ylab="density", n=1000, col="red", add=TRUE)
}


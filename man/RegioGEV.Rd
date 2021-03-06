% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/regional_gev.R
\name{RegioGEV}
\alias{RegioGEV}
\title{Regional (or local) parameter and quantile estimation}
\usage{
RegioGEV(x, p, j = 1, leftrim = 0, rightrim = 0, na.rm = TRUE, ...)
}
\arguments{
\item{x}{vector or matrix of observations (rows: observations, d columns: stations).}

\item{p}{a probability.}

\item{j}{quantile and parameter estimation for the jth
station (jth column of \code{x}). Irrelevant if is \code{x} is a vector.}

\item{leftrim}{integer indicating lower trimming parameter (\eqn{\ge 0}).}

\item{rightrim}{integer indicating upper trimming parameter (\eqn{\ge 0}).}

\item{na.rm}{Should missing values be removed?}

\item{...}{additional arguments, see \link[TLMoments]{TLMoments}.}
}
\value{
List of \itemize{
\item \code{quant} quantile calculation from an estimated GEV with a regional shape-parameter.
\item \code{param} estimated parameter vector from a GEV (using L-moments or trimmed L-moments).
\item \code{w} optimal or simple weighting (just returned if \code{x} is a matrix).}
}
\description{
Calculates regional (or local) parameters of a generalized extreme value (GEV) distribution
using (trimmed) L-moments
(see \link[TLMoments]{TLMoments} and \link[TLMoments]{parameters}) from a vector or matrix of observation.
Based on these parameters, a p-quantile of the GEV will be calculated for the jth station.
}
\details{
The optimal weights will be calculated as described in "Kinsvater, Fried and Lilienthal (2015):
Regional extreme value index estimation and a test of tail homogeneity,
Environmetrics, DOI: 10.1002/env.2376, Section 3.2". If it's not possible to calculate
optimal weights (negative eigenvaules of an estimated covarinace matrix), simple weights
will be calculated: \eqn{w_j=\frac{n_j}{sum_{j=1}^d n_j}}{w_j=n_j/sum_{j=1}^d n_j}
}
\examples{
library("evd")
# sample observations of 75 years at one station:
x <- rgev(75) # x is a vector
RegioGEV(x=x, p=0.95)

x2 <- c(NA, NA, x[1:60], NA, x[61:75]) # vector of observations with missing values
RegioGEV(x=x2, p=0.95) # NAs will be removed

# sample observations of 100 years at 4 stations:
set.seed(1053)
x <- matrix(rgev(400, 2, 1, 0.3), ncol=4)
RegioGEV(x=x, p=0.9, j=3, leftrim=0, rightrim=0) # optimal weighting
RegioGEV(x=x, p=0.9, j=3, leftrim=0, rightrim=1) # optimal weighting

# With missing values:
x2 <- x
x2[c(54, 89, 300)] <- NA
RegioGEV(x=x2, p=0.9, j=3, leftrim=0, rightrim=0)

# sample again observations of 100 years at 4 stations:
set.seed(958)
x <- matrix(rgev(400, 2, 1, 0.3), ncol=4)
RegioGEV(x=x, p=0.9, j=3, leftrim=0, rightrim=0) # simple weighting
}


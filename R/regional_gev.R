#' @title Regional (or local) parameter and quantile estimation
#' @description Calculates regional (or local) parameters of a generalized extreme value (GEV) distribution
#' using (trimmed) L-moments
#' (see \link[TLMoments]{TLMoments} and \link[TLMoments]{parameters}) from a vector or matrix of observation.
#' Based on these parameters, a p-quantile of the GEV will be calculated for the jth station.
#' @param x vector or matrix of observations (rows: observations, d columns: stations).
#' @param p a probability.
#' @param j quantile and parameter estimation for the jth
#' station (jth column of \code{x}). Irrelevant if is \code{x} is a vector.
#' @param leftrim	integer indicating lower trimming parameter (\eqn{\ge 0}).
#' @param rightrim	integer indicating upper trimming parameter (\eqn{\ge 0}).
#' @param na.rm Should missing values be removed?
#' @param ... additional arguments, see \link[TLMoments]{TLMoments}.
#' @details The optimal weights will be calculated as described in "Kinsvater, Fried and Lilienthal (2015):
#' Regional extreme value index estimation and a test of tail homogeneity,
#' Environmetrics, DOI: 10.1002/env.2376, Section 3.2". If it's not possible to calculate
#' optimal weights (negative eigenvaules of an estimated covarinace matrix), simple weights
#' will be calculated: \eqn{w_j=\frac{n_j}{sum_{j=1}^d n_j}}{w_j=n_j/sum_{j=1}^d n_j}
#' @return List of \itemize{
#' \item \code{quant} quantile calculation from an estimated GEV with a regional shape-parameter.
#' \item \code{param} estimated parameter vector from a GEV (using L-moments or trimmed L-moments).
#' \item \code{w} optimal or simple weighting (just returned if \code{x} is a matrix).}
#' @examples
#' library("evd")
#' # sample observations of 75 years at one station:
#' x <- rgev(75) # x is a vector
#' RegioGEV(x=x, p=0.95)
#'
#' x2 <- c(NA, NA, x[1:60], NA, x[61:75]) # vector of observations with missing values
#' RegioGEV(x=x2, p=0.95) # NAs will be removed
#'
#' # sample observations of 100 years at 4 stations:
#' set.seed(1053)
#' x <- matrix(rgev(400, 2, 1, 0.3), ncol=4)
#' RegioGEV(x=x, p=0.9, j=3, leftrim=0, rightrim=0) # optimal weighting
#' RegioGEV(x=x, p=0.9, j=3, leftrim=0, rightrim=1) # optimal weighting
#'
#' # With missing values:
#' x2 <- x
#' x2[c(54, 89, 300)] <- NA
#' RegioGEV(x=x2, p=0.9, j=3, leftrim=0, rightrim=0)
#'
#' # sample again observations of 100 years at 4 stations:
#' set.seed(958)
#' x <- matrix(rgev(400, 2, 1, 0.3), ncol=4)
#' RegioGEV(x=x, p=0.9, j=3, leftrim=0, rightrim=0) # simple weighting
#' @export
RegioGEV <- function(x, p, j=1, leftrim=0, rightrim=0, na.rm=TRUE, ...){
  if(!all(p < 1 & p > 0)) stop("p must be a probability: p in (0,1)")

  # local case (d=1):
  if(is.vector(x)){
    para <- TLMoments(x, leftrim=leftrim, rightrim=rightrim, na.rm=na.rm, ...) %>% parameters("gev")
    quan <- para %>% quantiles(p)
    return(list(quant=quan, param=para))
  }

  # regional case (d>1):
  if(is.matrix(x)){
    # estimate parameters from (T)L-moments:
    tlmom <- TLMoments(x, leftrim=leftrim, rightrim=rightrim, na.rm=na.rm, ...)
    para <- tlmom %>% parameters("gev")
    d <- ncol(x)

    # calculate weights w with method "xi_cov":

    Sigma_Inv <- solve(est_cov(para)[1:d*3, 1:d*3])
    ev <- eigen(Sigma_Inv)$values

    if(any(ev < 0)){
      n <- apply(x, 2, function(y) sum(!is.na(y))) # winter- and sommermaxima are observed pairwise or NA pairwise.
      w <- n/sum(n) # simple weighting if any eigenvalues < 0
      wname <- "w"
    }

    else{
      w <- rowSums(Sigma_Inv)/sum(Sigma_Inv) # optimal weighting
      wname <- "wopt"
    }

    # calculate the estimated regional xi
    xi <- w %*% para[3,]

    x1 <- x[,j]
    betas <- PWMs(x[,j], na.rm=TRUE)

    # calculate the estimated mu and sigma at site j
    if(leftrim==0 && rightrim==0){
      sig <- (2*betas[2]-betas[1])*xi/(gamma(1-xi)*(2^xi-1))
      mu <- betas[1]+ sig/xi * (1-gamma(1-xi))
    }
    if(leftrim==0 && rightrim==1){
      sig <- (4*betas[2]-betas[1]-3*betas[3])/(gamma(-xi)*(3^xi-2*2^xi+1))
      mu <- 2*(betas[1]-betas[2]) + sig/xi - sig*gamma(-xi)*(2^xi-2)
    }

    # calculate p-quantile
    quan <- as.vector(evd::qgev(p, mu, sig, xi))

    return(list(quant=quan, param=c(loc=mu, scale=sig, shape=xi), w=wname))
  }
}

#' @title Seasonal regional (or local) parameter and quantile estimation
#' @description Calculates regional (or local) parameters of a two-component GEV distribution (product of two GEVs)
#' using (trimmed) L-moments (see \link[TLMoments]{TLMoments} and \link[TLMoments]{parameters})
#' from two vectors or two matrices of observation, e.g. winter and summer observations
#' from one or more than one station.
#' Based on these two parameter vectors, a p-quantile of the two-component GEV will
#' be calculated for the jth station.
#' @param x1 vector or matrix of observations from season 1 (rows: observations, columns: stations).
#' @param x2 vector or matrix of observations from season 2 (rows: observations, columns: stations).
#' @param p a probability.
#' @param j quantile and parameter estimation for the jth station (jth column of \code{x}).
#' Irrelevant if is \code{x1} and \code{x2} are vectors.
#' @param leftrim	integer indicating lower trimming parameter (\eqn{\ge 0}).
#' @param rightrim	integer indicating upper trimming parameter (\eqn{\ge 0}).
#' @param na.rm Should missing values be removed?
#' @param ... additional arguments, see \link[TLMoments]{TLMoments}.
#' @return List of \itemize{
#' \item \code{quant} quantile calculation from an estimated two-component GEV with a
#' regional (or local) shape-parameters.
#' \item \code{param1} estimated parameter vector from season 1 from a GEV (using L-moments or trimmed L-moments).
#' \item \code{param2} estimated parameter vector from season 2 from a GEV (using L-moments or trimmed L-moments).}
#' @examples
#' library("evd")
#' # Seasonal observations of 80 years at one station:
#' x1 <- rgev(80, 2, 1, 0.2) # observations from season 1
#' x2 <- rgev(80, 3, 1, 0.3) # observations from season 2
#' RegioGEVSeas(x1=x1, x2=x2, p=0.95)
#'
#' # Missing values in both seasons in the same year(s):
#' x1a <- c(NA, x1, NA)
#' x2a <- c(NA, x2, NA)
#' RegioGEVSeas(x1a, x2a, p=0.99, leftrim=0, rightrim=0, na.rm=TRUE)
#'
#' # Missing values in both seasons in different year(s):
#' x1b <- x1
#' x1b[c(4,19)] <- NA
#' x2b <- x2
#' x2b[77] <- NA
#' RegioGEVSeas(x1b, x2b, p=0.99, leftrim=0, rightrim=0, na.rm=TRUE)
#'
#' # Seasonal observations of 100 years at 4 stations:
#' x1 <- matrix(rgev(400, 2, 1, 0.3), ncol=4)
#' x2 <- matrix(rgev(400, 4, 1, 0.2), ncol=4)
#' # estimate quantile for station 1 and 2 (consider the same shape-parameters):
#' RegioGEVSeas(x1, x2, p=0.99, j=1, leftrim=0, rightrim=0)
#' RegioGEVSeas(x1, x2, p=0.99, j=2, leftrim=0, rightrim=0)
#'
#' # With missing values:
#' x3 <- x1
#' x4 <- x2
#' x3[c(54, 89, 300)] <- NA
#' RegioGEVSeas(x3, x4, p=0.99, j=1, leftrim=0, rightrim=0)
#' @export
RegioGEVSeas <- function(x1, x2, p, j=1, leftrim=0, rightrim=0, na.rm=TRUE, ...){
  if(!all(p < 1 & p > 0)) stop("p must be a probability: p in (0,1)")

  param1 <- RegioGEV(x1, p=p, j=j, leftrim=leftrim, rightrim=rightrim, na.rm=na.rm, ...)$param
  param2 <- RegioGEV(x2, p=p, j=j, leftrim=leftrim, rightrim=rightrim, na.rm=na.rm, ...)$param
  quan <- qGEVxGEV(p=p, param1, param2)

  return(list(quant=quan, param1=param1, param2=param2))
}




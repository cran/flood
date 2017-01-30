#' @title Two-component generalized extreme value distribution (GEV)
#' @description Density, distribution function, quantile function and
#' random generation for a two-component GEV distribution (product of two GEVs).
#' @param x vector of quantiles.
#' @param q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations.
#' @param param1	three-dimensional vector (loc, scale, shape)' of a GEV from season 1.
#' @param param2	three-dimensional vector (loc, scale, shape)' of a GEV from season 2.
#' @details These functions use the parametrization of the \link[evd]{gev}-functions from the package 'evd'.
#' The distribution \eqn{F} of a two-component GEV is:
#' \eqn{F=F_1 \cdot F_2}{F=F_1 * F_2}, where \eqn{F_1} and \eqn{F_2} are two
#' distribution functions of a GEV.
#' @examples
#' # density and distribution function of a two-component GEV:
#' par(mfrow=c(3,1))
#' curve(dGEVxGEV(x, c(2,1,0.2), c(3,2,0.4)), from=0, to=20, ylab="Density", n=1000)
#' curve(pGEVxGEV(x, c(2,1,0.2), c(3,2,0.4)), from=0, to=20, ylab="Probability", n=1000)
#'
#' # quantiles of a two-component GEV:
#' qGEVxGEV(p=c(0.9, 0.99), c(2,1,0.2), c(3,2,0.4))
#'
#' # random numbers of a two-component GEV:
#' set.seed(23764)
#' rn <- rGEVxGEV(1000, c(2,1,0.1), c(3,2,0))
#' hist(rn, breaks=40, freq=FALSE, main="")
#' curve(dGEVxGEV(x, c(2,1,0.1), c(3,2,0)), from=0, to=20,
#' ylab="density", n=1000, col="red", add=TRUE)
#' @rdname twocompGEV
#' @export
dGEVxGEV <- function(x, param1, param2){
  fs <- evd::dgev(x, loc = param1[1], scale = param1[2], shape = param1[3])
  Fs <- evd::pgev(x, loc = param1[1], scale = param1[2], shape = param1[3])

  fw <- evd::dgev(x, loc = param2[1], scale = param2[2], shape = param2[3])
  Fw <- evd::pgev(x, loc = param2[1], scale = param2[2], shape = param2[3])

  fprodgev <- fs*Fw + fw*Fs
  return(fprodgev)
}

#' @rdname twocompGEV
#' @export
pGEVxGEV <- function(q, param1, param2){
  evd::pgev(q, loc = param1[1], scale = param1[2], shape = param1[3]) * evd::pgev(q, loc = param2[1], scale = param2[2], shape = param2[3])
}

#' @rdname twocompGEV
#' @export
qGEVxGEV <- function(p, param1, param2){
  if(!all(p < 1 & p > 0)) stop("p must be a probability: p in (0,1)")
  f <- function(q) evd::pgev(q, loc = param1[1], scale = param1[2], shape = param1[3]) * evd::pgev(q, loc = param2[1], scale = param2[2], shape = param2[3])
  sapply(p, function(y) uniroot(function(x) f(x)-y, interval=c(0.01, 1e10))$root)
}

#' @rdname twocompGEV
#' @export
rGEVxGEV <- function(n, param1, param2){
  u <- runif(n)
  qGEVxGEV(u, param1, param2)
}




#' @title Block maxima distribution
#' @description Calculates quantiles of a block maxima distribution.
#' @param p vector of probabilities.
#' @param b block length (in general \code{b} \eqn{\ge 2}).
#' @param param three-dimensional vector with location (mu), scale (sigma)
#' and shape (xi) parameter.
#' @details Formular of a block maxima distribution function:
#' \deqn{F_j(x)=F_j^{(b)}(x)=\left[2\cdot T_{1/\xi}\left(\left\{1+\xi\frac{x-\mu_j}{\sigma_j}\right\}\cdot T_{1/\xi}^{-1}\left(1-\frac{1}{2b}\right)\right)-1\right]^b,}{F_j(x)=F_j^(b)(x)=[2 T_{1/xi}({1+xi (x-mu_j)/(sigma_j)} T_{1/xi}^{-1}(1-1/(2b)))-1]^b,}
#' where \eqn{T_{\nu}}{T_nu} denote the t-distribution function with \eqn{\nu}{nu} degrees of freedom.
#' @return Quantile of a block maxima distribution.
#' @examples
#' qBM(p=c(0.75, 0.99), b=12, param=c(2,1,0.2))
#' @export
qBM <- function(p, b, param){
  if(!all(p < 1 & p > 0)) stop("p must be a probability: p in (0,1)")

  mu <- param[1]
  sigma <- param[2]
  xi <- param[3]

  return(mu + sigma/xi*( qt((1+p^(1/b))/2, df=1/xi)/qt((2-1/b)/2, df=1/xi) - 1))
}

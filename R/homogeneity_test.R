GANOVA <- function(est, Sig){
  d <- length(est)
  Sig_Inv <- solve(Sig)
  opt_weights <- 1/sum(Sig_Inv)*colSums(Sig_Inv)
  est_reg <- sum(opt_weights*est)
  stat <- t(est-est_reg)%*%Sig_Inv%*%(est-est_reg)
  pval <- 1-pchisq(stat, d-1)
  list(stat = stat, pval = pval)
}

#' Homogeneity test for the shape
#' @description A test for assumption \eqn{H_0}: "shape parameter is equal for all d GEV margins"
#' with test statistic based on (trimmed) L-moments.
#' @param x matrix of observations (rows: observations, d columns: stations).
#' @param leftrim	integer indicating lower trimming parameter (\eqn{\ge 0}).
#' @param rightrim	integer indicating upper trimming parameter (\eqn{\ge 0}).
#' @return p-value of the test.
#' @examples
#' library("evd")
#' # sample observations of 100 years at 5 stations:
#' set.seed(1053)
#' x19 <- matrix(rgev(400, 2, 1, 0.1), ncol=4) # 4 stations with the same shape
#' x10 <- rgev(100, 2, 1, 0.4) # one station with a different shape
#' x <- cbind(x19, x10)
#' xiAnova(x=x, leftrim=0, rightrim=1)
#' @importFrom stats pchisq qnorm qt runif uniroot
#' @importFrom TLMoments TLMoments parameters est_cov PWMs quantiles
#' @importFrom evd dgev pgev qgev
#' @importFrom copula An.biv
#' @importFrom magrittr %>%
#' @export
xiAnova <- function(x, leftrim=0, rightrim=1){
  d <- ncol(x)
  par_est <- parameters(TLMoments(x, leftrim=leftrim, rightrim = rightrim, na.rm=TRUE), "gev")
  xi_est <- par_est[3,]
  Sig_xi <- est_cov(par_est)[(1:d)*3, (1:d)*3]
  GANOVA(xi_est, Sig_xi)$pval
}

#' @title Estimated confidence intervals using sTL
#' @description Estimated regional (or local) (1-alpha)-confidence intervals for an estimated
#' quantile by using seasonal TL(0,1)-moments (sTL).
#' @param x1 vector or matrix of observations from season 1 (rows: observations, columns: stations).
#' @param x2 vector or matrix of observations from season 2 (rows: observations, columns: stations).
#' @param p a probability.
#' @param j quantile and parameter estimation for the jth station (jth column of \code{x}).
#' Irrelevant if is \code{x1} and \code{x2} are vectors.
#' @param alpha confidence level for confidence interval.
#' @return List of \itemize{
#' \item \code{ci} confidence interval.
#' \item \code{quant} estimated quantile from a two-component GEV using trimmed L-moments (leftrim=0, rightrim=1).
#' \item \code{var} variance of the estimated quantile.}
#' @examples
#' library("evd")
#' # Seasonal observations of 80 years at one station:
#' x1 <- rgev(80, 2, 1, 0.2) # observations from season 1
#' x2 <- rgev(80, 3, 1, 0.3) # observations from season 2
#' confInt_sTL(x1=x1, x2=x2, p=0.95, alpha=0.05)
#'
#' # Seasonal observations of 100 years at 4 stations:
#' x1 <- matrix(rgev(400, 2, 1, 0.3), ncol=4)
#' x2 <- matrix(rgev(400, 4, 1, 0.2), ncol=4)
#' confInt_sTL(x1=x1, x2=x2, j=2, p=0.95, alpha=0.05)
#' @export
confInt_sTL <- function(x1, x2, p, j=1, alpha=0.05){

  if(!all(p < 1 & p > 0)) stop("p must be a probability: p in (0,1)")

  ## local:
  if(is.vector(x1) && is.vector(x2)){
    n_som <- sum(!is.na(x1)) # number of summer observations
    n_win <- sum(!is.na(x2)) # number of winter observations
    n_ges <- max(n_som, n_win)

    para_som <- TLMoments(x1, leftrim=0, rightrim=1, na.rm=TRUE) %>% parameters("gev")
    Sigma_som <- n_som * (para_som %>% est_cov)

    para_win <- TLMoments(x2, leftrim=0, rightrim=1, na.rm=TRUE) %>% parameters("gev")
    Sigma_win <- n_win * (para_win %>% est_cov)

  }
  ## regional:
  if(is.matrix(x1) && is.matrix(x2)){
    n_som <- apply(x1, 2, function(y) sum(!is.na(y))) # number of summer observations
    n_win <- apply(x2, 2, function(y) sum(!is.na(y))) # number of winter observations
    n_ges <- max(n_som, n_win)

    nj_som <- n_som[j]
    nj_win <- n_win[j]
    d <- ncol(x1)

    Sigma_fun <- function(x, j, d, n){
      para <- TLMoments(x, leftrim=0, rightrim=1, na.rm=TRUE) %>% parameters("gev")

      Sigma <- n * (para %>% est_cov)
      Sigma_Inv <- solve(Sigma[1:d*3, 1:d*3])
      ev <- eigen(Sigma_Inv)$values
      if(any(ev < 0)) w <- n/sum(n) # easy weighting if eigenvalue < 0
      else w <- rowSums(Sigma_Inv)/sum(Sigma_Inv) # optimal weighting

      Trafo <- matrix(0, nrow=3, ncol=3*d)
      Trafo[1, 3*j-2] <- 1
      Trafo[2, 3*j-1] <- 1
      Trafo[3, 1:d*3] <- w

      Sigma_neu <- Trafo %*% Sigma %*% t(Trafo)
      paraj <- c(para[1:2,j], w %*% para[3,])

      return(list(paraj=paraj, Sigma=Sigma_neu))
    }

    sommer <- Sigma_fun(x1, j=j, d=d, n=nj_som)
    para_som <- sommer$paraj
    Sigma_som <- sommer$Sigma

    winter <- Sigma_fun(x2, j=j, d=d, n=nj_win)
    para_win <- winter$paraj
    Sigma_win <- winter$Sigma
  }

  # estimated quantile
  q_p <- qGEVxGEV(p, para_som, para_win)

  # density g (GEV):
  g_win <- dgev(q_p, para_win[1], para_win[2], para_win[3])
  g_som <- dgev(q_p, para_som[1], para_som[2], para_som[3])

  # distribution G (GEV):
  G_win <- pgev(q_p, para_win[1], para_win[2], para_win[3])
  G_som <- pgev(q_p, para_som[1], para_som[2], para_som[3])

  # Jacobi-matrix J (1x3) from G:
  G_mu <- function(x, mu, sig, xi) -1/sig * (1 + xi/sig * (x-mu))^(-1/xi - 1) * pgev(x, mu, sig, xi)

  G_sig <- function(x, mu, sig, xi) -(x-mu)/sig^2 * (1 + xi/sig * (x-mu))^(-1/xi - 1) * pgev(x, mu, sig, xi)

  G_xi <- function(x, mu, sig, xi){
    term1 <- -pgev(x, mu, sig, xi)
    term2 <- (1 + xi/sig * (x-mu))^(-1/xi)
    term3 <- log(1 + xi/sig * (x-mu))/xi^2  -  (x-mu)/(sig*xi * (1 + xi/sig * (x-mu)))
    erg <- prod(term1, term2, term3)
    return(erg)
  }

  J_win <- c(G_mu(q_p, para_win[1], para_win[2], para_win[3]),
             G_sig(q_p, para_win[1], para_win[2], para_win[3]),
             G_xi(q_p, para_win[1], para_win[2], para_win[3]))

  J_som <- c(G_mu(q_p, para_som[1], para_som[2], para_som[3]),
             G_sig(q_p, para_som[1], para_som[2], para_som[3]),
             G_xi(q_p, para_som[1], para_som[2], para_som[3]))


  varianz <- as.numeric((G_som^2 * J_win %*% Sigma_win %*% J_win + G_win^2 * J_som %*% Sigma_som %*% J_som)/((g_win * G_som + G_win*g_som)^2))

  z_alpha <- qnorm(1-alpha/2)

  ki <- q_p + c(-1, 1) * z_alpha * sqrt(varianz/n_ges)

  return(list(ci=ki, quant=q_p, var=varianz))
}
################################################################

#' @title Estimated confidence intervals using TL
#' @description Estimated regional (or local) (1-alpha)-confidence intervals for an estimated
#' quantile by using annual TL(0,1)-moments (TL).
#' @param x vector or matrix of annual observations.
#' @param p a probability.
#' @param j quantile and parameter estimation for the jth station (jth column of \code{x}).
#' Irrelevant if is \code{x} is a vector.
#' @param alpha confidence level for confidence interval.
#' @return List of \itemize{
#' \item \code{ci} confidence interval.
#' \item \code{quant} estimated quantile from a GEV using trimmed L-moments (\eqn{leftrim=0, rightrim=1}).}
#' @examples
#' library("evd")
#' # Seasonal observations of 80 years at one station:
#' x1 <- rgev(80, 2, 1, 0.2) # observations from season 1
#' x2 <- rgev(80, 3, 1, 0.3) # observations from season 2
#' x <- seas2ann(x1, x2) # calculaes annual maxima of the two seasons
#' confInt_TL(x=x, p=0.95, alpha=0.05)
#'
#' # Seasonal observations of 100 years at 4 stations:
#' x1 <- matrix(rgev(400, 2, 1, 0.3), ncol=4) # observations from season 1
#' x2 <- matrix(rgev(400, 2, 1, 0.2), ncol=4) # observations from season 2
#' x <- seas2ann(x1, x2) # calculaes annual maxima of the two seasons
#' confInt_TL(x=x, j=2, p=0.95, alpha=0.05)
#' @export
confInt_TL <- function(x, p, j=1, alpha=0.05){
  if(!all(p < 1 & p > 0)) stop("p must be a probability: p in (0,1)")

  if(is.vector(x)) x <- matrix(x, ncol=1) # local

  n <- apply(x, 2, function(y) sum(!is.na(y))) # number of observations
  n_ges <- max(n)
  nj <- n[j]
  d <- ncol(x)

  Sigma_fun <- function(x, j, d, n){
    para <- TLMoments(x, leftrim=0, rightrim=1, na.rm=TRUE) %>% parameters("gev")

    Sigma <- n * (para %>% est_cov)
    Sigma_Inv <- solve(Sigma[1:d*3, 1:d*3])
    ev <- eigen(Sigma_Inv)$values
    if(any(ev < 0)) w <- n/sum(n) # easy weighting, if eigenvalue < 0
    else w <- rowSums(Sigma_Inv)/sum(Sigma_Inv) # optimal weighting

    Trafo <- matrix(0, nrow=3, ncol=3*d)
    Trafo[1, 3*j-2] <- 1
    Trafo[2, 3*j-1] <- 1
    Trafo[3, 1:d*3] <- w

    Sigma_neu <- Trafo %*% Sigma %*% t(Trafo)
    paraj <- c(para[1:2,j], w %*% para[3,])

    return(list(paraj=paraj, Sigma=Sigma_neu))
  }

  jahr <- Sigma_fun(x, j=j, d=d, n=nj)
  para <- as.vector(jahr$paraj)
  Sigma <- jahr$Sigma


  # estimated quantile:
  q_p <- qgev(p, para[1], para[2], para[3])

  # Dg(theta):
  Ginv_mu <- 1

  Ginv_sig <- function(p, mu, sig, xi) 1/xi *((-log(p))^(-xi)-1)

  Ginv_xi <- function(x, mu, sig, xi){
    term1 <- -sig/xi^2 * ((-log(p))^(-xi)-1)
    term2 <- sig/xi * (-log(p))^(-xi) * log(-log(p))
    erg <- term1 - term2
    return(erg)
  }

  Dmat <- c(Ginv_mu,
            Ginv_sig(p, para[1], para[2], para[3]),
            Ginv_xi(p, para[1], para[2], para[3]))

  varianz <- as.numeric(t(Dmat) %*% Sigma %*% Dmat)

  z_alpha <- qnorm(1-alpha/2)

  ki <- q_p + c(-1, 1) * z_alpha * sqrt(varianz/n_ges)

  return(list(ci=ki, quant=q_p, var=varianz))
}

################################################################

#' @title Estimated confidence intervals using W
#' @description Estimated regional (or local) (1-alpha)-confidence intervals for an estimated
#' quantile by using the annual Weissman estimator (W).
#' @param x vector or matrix of annual observations.
#' @param p a probability. Should be between \eqn{1-k_j/n_j} and 1, where \eqn{n_j} is the
#' sample length of the \code{j}-th column.
#' @param j quantile and parameter estimation for the jth station (jth column of \code{x}).
#' Irrelevant if is \code{x1} and \code{x2} are vectors.
#' @param alpha confidence level for confidence interval.
#' @param ... additional arguments, see \link{RegioWeissman}
#' @return List of \itemize{
#' \item \code{ci} confidence interval.
#' \item \code{quant} estimated quantile using Weissman's estimator.}
#' @examples
#' library("evd")
#' # Seasonal observations of 80 years at one station:
#' x1 <- rgev(80, 2, 1, 0.2) # observations from season 1
#' x2 <- rgev(80, 3, 1, 0.3) # observations from season 2
#' x <- seas2ann(x1, x2) # calculaes annual maxima of the two seasons
#' confInt_W(x=x, p=0.95, alpha=0.05)
#'
#' # Seasonal observations of 100 years at 4 stations:
#' x1 <- matrix(rgev(400, 2, 1, 0.3), ncol=4) # observations from season 1
#' x2 <- matrix(rgev(400, 2, 1, 0.2), ncol=4) # observations from season 2
#' x <- seas2ann(x1, x2) # calculaes annual maxima of the two seasons
#' confInt_W(x=x, j=2, p=0.95, alpha=0.05)
#' @export
confInt_W <- function(x, p, j=1, alpha=0.05, ...){
  val <- RegioWeissman(x=x, j=j, p=p, alpha=alpha, ...)
  return(list(ci=val$CI, quant=val$est))
}




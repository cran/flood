################ Hill's estimator ######################
#' @title Hill's estimator
#' @description Estimation of heavy tails with Hill's estimator
#' @param x Vector or matrix of observations
#' @param k Number of relative excesses involved in the estimation of the extreme value
#' index gamma. If \code{k} is missing, it will be set to \itemize{
#'  \item \eqn{k=\left\lfloor 2n^{2/3}\right\rfloor}{k=floor(2*n^(2/3))}, where n is the
#' sample length of the vector \code{x} after removing missing values
#'  \item \eqn{k=\left\lfloor \frac{2n^{2/3}}{d^{1/3}}\right\rfloor}{k=floor(2*n^(2/3)/d^(1/3))}, where d is the number of
#'  columns of the matrix \code{x} and n the length of each column after removing missing values.}
#' @return Hill's estimator for each sample.
#' @examples
#' library("evd")
#' x1 <- rgev(100, loc = 2, scale = 1, shape=0.4)
#' hill(x1, k=20)
#' x2 <- rgev(100, loc = 2, scale = 1, shape=0.5)
#' hill(cbind(x1, x2), k = c(20, 25))
#' x2[c(4,8,39)] <- NA
#' hill(cbind(x1, x2), k=c(20, 25))
#' # if leaving out k, it will be set to floor(2*n^(2/3)/d^(1/3)) = c(34,33):
#' hill(cbind(x1, x2)) # is the same as:
#' hill(cbind(x1, x2), k=c(34,33))
#' @export
hill <- function(x, k){
  if(!missing(k)) k <- floor(k)

  if(is.vector(x)){
    x <- sort(x) # removes NAs
    n <- length(x)
    if(missing(k)) k <- floor(2*n^(2/3))
    H <- mean(log(x[n-1:k+1]/x[n-k]), na.rm=TRUE)
  }

  if(is.matrix(x)){
    x <- apply(x, 2, sort, na.last=FALSE)
    d <- ncol(x)
    n <- apply(x, 2, function(z) sum(!is.na(z)))
    if(!missing(k) && length(k)!=d) stop("wrong length of k compared to number of columns of x")
    if(missing(k)) k <- floor(2*n^(2/3)/d^(1/3))
    H <- numeric(d)
    for(i in 1:d){
      y <- x[!is.na(x[,i]),i]
      H[i] <- mean(log(y[n[i]-1:k[i]+1]/y[n[i]-k[i]]), na.rm=TRUE)
    }
  }
  return(H)
}

################ Regional EVI estimator ######################
#' @title Regional EVI estimator
#' @description Estimation of the positive extreme value index (EVI) based on multiple local Hill estimators. We assume heavy-tail homogeneity, i.e., all local EVI's are the same.
#' @param x Vector or matrix of observations
#' @param k Number of relative excesses involved in the estimation of the extreme value
#' index gamma. If \code{k} is missing, it will be set to \itemize{
#'  \item \eqn{k=\left\lfloor 2n^{2/3}\right\rfloor}{k=floor(2*n^(2/3))}, where n is the
#' sample length of the vector \code{x} after removing missing values
#'  \item \eqn{k=\left\lfloor \frac{2n^{2/3}}{d^{1/3}}\right\rfloor}{k=floor(2*n^(2/3)/d^(1/3))}, where d is the number of
#'  columns of the matrix \code{x} and n the length of each column after removing missing values.}
#' @param k.qu Tuning parameter for estimation of empirical variance; only needed if \code{type="opt"}.
#' @param type Choose either \code{"evopt"} if extreme value dependent, \code{"ind"} if independent or \code{"opt"} for arbitrarily dependent components.
#' @param alpha Confidence level for confidence interval.
#' @param ci Either \code{"nonlog"} for standart or \code{"log"} for non-standart confidence interval based on log-transformed hill estimates.
#' @return List of \itemize{
#'  \item \code{est} a weighted average of local Hill estimates.
#'  \item \code{Sigma} an estimate of the corresponding variance matrix.
#'  \item \code{CI} a confidence interval.}
#' @examples
#' library("evd")
#' x1 <- rgev(150, loc = 2, scale = 1, shape=0.4)
#' hill(x1, k=20)
#' x2 <- rgev(100, loc = 2.5, scale = 1, shape=0.4)
#' x2 <- c(x2, rep(NA, 50))
#' x <- cbind(x1, x2)
#' k <- c(40, 30)
#' RegioHill(x, k)
#' @export
RegioHill <- function(x, k, k.qu=20, type="evopt", alpha=0.05, ci="nonlog"){

  if(!missing(k)) k <- floor(k)

  ### purely local:
  if(is.vector(x)){
    if(missing(k)) k <- floor(2*length(x)^(2/3))
    ## Hill-estimator:
    HEst <- hill(x=x, k=k)

    ## confidence interval:
    if(ci=="nonlog") confint <- HEst+c(-1,1)*HEst*qnorm(1-alpha/2)*sqrt(1/k)
    if(ci=="log") confint <- HEst*exp(c(-1,1)*qnorm(1-alpha/2)*sqrt(1/k))

    return(list(est=HEst, CI=confint))
  }

  ### regional:
  if(is.matrix(x)){
    d <- ncol(x)
    if(missing(k)){
      n <- apply(x, 2, function(z) sum(!is.na(z)))
      k <- floor(2*n^(2/3)/d^(1/3))
    }

    ## Estimation of variance matrix Sigma and computation of weights:
    if(type=="ind"){
      Sigma <- diag(k[1]/k)
      w <- k/sum(k)
    }

    if(type=="evopt"){
      Sigma <- Sigma_hill_hat(x=x, k=k, k.qu=k.qu, ev=TRUE)
      SigmaInv <- solve(Sigma)
      w <- 1/sum(SigmaInv)*colSums(SigmaInv)
    }

    if(type=="opt"){
      Sigma <- Sigma_hill_hat(x=x, k=k, k.qu=k.qu, ev=FALSE)
      SigmaInv <- solve(Sigma)
      w <- 1/sum(SigmaInv)*colSums(SigmaInv)
    }

    ## Hill-estimator vector:
    HEst <- hill(x=x, k=k)
    gammaEst <- as.vector(w%*%HEst)

    ## confidence interval:
    if(ci=="nonlog") confint <- gammaEst+c(-1,1)*gammaEst*qnorm(1-alpha/2)*sqrt(t(w)%*%Sigma%*%w/k[1])
    if(ci=="log") confint <- gammaEst*exp(c(-1,1)*qnorm(1-alpha/2)*sqrt(t(w)%*%Sigma%*%w/k[1]))

    return(list(est=gammaEst, Sigma=Sigma, w=w, CI=confint))
  }
}

### Auxilliary function: Estimation of limiting variance matrix
Sigma_hill_hat <- function(x, k, k.qu=20, ev=TRUE){
  if(is.vector(x)) M <- 1

  if(is.matrix(x)){
    d <- ncol(x)
    n <- apply(x, 2, function(z) sum(!is.na(z)))

    if(!missing(k) && length(k)!=ncol(x)) stop("Invalid length of k")
    if(!missing(k)) k <- floor(k)
    if(missing(k)) k <- floor(2*n^(2/3)/d^(1/3))

    ck <- k[1]/k
    tau <- n/max(n)

    M <- matrix(NA, ncol=d, nrow=d)
    diag(M) <- ck

    for(l in 1:d){
      for(m in (1:d)[-which(1:d==l)]){
        tauMin <- min(tau[m], tau[l])
        vorfak <- ck[l]*ck[m]*tauMin

        ind <- which(!is.na(x[,l]+x[,m]))
        Dlm <- cbind(x[ind,l], x[ind,m])
        Nlm <- nrow(Dlm)
        arg <- c(1/(tau[l]*ck[l]), 1/(tau[m]*ck[m]))

        if(ev) Lambda <- sum(arg) * (1-min(1, copula::An.biv(x=Dlm, w=arg[2]/(arg[1]+arg[2]), "CFG")))
        else{
          DlmSort <- apply(Dlm, 2, sort)
          cut <- c(ceiling(Nlm-k.qu*arg[1]), ceiling(Nlm-k.qu*arg[2]))
          Lambda <- 1/k.qu * sum((Dlm[,1]>DlmSort[cut[1],1]) * (Dlm[,2]>DlmSort[cut[2],2]))
        }

        M[l,m] <- vorfak * Lambda
      }
    }
  }

  return(M)
}



################ Test of heavy-tail homogeneity ######################
#' @title Heavy-tail ANOVA
#' @description A test of heavy-tail homogeneity, that is, equality of the positive extreme value index for all d columns of \code{x}.
#' @param x Matrix of observations
#' @param k Number of relative excesses involved in the estimation of the extreme value
#' index gamma. If \code{k} is missing, it will be set to
#' \eqn{k=\left\lfloor \frac{2n^{2/3}}{d^{1/3}}\right\rfloor}{k=floor(2*n^(2/3)/d^(1/3))}, where d is the number of
#'  columns of the matrix \code{x} and n the length of each column after removing missing values.
#' @param k.qu Tuning parameter for estimation of empirical variance; only needed if \code{type="opt"}.
#' @param type Choose either \code{"evopt"} if extreme value dependent, \code{"ind"} if independent or \code{"opt"} for arbitrarily dependent components.
#' @param cf If \code{TRUE}, a correctur factor is used, which improves the size at the cost of power.
#' @return Test statistic and p-value.
#' @examples
#' library("evd")
#' set.seed(6754)
#' x1 <- rgev(150, loc = 2, scale = 1, shape=0.4)
#' x2 <- rgev(150, loc = 2.5, scale = 1, shape=0.1) # H_0 violated because of different shapes
#' x <- cbind(x1, x2)
#' TailAnova(x)
#'
#' x1 <- rgev(150, loc = 2, scale = 1, shape=0.3)
#' x2 <- rgev(150, loc = 2.5, scale = 1, shape=0.3) # H_0 not violated because of same shapes
#' x <- cbind(x1, x2)
#' TailAnova(x)
#' @export
TailAnova <- function(x, k, k.qu=20, type="evopt", cf=TRUE){
  d <- ncol(x)
  n <- apply(x, 2, function(z) sum(!is.na(z)))

  if(missing(k)) k <- floor(2*n^(2/3)/d^(1/3))
  if(!missing(k)) k <- floor(k)

  HEst <- hill(x=x, k=k)
  weight.erg <- RegioHill(x=x, k=k, k.qu=k.qu, type=type, ci="nonlog")
  w <- weight.erg$w
  Sigma <- weight.erg$Sigma
  gammaEst <- weight.erg$est
  teststat <- as.vector(k[1]/gammaEst^2 * t(HEst-gammaEst) %*% solve(Sigma) %*% (HEst-gammaEst))
  if(cf) teststat <- teststat * (1-d/(5*min(n)))
  p.val <- 1-pchisq(teststat, d-1)

  return(list(stat=teststat, pval=p.val))
}


####################### Regional Weissman #######################
#' @title Quantile estimation: Weissman's extrapolation
#' @description Estimation of the p-quantile based on multiple local Hill estimators and Weissman's extrapolation formula. We assume heavy-tail homogeneity, i.e., all local EVI's are the same.
#' @param x Vector or matrix of observations
#' @param j The number of the target site, i.e., if \code{j=2} the p-quantile of the second column of \code{x} is estimated.
#' @param p The probability of interest; should be between \eqn{1-k_j/n_j} and 1, where \eqn{n_j} is the sample length of the j-th column.
#' @param k Number of relative excesses involved in the estimation of the extreme value
#' index gamma. If \code{k} is missing, it will be set to \itemize{
#'  \item \eqn{k=\left\lfloor 2n^{2/3}\right\rfloor}{k=floor(2*n^(2/3))}, where n is the
#' sample length of the vector \code{x} after removing missing values
#'  \item \eqn{k=\left\lfloor \frac{2n^{2/3}}{d^{1/3}}\right\rfloor}{k=floor(2*n^(2/3)/d^(1/3))}, where d is the number of
#'  columns of the matrix \code{x} and n the length of each column after removing missing values.}
#' @param k.qu Tuning parameter for estimation of empirical variance; only needed if \code{type="opt"}.
#' @param type Choose either \code{"evopt"} if extreme value dependent, \code{"ind"} if independent or \code{"opt"} for arbitrarily dependent components.
#' @param alpha Confidence level for confidence interval.
#' @return List of \itemize{
#'  \item \code{est} Point estimate of p-quantile of column j
#'  \item \code{CI} the corresponding alpha-confidence interval
#'  \item \code{EVI} estimate of the extreme value index
#'  \item \code{k} tail sample size
#'  \item \code{p} a probability
#'  \item \code{u.kn} (n-k)-th largest observation, where n is the sample length at station j after removing missing values
#'  \item \code{description} a short description.}
#' @examples
#' library("evd")
#' # sample observations of 75 years at one station:
#' x <- rgev(75, 0, 1, 0) # x is a vector
#' RegioWeissman(x=x, p=0.95)
#'
#' x2 <- c(NA, NA, x[1:60], NA, x[61:75]) # vector of observations with missing values
#' RegioWeissman(x=x2, p=0.95) # NAs will be removed
#'
#' # sample observations of 100 years at 4 stations:
#' set.seed(1053)
#' x <- matrix(rgev(400, 2, 1, 0.3), ncol=4)
#' RegioWeissman(x=x, p=0.9, j=3)
#'
#' # With missing values:
#' x2 <- x
#' x2[sample(100, 12),2] <- NA
#' RegioWeissman(x=x2, p=0.9, j=3)
#' @export
RegioWeissman <- function(x, j=1, p, k, k.qu=20, type="evopt", alpha=0.05){
  ## case d=1 (local):
  if(is.vector(x)){
    x <- sort(x)
    n <- length(x)
    if(missing(k)) k <- floor(2*n^(2/3))
    if(!missing(k)) k <- floor(k)
    if(missing(p)) p <- 1-1/n
    if(p < 1-k/n || p >= 1) stop("Error: p not in [1-k/n, 1)")

    gammaEst <- hill(x=x, k=k)
    u.kn <- x[n-k]
    qu <- u.kn *(k/(n*(1-p)))^gammaEst
    conf.int <- qu * (1 + c(-1,1) * qnorm(1-alpha/2)* sqrt(gammaEst^2/k) * log(k/(n*(1-p))))

    return(list(est=qu, CI=conf.int, EVI=gammaEst, k=k, p=p, u.kn=u.kn, description="Quantile estimation method: Local Weissman extrapolation"))
  }

  ## case d>1 (regional):
  if(is.matrix(x)){
    n <- apply(x, 2, function(z) sum(!is.na(z)))
    d <- ncol(x)
    if(missing(k)) k <- floor(2*n^(2/3)/d^(1/3))
    if(!missing(k)) k <- floor(k)
    if(missing(p)) p <- 1-1/n[j]
    if(p < 1-k[j]/n[j] || p >= 1) stop('Error: p not in [1-k[j]/n[j], 1)')

    ind.j <- which(is.na(x[,j]))
    if(length(ind.j > 0)) x.j <- sort(x[-ind.j,j])
    else x.j <- sort(x[,j])

    w1 <- RegioHill(x=x, k=k, k.qu=k.qu, type=type) # joint Hill-estimator is w1$est
    u.kn <- x.j[n[j]-k[j]]

    qu <- u.kn *(k[j]/(n[j]*(1-p)))^w1$est
    conf.int <- qu * (1 + c(-1,1) * qnorm(1-alpha/2) * sqrt((w1$est^2 * w1$w %*% w1$Sigma %*% w1$w)/k[1]) * log(k[j]/(n[j]*(1-p))))

    return(list(est=qu, CI=conf.int, EVI=w1$est, k=k, p=p, u.kn=u.kn, description="Quantile estimation method: Regional Weissman extrapolation"))
  }
}

#### Seasonal (regional or local) Weissman
#' @title Quantile estimation: Weissman's extrapolation for seasonal data
#' @description Estimation of the p-quantile based on multiple local Hill estimators and Weissman's extrapolation
#' formula with seasonality. We assume heavy-tail homogeneity within each season, i.e., all local EVI's are the same.
#' @param x List of 2 elements: each element is consistent a vector or matrix of observations
#' @param j The number of the target site, i.e., if \code{j=2} the p-quantile of the second station is estimated.
#' @param p The probability of interest; should be between \eqn{1-k_j/n_j} and 1, where \eqn{n_j} is the
#' sample length of the \code{j}-th column.
#' @param ... additional arguments, see \link{RegioWeissman}
#' @return Point estimate of seasonal p-quantile of column \code{j} and a short description.
#' @examples
#' library("evd")
#' # Local & seasonal (observations of 80 years at one station):
#' x1 <- rgev(80, 2, 1, 0.2) # observations from season 1
#' x2 <- rgev(80, 3, 1, 0.3) # observations from season 2
#' x <- list(x1, x2)
#' RegioWeissmanSeas(x=x, j=1, p=0.99)
#'
#' x1 <- matrix(rgev(400, 3, 1, 0.1), ncol=4) # observations of season 1 at 4 stations
#' x2 <- matrix(rgev(400, 2, 1, 0.3), ncol=4) # observations of season 2 at 4 stations
#' x <- list(x1, x2)
#' RegioWeissmanSeas(x=x, j=1, p=0.99)
#' @export
RegioWeissmanSeas <- function(x, j=1, p, ...){
  if(is.vector(x[[1]])){
    x <- lapply(x, matrix, ncol=1)
    descr <- "Quantile estimation method: Local+Seasonal Weissman extrapolation"
  }
  else descr <- "Quantile estimation method: Regional+Seasonal Weissman extrapolation"

  est <- lapply(x, RegioWeissman, j=j, p=p, ...)

  s1 <- est[[1]] # Regional/Local Weissman estimation for season 1
  s2 <- est[[2]] # Regional/Local Weissman estimation for season 2

  # k, n and u for each season at the jth station
  n1 <- sum(!is.na(x[[1]][,j]))
  n2 <- sum(!is.na(x[[2]][,j]))
  k1 <- s1$k[j]
  k2 <- s2$k[j]

  u.kn1 <- s1$u.kn
  u.kn2 <- s2$u.kn

  u.m <- max(c(u.kn1, u.kn2))

  Ffun <- function(y){
    (1 - k1/n1 * (y/u.kn1)^(-1/s1$EVI)) * (1 - k2/n2 * (y/u.kn2)^(-1/s2$EVI))
  }

  Q <- uniroot(function(y) Ffun(y) - p, interval=c(u.m, 1e100))$root

  return(list(est=Q, description=descr))
}



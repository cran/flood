#' Annual maxima from seasonal maxima
#' @description Calculates annual maxima from seasonal maxima of two seasons.
#' @param x1 vector or matrix of observations from season 1 (rows: observations, columns: stations).
#' @param x2 vector or matrix of observations from season 2 (rows: observations, columns: stations).
#' @param na.rm logical. \code{TRUE}: the annual maximum will be calculated even if one
#' observation of the two seasons is a missing value, e.g. winter maximum is 58 and summer
#' maximum is \code{NA} the annual maximum is 58. If both observations are NA, the annual maximum is set to \code{NA}, too.
#' \code{FALSE}: the annual maximum will be set to \code{NA} if one
#' observation of the two seasons is a missing value, e.g. winter maximum is 58 and summer maximum
#' is \code{NA} the annual maximum is \code{NA}.
#' @return Matrix of annual observations (rows: observations, columns: stations).
#' @examples
#' set.seed(28379)
#' x1 <- matrix(round(rnorm(8, 20, 25)), ncol=2)
#' x1[2] <- NA
#' x2 <- matrix(round(rnorm(8, 20, 25)), ncol=2)
#' x2[c(2,5,6)] <- NA
#' x1
#' x2
#' seas2ann(x1,x2,TRUE)
#' seas2ann(x1,x2,FALSE)
#' @export
seas2ann <- function(x1, x2, na.rm=TRUE){
  if(is.vector(x1) && is.vector(x2)){
    x1 <- matrix(x1, ncol=1)
    x2 <- matrix(x2, ncol=1)
  }

  x_array <- array(0, dim=c(nrow(x1), ncol(x1),2))
  x_array[,,1] <- x1
  x_array[,,2] <- x2

  x_ann <- suppressWarnings(apply(x_array, 1:2, max, na.rm=na.rm))
  if(na.rm==TRUE && any(x_ann==-Inf)) x_ann[which(x_ann==-Inf)] <- NA

  return(x_ann)
}

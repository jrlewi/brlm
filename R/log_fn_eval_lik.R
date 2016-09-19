#' Full likelihood under Gaussian model
#'
#' @param ystst a data vector such that T(ystst)=T(y_obs)
#' @param X design matrix
#' @param Beta regression coeff
#' @param sig sigma, standard deviation of residuals
#' @return the log likelihood
#' @details designed for use within \code{fn.one.rep.y}
log.fn.eval.lik <- function(ystst, X,Beta, sig)
{
  lik <- sum(dnorm(ystst,mean=X%*%Beta,sd=sig, log=TRUE))
  return(lik)
}




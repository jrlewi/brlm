#' Full likelihood under Gaussian model or T model
#'
#' @param ystst a data vector such that T(ystst)=T(y_obs)
#' @param X design matrix
#' @param Beta regression coeff
#' @param sig sigma, standard deviation of residuals
#' @param base_model normal or t
#' @param nu degrees of freedom if base_model = "t"
#' @return the log likelihood
#' @details designed for use within \code{fn.one.rep.y}
log.fn.eval.lik <- function(ystst, X,Beta, sig, base_model = "normal", nu = NULL)
{
  if(base_model == "normal"){
  lik <- sum(dnorm(ystst,mean=X%*%Beta,sd=sig, log=TRUE))
  }
  if(base_model == "t"){
    lik <- sum(
      dt((ystst - X%*%Beta)/sig, df = nu, log = TRUE) - log(sig)
      )
  }
  return(lik)
}




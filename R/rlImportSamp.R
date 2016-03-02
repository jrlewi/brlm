#' Fitting restricted likelihood model using importance sampling technique.

#'
#' Full model: \deqn{
#' \beta~N(\mu_0, \Sigma_0),
#' \sigma^2~IG(\alpha, \beta),
#' y~N(X\beta, \sigma^2)
#' }
#' For the restricted likelihood, conditioning is done on a pair of location and scale statistics \eqn{T(y)=(b(y), s(y))}.
#' @inheritParams rlDirectEval
#' @param x	a matrix or data frame containing the explanatory variables. The matrix should include a vector of 1's if intercept is desired.
#' @param instDist Placeholder to allow for more user specific instrumental distributions in the future. Currently not used.
#' @param sdInstDist Vector of length 2 defining the scale for the instrumental distribution when instDist=NULL. In this case, the instrumental distribution for \eqn{\mu} and \eqn{\sigma^2} are independent normal and log normal distributions, respectively. The first value is then the standard deviation of the normal. The second is the standard deviation of \eqn{\log\sigma^2} (i.e. the \code{sdlog} argument in \code{\link[stats]{Lognormal}}). If left as NULL a mutltiple of 5 is used on the asymptotic variance covariance matrix. Importance sample weights should be examined to evaluate the appropriatness of this choice.
#' @param Nins number of samples from the instrumental distribution
#'
rlImportSamp<-function(x,y, psi,..., scale.est='Huber',k2=1.345, mu0, Sigma0, alpha, beta,instDist=NULL, sdInstDist =NULL, smooth=1,N,Nins, maxit=1000){

  if (!is.function(psi)){
    psi <- get(psi, mode = "function")
  }

  arguments<-list(...)
  if (length(arguments)) {
    pm <- pmatch(names(arguments), names(formals(psi)), nomatch = 0L)
    if (any(pm == 0L)){
      warning("some of ... do not match")
    }
    pm <- names(arguments)[pm > 0L]
    formals(psi)[pm] <- unlist(arguments[pm])
  }

  n<-length(y)
  # find summary statistics #----
  rlm.obs<-rlm(x,y,psi=psi, scale.est=scale.est, k2=k2,maxit=maxit)
  t.obs<-rlm.obs$coefficients
  s.obs<-rlm.obs$s
  log.s.obs<-log(s.obs)

#look at belguim calls data to finish this one off ----


}

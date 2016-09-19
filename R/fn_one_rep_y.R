#' Perform a single Metropolis-Hastings step
#'
#' Performs a single Metropolis-Hastings step for the restricted likelihood method
#'@param y.curr current data vector
#'@param Beta vector of regression coefficient parameters
#'@param sig sigma-the standard deviation parameter
#'@param l1obs conditioning statistic-estimate of Beta from observed data
#'@param s1obs conditioning statistic-estimate of sigma from observed data
#'@param X regression design matrix (i)
#'@param log.prop.den.curr value of the proposal density for the current data vector (free of other parameters so it can be saved from iteration to iteration to save computation time)
#'@param proj the projection matrix onto the deviation space (aka the least squares residual space)
#'@param fn.psi custom psi function returning psi(x) or psi'(x) defining M-estimator for location ; default is fn.psi.huber
#'@param fn.chi custom chi function returning chi(x) or chi'(x) defining M-estimator for scale; default (and currently only capability) is fn.chi.prop2
#'@param psi psi function returning psi(x)/x. Must match fn.spi
#'@param scaleEst 'Huber' (default and currently only capability).
#'@param maxit max iterations for the estimation of the summary statistics
#' @return list y.curr: the new (if accepted) or same (if rejected) data vector, a: indicator of acceptance of proposed value 1 if accepted, 0 if rejected, log.prop.den.curr: proposal density on log scale of the returned data vector, rat: MH ratio, log.rat1: log ratio of proposal densities

fn.one.rep.y<-function(y.curr,Beta,sig,l1obs, s1obs,X, log.prop.den.curr, proj,fn.psi,fn.chi,psi, scaleEst,maxit=400)
{
  n<-length(y.curr)
  p<-ncol(X)
  # propose new y vector
  y.prop <- rnorm(n)
  y.prop <- fn.comp.ystst(y.prop,X,l1obs, s1obs, psi,scaleEst,maxit)
  # compute acceptance ratio in steps
  log.rat1 <- log.fn.eval.lik(y.prop,X,Beta, sig)-log.fn.eval.lik(y.curr,X,Beta, sig)

  log.prop.den.prop<-log.prop.den(y.prop,X, proj, l1obs, s1obs,fn.psi,fn.chi, n,p)

  log.rat2 <- log.prop.den.curr - log.prop.den.prop
  ratPre<-exp(log.rat1+log.rat2)
  rat <- min(c(1,ratPre))
  tmp <- runif(1)
  if (tmp < rat) {y.curr <- y.prop; a<-1;log.prop.den.curr<-log.prop.den.prop} else {
    a<-0
  }
  return(list(y.curr,a,log.prop.den.curr,rat, log.rat1))
}

#' @rdname fn.one.rep.y
#' @inheritParams fn.attenuation2
#'@details Two equivalent versions, the second uses \code{fn.attenuation2} and is faster. This function is designed for use within an MCMC sampler. Many of the inputs are iteration dependent. The input \code{log.prop.den.curr} is designed as a time saver. Since this is free of other parameters it can be saved from iteration to iteration to save computation time. Note the use of both fn.psi (custome in the \link{brlm} package) and psi (from the \link{MASS} package). This is poor coding and should be corrected.
fn.one.rep.y2<-function(y.curr,Beta,sig,l1obs, s1obs,X, log.prop.den.curr, proj, Qt,fn.psi,fn.chi,psi, scaleEst,maxit=400)
{
  n<-length(y.curr)
  p<-ncol(X)
  # propose new y vector
  y.prop <- rnorm(n)
  y.prop <- fn.comp.ystst(y.prop,X,l1obs, s1obs, psi,scaleEst,maxit)
  # compute acceptance ratio in steps
  log.rat1 <- log.fn.eval.lik(y.prop,X,Beta, sig)-log.fn.eval.lik(y.curr,X,Beta, sig)

  #log.prop.den.curr <-log(fn.attenuation(y.curr,X, proj=proj, l1obs,s1obs))-(n-p-1)*log(fn.radius(y.curr, X))
  #log.prop.den.prop<-log.prop.den(y.prop,X, proj, l1obs, s1obs,fn.psi,fn.chi, n,p)
  log.prop.den.prop<-log.prop.den2(y.prop,X, proj,Qt, l1obs, s1obs,fn.psi,fn.chi, n,p)

  log.rat2 <- log.prop.den.curr - log.prop.den.prop
  ratPre<-exp(log.rat1+log.rat2)
  #print(rat)
  rat <- min(c(1,ratPre))
  tmp <- runif(1)
  if (tmp < rat) {y.curr <- y.prop; a<-1;log.prop.den.curr<-log.prop.den.prop} else {
    a<-0
  }
  #  update if needed
  return(list(y.curr,a,log.prop.den.curr,rat, log.rat1))
}




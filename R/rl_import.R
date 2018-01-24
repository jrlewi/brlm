#' Fitting restricted likelihood model using importance sampling technique.

#'
#' Full model: \deqn{
#' \beta~N(\mu_0, \Sigma_0),
#' \sigma^2~IG(\alpha, \beta),
#' y~N(X\beta, \sigma^2)
#' }
#' For the restricted likelihood, conditioning is done on a pair of location and scale statistics \eqn{T(y)=(b(y), s(y))}. Current implementation allows for these to be a pair of M-estimators as implemented in \code{\link[MASS]{rlm}}
#' \code{statistic} is a function of \code{X}, \code{y}  and outputing a p+1-dimensional vector with the p location statistics (for \code{beta}) and 1 scale statistic to be conditioned on.
#'
#'
#' @inheritParams rl_direct
#' @param X	a matrix or data frame containing the explanatory variables. The matrix should include a vector of 1's if intercept is desired.
#' @param statistic , character name of a function computing the location and scale statistic. See details for specification of this function.
#' @param scale vector length p+1 defining the standard deviations for the normal instrumental distributions for \eqn{\beta} and \eqn{\sigma^2} (with \eqn{\sigma^2} truncated at 0). Importance sample weights should be examined to evaluate the appropriatness of this choice.
#' @param Nins number of samples from the instrumental distribution
#' @return list with elements \code{impSamps}, \code{w}, \code{fit}. These are the \code{Nins} by \code{length(beta)})+1} matrix of importance samples, the corresponding weights, and the observed statistic.
#' @author John R. Lewis \email{lewis.865@@osu.edu}
#' @export
rl_import<-function(y,X, statistic , mu0, Sigma0, alpha, beta,scale, smooth=1,N,Nins){


  if (!is.function(statistic)){
    statistic <- get(statistic, mode = "function")
  }
  n <- length(y)
  p <- ncol(X)
  if(length(scale) != p+1){stop('length of scale not p+1')}

  # find summary statistics ----
  obs_stat<-statistic(y,X)
  if(length(obs_stat) != p+1){stop('statistic value not length p+1')}
  t.obs<-obs_stat[1:p]
  s.obs<-obs_stat[p+1]
  log.s.obs<-log(s.obs)
  meanInstBeta <-b.obs
  sigmaInstBeta <- diag(scale[1:p])
  meanInstSigma2 <- s.obs
  sdInstSigma2 <- scale[p+1]


  # rlm.estimators<-function(Y){
  #   rlm.out<-rlm(X,Y,psi=psi, scale.est=scale.est, k2=k2,maxit=maxit)
  #   return(c(as.numeric(rlm.out$coefficients), log(as.numeric(rlm.out$s))))
  # }
  #computing the estimators----
  estimators<-function(Y){
    out <- statistic(X,Y)
    c(out[1:p], log(out[p+1]))
  }



  Y <- matrix(rnorm(N*n,0,1),nrow=N,ncol=n)
  rlm.matrix <- apply(Y, MARGIN = 1, FUN=estimators)
  h <- apply(rlm.matrix,1, ks::hpi, binned=TRUE)
  H <- diag(smooth*h)

  kernel_density<-function(x1,x2){
    mean(dmvnorm(t(rlm.matrix),mean=c(x1,x2),sigma=H^2))
  }

  fn.posterior.estimate<-function(parms){
    l <- length(parms)
    betas <- parms[1:(l-1)]
    sigma2 <- parms[l]
    kDensityEst <- log(kernel_density((b.obs-betas)/sqrt(sigma2),log.s.obs-.5*log(sigma2))/(sigma2^(.5*(l-1)))) + dmvnorm(betas,mu0,Sigma0,log=TRUE) + log(dinvgamma(sigma2,alpha,beta))
    return(kDensityEst)
  }


  fn.logInsLikelihood <- function(parms){
    l <- length(parms)
    betas <- parms[1:(l-1)]
    sigma2 <- parms[l]
    mvtnorm::dmvnorm(betas, mean=meanInstBeta, sigma=sigmaInstBeta, log=TRUE)+dtnorm(sigma2,mean=meanInstSigma2, sd=sdInstSigma2,log=TRUE, lower=0)
  }

  fn.sample.instBeta <- function(){
    mvtnorm::rmvnorm(Nins, mean=meanInstBeta, sigma=sigmaInstBeta)
  }

  fn.sample.instSigma2<-function(){
    rtnorm(Nins,mean=meanInstSigma2, sd=sdInstSigma2, lower=0)
  }

  insBetas<-fn.sample.instBeta()
  insSigma2<-fn.sample.instSigma2()


  insS<-cbind(insBetas,insSigma2)
  logInsLikelihood<-apply(insS,1,FUN=fn.logInsLikelihood)

  logPost<-apply(insS,1,fn.posterior.estimate)
  wgts<-exp(logPost-logInsLikelihood)
  wgts<-wgts/sum(wgts)

  out<-list(impSamps = insS, w = wgts, fit = obs_stat)
  out
}


#' MCMC algorithm for Bayesian Linear Model
#'
#' Fits the linear regression model
#' \eqn{\beta~N(\mu_0, \Sigma_0)},
#' \eqn{\sigma^2~IG(a0, b0)},
#' \eqn{y~N(X\beta, \sigma^2 I)}
#'
#' Uses Gibbs sampling to sample from the posterior under the above linear regression model.
#'
#' @param  y vector of repsonse
#' @param X design matrix for regression. If intercept is desired, a column of 1's is needed
#' @param mu0 prior mean for beta
#' @param Sigma0 prior var-cov matrix for \eqn{\beta}
#' @param a0 prior shape parameter for \eqn{\sigma^2}
#' @param a0,b0 prior scale parameter for \eqn{\sigma^2}
#' @param sigma2Init initial value for \eqn{\sigma^2} in MCMC
#' @param nkeep number of iterations to keep
#' @param nburn number of iterations to toss
#' @return list with mcmc sample and mean fitted values
#' @export
bayesLm<-function(y, X,mu0, Sigma0, a0, b0,sigma2Int, nkeep=1e4, nburn=1e3){
  p<-ncol(X)
  n<-length(y)
  total<-nkeep+nburn
  betaSamples<-array(NA, dim=c(total,p))
  sigma2Samples<-numeric(total)
  #set inital values to current values
  sigma2Cur<-sigma2Int
  Sigma0Inv<-solve(Sigma0)

  #################################
  #Sampling Functions for Gibbs Sampler
  #################################
  #[sigma^2|---]=IG(a0+n/2, b0+.5*t(Y-Xbeta)(Y-Xbeta))
  #beta|----]=N(mun, Lambdan)
  #mun=Lambdan%*%(t(X)%*%Y+Sigma0^-1mu0)/sigma^2
  #Lambdan=sigma^2(t(X%*%X)+Sigma0^-1)^-1; (t(X%*%X)+Sigma0^-1)^-1=varCov_0

  #a0 b0 must be specified
  sampleSigma2<-function()
  {
    an<-a0+.5*n
    resid<-y-X%*%betaCur
    bn<-as.numeric(b0+.5*crossprod(resid,resid))
    return(rinvgamma(1,an,bn))
  }


  sampleBeta<-function(){
    varCov.n<-solve(t(X)%*%X/sigma2Cur+Sigma0Inv)
    mu.n<-varCov.n%*%(t(X)%*%y/sigma2Cur+Sigma0Inv%*%mu0)
    return(mvrnorm(1,mu.n,varCov.n))
  }
  ################################################
  #Start the gibbs sampler
  ################################################
  for(i in 1:total){
    #step one: [beta|sigma2,y]~N()
    betaCur<-sampleBeta()
    #step 2 [sigma2|---]~IG
    sigma2Cur<-sampleSigma2()
    betaSamples[i,]<-betaCur
    sigma2Samples[i]<-sigma2Cur
  }

  colnames(betaSamples)<-sapply(seq(1:p), FUN=function(x) paste('beta',x,sep=''))
  mcmcBetaSigma2<-mcmc(cbind(betaSamples[(nburn+1):total,],sigma2Samples[(nburn+1):total]))
  colnames(mcmcBetaSigma2)[p+1]<-'sigma2'
  fittedValues<-X%*%summary(mcmcBetaSigma2)$statistics[1:p,1]
  out<-list()
  out$mcmc<-mcmcBetaSigma2
  out$fitted.values<-fittedValues
  out
}

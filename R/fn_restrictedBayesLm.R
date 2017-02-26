#' MCMC algorithm for Restricted Likelihood Bayesian Linear Model
#'
#' Fits the linear regression model
#' \eqn{\beta~N(\mu_0, \Sigma_0)},
#' \eqn{\sigma^2~IG(a0, b0)},
#' \eqn{y~N(X\beta, \sigma^2 I)}
#'
#' Uses Gibbs sampling to sample from the restricted posterior under the above linear regression model conditioned on estimators \code{regEst} and \code{scaleEst} for \eqn{\beta} and \eqn{\sigma}, respectively. This is just the Gibbs sampler for the full posterior (implemented in \code{brlm::bayesLm} augmented with a step to sample new data given all of the parameters and the observed conditioning statistics
#'
#'
#' @inheritParams bayesLm
#' @param regEst The estimator of \eqn{\beta}  to condition on. Current options available are \code{Huber} and \code{Tukey} for Huber's and Tukey's (bisquare) estimators, respectively.
#' @param scaleEst The estimator of \eqn{\sigma} to condition on. Currently Huber's proposal 2 is the only option available. Specified with \code{Huber}.
#' @param maxit maximum number of iterations for use in \code{\link[MASS]{rlm}} when computing conditioning statistic during each proposal within the MCMC step for [y | rest, observed statistics]
#' @return list with mcmc sample, mean fitted values, robust regression object from \code{rlm}, and observed conditioning statistics.
#' @export
restrictedBayesLm<-function(y, X,regEst='Huber',scaleEst='Huber',mu0, Sigma0, a0, b0,sigma2Int, nkeep=1e4, nburn=1e3, maxit=400){
  #y is the response
  #X is the design Matrix
  # sigma2Int is the initial value for sigma2
  #mu0 prior mean of beta
  #Sigma0 is the var cov matrix of beta b~N(mu0,Sigma0)
  #a0, b0 prior parameters for sigma2
  #regEst: regression estimates used: options are 'Huber' and 'Tukey'
  #scaleEst='Huber' only option is huber's proposal two
  p<-ncol(X)
  n<-length(y)
  total<-nkeep+nburn
  Sigma0Inv<-solve(Sigma0)
  betaSamples<-array(NA, dim=c(total,p))
  sigma2Samples<-numeric(total)
  yAccept<-numeric(total)
  #set inital values to current values
  sigma2Cur<-sigma2Int

  Q<-qr.Q(qr(X))
  projMatrix<-diag(n)-tcrossprod(Q,Q) #Q%*%t(Q)

  ############################
  #define the psi and chi functions
  ############################
  if(regEst=='Huber') {
    psi<-get('psi.huber') #internal
    fn.psi<-get('fn.psi.huber')

  } else {
    if(regEst=='Tukey'){
      psi<-get('psi.bisquare') #internal
      fn.psi<-get('fn.psi.bisquare')
    } else {stop("only set up for Huber or Tukey regression estimates")}}

  if(scaleEst!='Huber'){
    stop('scale estimate only set up for Hubers Prop2 ')
  }
  fn.chi<-fn.chi.prop2
  ################################


  #run the robust regression
  robust<-rlm(X,y,psi=psi, scale.est=scaleEst, maxit=maxit)
  #condition on these estimates
  l1obs<-robust$coefficients
  s1obs<-robust$s

  #################################
  #Sampling Functions for Gibbs Sampler.
  #################################
  #[sigma^2|---]=IG(a0+n/2, b0+.5*t(Y-Xbeta)(Y-Xbeta))
  #[beta|----]=N(mun, Lambdan)
  #mun=Lambdan%*%(t(X)%*%Y+Sigma0^-1mu0)/sigma^2
  #Lambdan=sigma^2(t(X%*%X)+Sigma0^-1)^-1; (t(X%*%X)+Sigma0^-1)^-1=varCov_0

  #a0 b0 must be specified
  sampleSigma2<-function(data)
  {
    an<-a0+.5*n
    resid<-data-X%*%betaCur
    bn<-as.numeric(b0+.5*crossprod(resid,resid))
    return(rinvgamma(1,an,bn))
  }

  sampleBeta<-function(data){
    varCov.n<-solve(t(X)%*%X/sigma2Cur+Sigma0Inv)
    mu.n<-varCov.n%*%(t(X)%*%data/sigma2Cur+Sigma0Inv%*%mu0)
    return(mvrnorm(1,mu.n,varCov.n))
  }
  #choose starting y data
  y.prop<-rnorm(n)
  y.curr<-fn.comp.ystst(y.prop,X,l1obs,s1obs,psi,scaleEst,maxit)
  log.prop.den.curr <-log.prop.den(y.curr,X, projMatrix, l1obs, s1obs,fn.psi, fn.chi,n,p)
  ################################################
  #Start the gibbs sampler
  ################################################
  for(i in 1:total){
    #step one: [beta|sigma2,y]~N()
    betaCur<-sampleBeta(y.curr)
    #step 2 [sigma2|---]~IG
    sigma2Cur<-sampleSigma2(y.curr)
    betaSamples[i,]<-betaCur
    sigma2Samples[i]<-sigma2Cur
    #step 3: augmented data
    ySample<-fn.one.rep.y(y.curr,betaCur,sqrt(sigma2Cur),l1obs, s1obs,X, log.prop.den.curr, projMatrix,fn.psi,fn.chi, psi,scaleEst,maxit)
    y.curr<-ySample[[1]]
    yAccept[i]<-ySample[[2]]
    log.prop.den.curr<-ySample[[3]]
  }

  colnames(betaSamples)<-sapply(seq(1:p), FUN=function(x) paste('beta',x,sep=''))
  mcmcBetaSigma2<-mcmc(cbind(betaSamples[(nburn+1):total,],sigma2Samples[(nburn+1):total]))
  colnames(mcmcBetaSigma2)[p+1]<-'sigma2'
  fittedValues<-X%*%summary(mcmcBetaSigma2)$statistics[1:p,1]
  out<-list()
  out$mcmc<-mcmcBetaSigma2
  out$fitted.values<-fittedValues
  out$yAccept<-yAccept
  out$robust<-robust
  out$coef<-l1obs
  out$s<-s1obs
  out
}

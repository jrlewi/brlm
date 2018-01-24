#' MCMC algorithm for Heavy Tailed (t-distribution) Bayesian Linear Model
#'
#' Fits the heavy-tailed linear regression model
#' \eqn{y_i~N(x_i'beta,v_i)},
#' \eqn{v_i~inv-chi^2(nu,sigma2)}
#' \eqn{\beta~N(mu0,Sigma_0)}
#' \eqn{sigma2~IG(a0,b0)}.
#'
#' Uses Gibbs sampling to sample from the posterior under the above linear regression model. This model is equivalent to a linear regression with t-distributed errors.
#'
#' @inheritParams  bayesLm
#' @param parInit vector of initial values for \eqn{\beta, \sigma^2}. Default is \code{NULL} in which case initial values are chosen using a robust regression using the default settings in \code{\link[MASS]{rlm}}
#' @param nu degrees of freedom for t-distribution - assumed fixed
#' @param rwTune tuning parameter for the Metropolis-Hastings random walk used to sample \eqn{\sigma^2}
#' @return list with mcmc sample, mean fitted values, and logical acceptance vector for \eqn{\sigma^2} proposals.
#' @export
bayesTdistLm<-function(y, X,mu0, Sigma0 , a0, b0,parInit=NULL,nu, nkeep=1e4, nburn=1e3,rwTune=NULL){
  #y is the response
  #X is the design Matrix
  # parInit is the initial values for c(beta,sigma2)
  #mu0 prior mean of beta
  #Sigma0 is the var cov matrix of beta b~N(mu0,Sigma0)
  #a0, b0 prior parameters for sigma2
  #rwTune=MH rando  m walk for sigma2 st.dev

  p<-ncol(X)
  n<-length(y)
  total<-nkeep+nburn
  Sigma0Inv<-solve(Sigma0)

  if(is.null(parInit)){
    fit1<-rlm(X,y, maxit=400)
    parInit<-c(coef(fit1),summary(fit1)$sigma^2)
  }
  betaCur<-parInit[1:p]
  sigma2Cur<-parInit[p+1]

  if(is.null(rwTune)){
    fit1<-rlm(X,y, maxit=400)
    rwTune<-summary(fit1)$sigma/5
  }

  #Vectors to save samples
  betaSamples<-array(NA, dim=c(total,p))
  sigma2Samples<-numeric(total)
  acceptSigma2<-numeric(total)

  ############
  #Sampling functions
  ############
  #for V
  sampleV<-function(){
    as.vector(((y-X%*%betaCur)^2+nu*sigma2Cur)/rchisq(n,nu+1))
  }
  #for beta
  sampleBeta<-function()
  {
    vInv<-diag(vCur^-1)
    Var<-solve(t(X)%*%vInv%*%X+Sigma0Inv)
    Mean<-Var%*%(t(X)%*%vInv%*%y+Sigma0Inv%*%mu0)
    return(mvrnorm(1,mu=Mean, Sigma=Var))
  }

  #Sampling for sigma2
  #randomWalkMH for sigma2
  #log posterior of [sigma2|everything else]
  lpSigma2<-function(par){
    if(par>0){
      return((-a0-1+0.5*n*nu)*log(par)-b0/par-.5*nu*par*sum(1/vCur))}
    else return(-Inf)
  }
  mhStepSigma2<-function(){
    sigma2Prop<-rnorm(1, sigma2Cur, sd=rwTune)
    logMHrat<-lpSigma2(sigma2Prop)-lpSigma2(sigma2Cur)
    aprob<-min(1, exp(logMHrat))
    if(runif(1)<aprob){
      return(c(sigma2Prop,1))
    } else{
      return(c(sigma2Cur,0))
    }
  }
  ############################

  ############################
  #Start the Gibbs Sampler
  ############################

  for(i in 1:total){
    vCur<-sampleV()
    betaCur<-sampleBeta()
    outSigma2<-mhStepSigma2()
    sigma2Cur<-outSigma2[1]
    betaSamples[i,]<-betaCur
    sigma2Samples[i]<-sigma2Cur
    acceptSigma2[i]<-outSigma2[2]
  }
  #output
  colnames(betaSamples)<-sapply(seq(1:p), FUN=function(x) paste('beta',x,sep=''))
  mcmcBetaSigma2<-mcmc(cbind(betaSamples[(nburn+1):total,],sigma2Samples[(nburn+1):total]))
  colnames(mcmcBetaSigma2)[1:p] <- paste0('beta', 1:p)
  colnames(mcmcBetaSigma2)[p+1]<-'sigma2'
  fittedValues<-X%*%summary(mcmcBetaSigma2)$statistics[1:p,1]
  out<-list()
  out$mcmc<-mcmcBetaSigma2
  out$fitted.values<-fittedValues
  out$acceptSigma2<-acceptSigma2
  out
}

#' Fitting restricted likelihood model using importance sampling technique.

#'
#' Full model: \deqn{
#' \beta~N(\mu_0, \Sigma_0),
#' \sigma^2~IG(\alpha, \beta),
#' y~N(X\beta, \sigma^2)
#' }
#' For the restricted likelihood, conditioning is done on a pair of location and scale statistics \eqn{T(y)=(b(y), s(y))}. Current implementation allows for these to be a pair of M-estimators as implemented in \code{\link[MASS]{rlm}}
#'
#' Fill in some details
#'
#' @inheritParams rlDirectEval
#' @param X	a matrix or data frame containing the explanatory variables. The matrix should include a vector of 1's if intercept is desired.
#' @param instDist Placeholder to allow for more user specific instrumental distributions. Currently not used.
#' @param sdInstDist Vector of length 2 defining the scale for the instrumental distribution when instDist=NULL. In this case, the instrumental distribution for \eqn{\mu} and \eqn{\sigma^2} are independent normal and log normal distributions, respectively. The first value is then the standard deviation of the normal. The second is the standard deviation of \eqn{\log\sigma^2} (i.e. the \code{sdlog} argument in \code{\link[stats]{Lognormal}}). If left as NULL a mutltiple of 5 is used on the asymptotic variance covariance matrix. Importance sample weights should be examined to evaluate the appropriatness of this choice.
#' @param Nins number of samples from the instrumental distribution
#' @return list with elements \code{impSamps}, \code{w}, \code{fit}. These are the \code{Nins} by \code{length(beta)})+1} matrix of importance samples, the corresponding weights, and the fitted robust regression.
#' @author John R. Lewis \email{lewis.865@@osu.edu}
  #' @examples
  #' Example using the Belgium Phone calls data
  #' Conditioning statistics are simultaneous M-estimators of
  #' coefficients and scale. In the example below, Tukey's
  #' bisquare is used for the coefficients and median absolute
  #' deviation (MAD) is used for the scale. Note the sample sizes are dramatically reduces to run quickly.
  #' nForPrior<-3 #number of data points to use formulate prior
  #' data(MASS::phones)
  #' phones<-data.frame(cbind(phones$year, phones$calls))
  #' names(phones)<-c("year", "calls")
  #' phones$year<-phones$year-mean(phones$year) #center
  #' phonesFit<-phones[-c(1:nForPrior),]
  #' n<-nrow(phonesFit)
  #' X<-cbind(rep(1, n), phonesFit$year)
  #' y<-phonesFit$calls
  #' phonesPrior<-phones[c(1:nForPrior),]
  #' #formulate priors
  #' priorFit<-lm(calls~year, data=phonesPrior)
  #' summary(priorFit)
  #' mu0<-coef(priorFit)
  #' sigma2Hat<-summary(priorFit)$sigma^2
  #' alpha<-2
  #' beta<-(alpha-1) #sigma2Hat*(alpha-1)
  #' # choosing how much to inlfate fitted variance covariance
  #' matrix from the prior data fit
  #' inflate<-10
  #' Sigma0<-vinflate[inflate]*vcov(priorFit)
  #' restFit<-rlImportSamp(X=X,y=y, psi=psi.bisquare, scale.est='MAD', mu0=mu0, Sigma0=Sigma0, alpha=alpha, beta=beta, smooth=1,N=1e1,Nins=1e1, maxit=10000)
  #' plot(tst$w, cex=.1, pch=19)
  #' plot(density(tst$impSamps[,1],weights=tst$w))
  #' abline(v=b.obs[1], col=2)
  #' plot(density(tst$impSamps[,2],weights=tst$w))
  #' abline(v=b.obs[2], col=2)
  #' plot(density(tst$impSamps[,3],weights=tst$w))
  #' abline(v=s.obs, col=2)
#' @export
rlImportSamp<-function(X,y, psi, scale.est='Huber',k2=1.345, mu0, Sigma0, alpha, beta,instDist=NULL, sdInstDist =NULL, smooth=1,N,Nins, maxit=1000,...){

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
  rlm.obs<-rlm(X,y,psi=psi, scale.est=scale.est, k2=k2,maxit=maxit)
  b.obs<-rlm.obs$coefficients
  s.obs<-rlm.obs$s
  log.s.obs<-log(s.obs)
print(s.obs)
  if(!is.null(instDist)){stop('Currently instDist must be NULL')}

  meanInstBeta<<-b.obs
  sigmaInstBeta<<-5*vcov(rlm.obs)
  meanInstSigma2<<-s.obs^2
  sdInstSigma2<<-2*s.obs


  rlm.estimators<-function(Y){
    rlm.out<-rlm(X,Y,psi=psi, scale.est=scale.est, k2=k2,maxit=maxit)
return(c(as.numeric(rlm.out$coefficients), log(as.numeric(rlm.out$s))))
  }

Y<-matrix(rnorm(N*n,0,1),nrow=N,ncol=n)
rlm.matrix<-apply(Y,MARGIN=1, FUN=rlm.estimators)
h<-apply(rlm.matrix,1, hpi, binned=TRUE)
H<-diag(smooth*h)

kernel_density<-function(x1,x2){
  mean(dmvnorm(t(rlm.matrix),mean=c(x1,x2),sigma=H^2))
}

fn.posterior.estimate<-function(parms){
  l<-length(parms)
  betas<-parms[1:(l-1)]
  sigma2<-parms[l]
  kDensityEst<-log(kernel_density((b.obs-betas)/sqrt(sigma2),log.s.obs-.5*log(sigma2))/(sigma2^(.5*(l-1))))+dmvnorm(betas,mu0,Sigma0,log=TRUE)+log(dinvgamma(sigma2,alpha,beta))
  return(kDensityEst)
}


fn.logInsLikelihood<-function(parms){
  l<-length(parms)
  betas<-parms[1:(l-1)]
  sigma2<-parms[l]
  mvtnorm::dmvnorm(betas, mean=meanInstBeta, sigma=sigmaInstBeta, log=TRUE)+msm::dtnorm(sigma2,mean=meanInstSigma2, sd=sdInstSigma2,log=TRUE, lower=0)
}

fn.sample.instBeta<-function(){
  mvtnorm::rmvnorm(Nins, mean=meanInstBeta, sigma=sigmaInstBeta)
}

fn.sample.instSigma2<-function(){
  msm::rtnorm(Nins,mean=meanInstSigma2, sd=sdInstSigma2, lower=0)
}

insBetas<-fn.sample.instBeta()
insSigma2<-fn.sample.instSigma2()


insS<-cbind(insBetas,insSigma2)
logInsLikelihood<-apply(insS,1,FUN=fn.logInsLikelihood)

logPost<-apply(insS,1,fn.posterior.estimate)
wgts<-exp(logPost-logInsLikelihood)
wgts<-wgts/sum(wgts)

out<-list(impSamps=insS, w=wgts, fit=rlm.obs)
out
}


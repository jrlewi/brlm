#' Fitting restricted likelihood location-scale model using direct evaluation.
#'
#' Full model: \deqn{
#' \mu~N(\eta, \tau^2)
#' \sigma^2~IG(\alpha, \beta)
#' y_1, \dots, y_n~N(\mu, \sigma^2)
#' }
#' For the restricted likelihood, conditioning is done on a pair of location and scale statistics \eqn{T(y)=(l(y), s(y))} with the property that \eqn{l(\sigma y+\mu)=\sigma l(y)+\mu} and \eqn{s(\sigma y+\mu)=\sigma s(y)}.Current implementation allows for these to be a pair of M-estimators as implemented in \code{\link[MASS]{rlm}}
#'
#' Direct evaluation uses kenerel density estimation with a Gaussian kernel to estimate the restricted likelihood. \code{N} specifies the number of samples of the statistics to generate for the kernel density estimate. \code{\link[ks]{hpi}} is used to specify the initial bandwidths independently for the location and the scale. These can be multiplied by \code{smooth} to 'oversmooth' or 'undersmooth' the estimate. Oversmoothing may result in more stable estimates.
#'
#' A simple Riemann sum is used to determine the normalizing constant. Parameters for this numerical integration are specified by \code{mu_lims,sigma2_lims, length_mu}, and \code{length_sigma2}.

#' @param y  vector of data
#' @param psi the \code{psi} function used to define the location estimator. It is a function (possibly give by name) that is described in \code{\link[MASS]{rlm}}. Tuning constants are passed via \code{...}
#' @param scale.est,k2 specification of the scale estimator (issued by name- either \code{'MAD', 'Huber'}, or \code{'proposal 2'}) and the tuning constant for the Huber proposal 2 scale estimation. See these parameters in \code{\link[MASS]{rlm}} for more information.
#'
#' @param eta,tau prior mean and standard deviation for \eqn{mu}
#' @param alpha,beta prior shape and scale for \eqn{sigma^2}
#' @param mu_lims,sigma2_lims vectors of length 2 defining the limits for the numerical integration necessary to find the normalizing constant in Bayes' Theorem
#' @param length_mu,length_sigma2 length of the grid in each parameter to do the numerical integration (standard Riemann sum)
#' @param smooth scalar or vector of length 2 that will scale the initial bandwidth computed by \code{\link[ks]{hpi}}.
#'
#'@param N number of samples of the statistics to use for the kernel density estimation
#'@param maxit the limit on the number of IWLS iterations. Same as in \code{\link[MASS]{rlm}}
#' @param ... additional arguments to be passed to \code{psi}
#' @return A list of length 4: the joint posterior, the two marginals, and the bandwidths used for the kernel density estimate

#'@export
rlDirectEval<-function(y, psi,..., scale.est='Huber',k2=1.345, eta, tau, alpha, beta, mu_lims, sigma2_lims,length_mu, length_sigma2, smooth=1,N, maxit=1000){


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
rlm.obs<-rlm(y~1,psi=psi, scale.est=scale.est, k2=k2,maxit=maxit)
t.obs<-rlm.obs$coefficients
s.obs<-rlm.obs$s
log.s.obs<-log(s.obs)

#define the grid
#delta.mu<-.02  #distance between points on mu grid

mu.grid<-seq(mu_lims[1], mu_lims[2],length.out=length_mu)
delta.mu<-diff(mu.grid[1:2])

sigma2.grid<-seq(sigma2_lims[1], sigma2_lims[2],length.out=length_sigma2)
delta.sigma2<-diff(sigma2.grid[1:2])

#expand mu.grid and sigma.grid so that there are two vectors. Every possible pair of (mu,sigma2) is a single Row
#use expand.grid
expandGrid<-expand.grid(mu.grid,sigma2.grid)
#first column is the expanded mu vector, second column is the expanded sigma2 column



#computing the estimators----
estimators<-function(x){
  out<-rlm(x~1,psi=psi, scale.est=scale.est, k2=k2,maxit=maxit)
  return(c(as.numeric(out$coefficients), log(as.numeric(out$s))))
}

X<-matrix(rnorm(N*n,0,1),nrow=N,ncol=n) #each row is a sample from N(0,1)
est.matrix<-apply(X,MARGIN=1, FUN=estimators)
h1<-hpi(est.matrix[1,], binned=TRUE)
h2<-hpi(est.matrix[2,], binned=TRUE)
H<-diag(smooth*c(h1,h2))

#kd function using H as the bandwidth and est.matrix as the data
#using Guassian Kernel
kernel_density<-function(x1,x2, H){
  mean(dmvnorm(t(est.matrix),mean=c(x1,x2),sigma=H^2))
}

#unormalized posterior estimate function: estimates the posterior (unormalized) on the log scale
#using Gaussian
log.posterior<-function(mu,sigma2,H){
  #unormalized kernel density estimate of the postorior [mu, sigma2|t.obs,s.obs] on the log scale
  kDensityEst<-log(kernel_density((t.obs-mu)/sqrt(sigma2),log.s.obs-.5*log(sigma2),H)/sqrt(sigma2))+log(dnorm(mu,mean=eta,sd=tau))+log(dinvgamma(sigma2,alpha,beta))
  return(kDensityEst)
}

#-compute posterior
log.post.est<-mapply(FUN=log.posterior,expandGrid[,1],expandGrid[,2], MoreArgs=list(H=H))
normalizing.constant<-sum(exp(log.post.est))*delta.mu*delta.sigma2
post.density.final<-cbind(expandGrid[,1],expandGrid[,2],exp(log.post.est)/normalizing.constant)

#Marginalizing----
#Marginal for mu
post.mu.pdf<-sapply(mu.grid,FUN=function(x){delta.sigma2*sum(post.density.final[expandGrid[,1]==x,3])})
#Marginal for sigma2
post.sigma2.pdf<-sapply(sigma2.grid,FUN=function(x){delta.mu*sum(post.density.final[expandGrid[,2]==x,3])})


margMu<-cbind(mu.grid,post.mu.pdf)
colnames(margMu)<-c('mu', 'posterior')
margSigma2<-cbind(sigma2.grid,post.sigma2.pdf)
colnames(margSigma2)<-c('sigma2', 'posterior')
out<-list('jointPost'=post.density.final, "muPost"=margMu,"sigma2Post"=margSigma2, 'H'=diag(H))
out
}



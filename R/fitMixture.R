#' A standard Bayesian mixture model using a mixture of normal distributions for the likelihood.
#'
#' Fiting the following mixture model: \deqn{
#' \theta\sim N(\mu, \tau^2),
#' \sigma^2 \sim IG(\alpha, \beta),
#' y_i~(1-p)N(\theta, \sigma^2)+pN(\theta, c*\sigma^2)
#' }
#' where p and c are fixed by the user
#'
#'  This implementation is fit on a grid of \eqn{\theta} and \eqn{\sigma^2} using Riemann sum numerical integration to find the marginal distribution of the data in a similar way as in \code{\link{fitOrderStat}}. This is rather naive but is relatively quick - it only works well if the choice of grid on which the Riemann sum is done is chosen carefully. The grid is specified by \code{theta.lims}, \code{sigma2.lims}, \code{length.theta}, and \code{length.sigma2}. The grid should cover the region of non-negligble posterior mass. The number of grid points (\code{length.theta}/\code{length.sigma2}) must be large enough for good precision. However, larger values increase computation time.
#'
#' @param p,c fixed values for p and c in above model
#' @inheritParams fitOrderStat
#' @seealso \code{\link{fitOrderStat}}
#' @return A list of length 4: the joint posterior, the two marginals, and the posterior means
#' @export

fitMixture<-function(y, theta.lims,length.theta, sigma2.lims,length.sigma2, mu, tau, alpha, beta,p,c){

n<-length(y)
logPrior<-function(theta, sigma2, mu, tau, alpha, beta){
  dnorm(theta, mu, tau, log=TRUE)+log(dinvgamma(sigma2,alpha, beta))
  }

logMixtureLike<-function(y,theta,sigma2,p,c){
    sum(log((1-p)*dnorm(y,theta,sqrt(sigma2))+p*dnorm(y,theta,sqrt(c*sigma2))))
  }

  grid.theta<-seq(theta.lims[1],theta.lims[2], length.out=length.theta)
  dtheta<-grid.theta[2]-grid.theta[1]
  grid.sigma2<-seq(sigma2.lims[1],sigma2.lims[2], length.out=length.sigma2)
  dsigma2<-grid.sigma2[2]-grid.sigma2[1]
  grd<-expand.grid(grid.theta,grid.sigma2)

fn.log.post<-function(theta, sigma2){logPrior(theta, sigma2, mu, tau, alpha, beta)+logMixtureLike(y,theta,sigma2,p,c)}

  logUnormPost<-mapply(fn.log.post,grd[,1], grd[,2])
  post<-exp(logUnormPost+2*length(y))
  post<-post/(sum(post)*dtheta*dsigma2)
  margTheta<-sapply(grid.theta, FUN=function(x) sum(post[grd[,1]==x]))*dsigma2
  meanTheta<-sum(grid.theta*margTheta)*dtheta
  margSigma2<-sapply(grid.sigma2, FUN=function(x) sum(post[grd[,2]==x]))*dtheta
  meanSigma2<-sum(grid.sigma2*margSigma2)*dsigma2
  postMeans<-c(meanTheta, meanSigma2)
  posterior<-cbind(grd, post)
  colnames(posterior)<-c('theta','sigma2', 'posterior')
  margTheta<-cbind(grid.theta,margTheta)
  colnames(margTheta)<-c('theta','posterior')
  margSigma2<-cbind(grid.sigma2, margSigma2)
  colnames(margSigma2)<-c('sigma2','posterior')
  out<-list('jointPost'=posterior, "thetaPost"=margTheta, "sigma2Post"=margSigma2, 'postMeans'=postMeans)
  out
}

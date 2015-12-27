#' Order statistics likelihood as the restricted likelihood in a standard location-scale Bayesian model.
#'
#' Conditioning on a set of order statistics under a Guassian location and scale model.
#' The full model is: \deqn{
#' \theta\sim N(\mu, \tau^2),
#' \sigma^2 \sim IG(\alpha, \beta),
#' y_i~N(\theta, \sigma^2)
#' }
#' Conditioning is on a middle set of order statistics \eqn{(Y_{(k+1)}, \dots Y_{(n-k)})} so that if \eqn{k=0}, the full normal theory model is fit
#'
#' Fit is done a grid of \eqn{\theta} and \eqn{\sigma^2} using Riemann sum numerical integration to find the marginal distribution of the data. This is rather naive but is relatively quick - it only works well if the choice of grid on which the Riemann sum is done is chosen carefully. The grid is specified by \code{theta.lims}, \code{sigma2.lims}, \code{length.theta}, and \code{length.sigma2}. The grid should cover the region of non-negligble posterior mass. The number of grid points (\code{length.theta}/\code{length.sigma2}) must be large enough for good precision. However, larger values increase computation time.
#'
#' @param y  vector of data
#' @parm k numeric: choice of k determining the set of order statistics on which to condition.
#' @param theta.lims,sigma2.lims vectors of length 2 specifying the lower and upper limits for the numerical integration. Should span region of non-negligble posterior probability, otherwise the posterior will be incorrect
#' @param length.theta,length.sigma2 Number of grid points for each parameter used in the numerical integration
#' @param mu,tau mean and standard deviation of the normal prior distribution on \eqn{\theta}
#' @param alpha,beta prior shape and scale of inverse gamma prior on \eqn{\sigma^2}
#' @return A list of length 4: the joint posterior, the two marginals, and the posterior means
#' @export
fitOrderStat<-function(y,k, theta.lims,length.theta, sigma2.lims,length.sigma2, mu, tau, alpha, beta){
n<-length(y)

stopifnot(k+1<=n-k)
logPrior<-function(theta, sigma2, mu, tau, alpha, beta){
dnorm(theta, mu, tau, log=TRUE)+log(dinvgamma(sigma2, alpha, beta))}
logOrderLike<-function(theta,sigma2){
    n<-length(y)
    y<-sort(y)[(k+1):(n-k)]
    u<-y[1]
    v<-y[length(y)]
    sum(dnorm(y, theta, sqrt(sigma2), log=TRUE))+k*pnorm(u, theta, sqrt(sigma2), log=TRUE)+k*pnorm(v,theta,sqrt(sigma2), log=TRUE, lower.tail=FALSE)
  }

    grid.theta<-seq(theta.lims[1],theta.lims[2], length.out=length.theta)
    dtheta<-grid.theta[2]-grid.theta[1]
    grid.sigma2<-seq(sigma2.lims[1],sigma2.lims[2], length.out=length.sigma2)
    dsigma2<-grid.sigma2[2]-grid.sigma2[1]
    grd<-expand.grid(grid.theta,grid.sigma2)

    fn.log.post<-function(theta, sigma2){logPrior(theta, sigma2, mu, tau, alpha, beta)+logOrderLike(theta, sigma2)}

    logUnormPost<-mapply(fn.log.post,grd[,1], grd[,2])
    mn<-min(logUnormPost)
    post<-exp(logUnormPost) #+2*length(data)) # attempt to handle numerical issues
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

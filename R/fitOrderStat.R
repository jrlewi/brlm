#' Restricted likelihood Bayesian model using a middle set of order statistics to define the restricted likelihood.
#'
#' Conditioning on a set of order statistics under a Guassian location and scale model.
#' The full model is: \deqn{
#' \mu\sim N(\eta, \tau^2),
#' \sigma^2 \sim IG(\alpha, \beta),
#' y_i~N(\mu, \sigma^2)
#' }
#' Conditioning is on a middle set of order statistics \eqn{(Y_{(k+1)}, \dots Y_{(n-k)})}. That is, the likelihood used in the Bayesian update is the likelihood of these order statistics. When \eqn{k=0}, the model is the full normal theory model written above.
#'
#' Fit is done a grid of \eqn{\mu} and \eqn{\sigma^2} using Riemann sum numerical integration to find the marginal distribution of the data. This is rather naive but is relatively quick - it only works well if the choice of grid on which the Riemann sum is done is chosen carefully. The grid is specified by \code{mu_lims}, \code{sigma2_lims}, \code{length_mu}, and \code{length_sigma2}. The grid should cover the region of non-negligble posterior mass. The number of grid points (\code{length_mu}/\code{length_sigma2}) etast be large enough for good precision. However, larger values increase computation time.
#'
#' @param y  vector of data
#' @param k numeric: choice of k determining the set of order statistics on which to condition.
#' @param mu_lims,sigma2_lims vectors of length 2 specifying the lower and upper limits for the numerical integration. Should span region of non-negligble posterior probability, otherwise the posterior will be incorrect
#' @param length_mu,length_sigma2 Number of grid points for each parameter used in the numerical integration
#' @param eta,tau mean and standard deviation of the normal prior distribution on \eqn{\mu}
#' @param alpha,beta prior shape and scale of inverse gamma prior on \eqn{\sigma^2}
#' @return A list of length 4: the joint posterior, the two marginals, and the posterior means
#' @seealso \code{\link{fitMixture}}
#' @export
fitOrderStat<-function(y,k, mu_lims,length_mu, sigma2_lims,length_sigma2, eta, tau, alpha, beta){
n<-length(y)

stopifnot(k+1<=n-k)
logPrior<-function(mu, sigma2, eta, tau, alpha, beta){
dnorm(mu, eta, tau, log=TRUE)+log(dinvgamma(sigma2, alpha, beta))}
logOrderLike<-function(mu,sigma2){
    y<-sort(y)[(k+1):(n-k)]
    u<-y[1]
    v<-y[length(y)]
    sum(dnorm(y, mu, sqrt(sigma2), log=TRUE))+k*pnorm(u, mu, sqrt(sigma2), log=TRUE)+k*pnorm(v,mu,sqrt(sigma2), log=TRUE, lower.tail=FALSE)
  }

    grid.mu<-seq(mu_lims[1],mu_lims[2], length.out=length_mu)
    dmu<-grid.mu[2]-grid.mu[1]
    grid.sigma2<-seq(sigma2_lims[1],sigma2_lims[2], length.out=length_sigma2)
    dsigma2<-grid.sigma2[2]-grid.sigma2[1]
    grd<-expand.grid(grid.mu,grid.sigma2)

    fn.log.post<-function(mu, sigma2){logPrior(mu, sigma2, eta, tau, alpha, beta)+logOrderLike(mu, sigma2)}

    logUnormPost<-mapply(fn.log.post,grd[,1], grd[,2])
    mn<-min(logUnormPost)
    post<-exp(logUnormPost+2*n) # attempt to handle numerical issues
    post<-post/(sum(post)*dmu*dsigma2)
    margMu<-sapply(grid.mu, FUN=function(x) sum(post[grd[,1]==x]))*dsigma2
    meanmu<-sum(grid.mu*margMu)*dmu
    margSigma2<-sapply(grid.sigma2, FUN=function(x) sum(post[grd[,2]==x]))*dmu
    meanSigma2<-sum(grid.sigma2*margSigma2)*dsigma2
    postMeans<-c(meanmu, meanSigma2)
    posterior<-cbind(grd, post)
    colnames(posterior)<-c('mu','sigma2', 'posterior')
    margMu<-cbind(grid.mu,margMu)
    colnames(margMu)<-c('mu','posterior')
    margSigma2<-cbind(grid.sigma2, margSigma2)
    colnames(margSigma2)<-c('sigma2','posterior')
    out<-list('jointPost'=posterior, "muPost"=margMu, "sigma2Post"=margSigma2, 'postMeans'=postMeans)
    out
    }

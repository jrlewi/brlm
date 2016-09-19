#' Compute y**
#' Computes y** from a proposed y* value

#' @param yst data vector to be transformed to ystst satisfying the observed summary statistics
#' @param X design matrix
#' @param l1obs, s1obs: observed statsitics
#' @param psi: psi function returning psi(x)/x or psi'(x). default is psi.huber from \link{MASS}
#' @param scale.est: defines the scale estimator. Currently supports only 'Huber'
#' @param maxit number of iterations to find T(yst)

#' @return ystst a data vector such that T(ystst)=T(y_obs)

fn.comp.ystst<-function(yst,X,l1obs,s1obs,psi=psi.huber, scale.est='Huber', maxit=400){
  fit.yst<-rlm(X,yst,psi=psi, scale.est=scale.est, maxit=maxit)
  l1yst<-fit.yst$coefficients #(fit.yst)
  s1yst<-fit.yst$s
  res<-yst-X%*%l1yst
  ystst<-s1obs*res/s1yst+X%*%l1obs
  return(ystst)
}

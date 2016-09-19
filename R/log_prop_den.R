#' Caluclate the proposal dens on the log scale
#' Computes proposal density by computing the necessary Jacobian of the transformation from y* to y**
#' @param  y.prop proposed data vector satisfying T(y.prop)=T(y_obs)
#' @param  X design matrix
#' @param  l1obs, s1obs observed statsitics
#' @param  fn.psi custom psi function returning psi(x) or psi'(x) defining M-estimator for location ; default is fn.psi.huber
#' @param  fn.chi custom chi function returning chi(x) or chi'(x) defining M-estimator for scale; default is fn.chi.prop2
#' @param  n \code{length(y_obs)=length(y.prop)}
#' @param  p number of regression coefficients
#' @return The proposal density for y.prop on the log scale (assuming original vector sampled uniformly on unit sphere)

log.prop.den <-function(y.prop,X, proj, l1obs, s1obs,fn.psi,fn.chi,n,p){log(fn.attenuation(y.prop,X,proj, l1obs,s1obs,fn.psi, fn.chi))-(n-p-1)*log(fn.radius(y.prop,X))}

#' @rdname log.prop.den
#' @inheritParams fn.attenuation2
#' @details Designed for use within \code{fn.one.rep.y}. Two equivalent versions, the second uses \code{fn.attenuation2} and is faster.
log.prop.den2 <-function(y.prop,X, proj,Qt, l1obs, s1obs,fn.psi,fn.chi,n,p)
{
  log(fn.attenuation2(y.prop,X,proj,Qt, l1obs,s1obs,fn.psi, fn.chi))-(n-p-1)*log(fn.radius2(y.prop,proj))

}


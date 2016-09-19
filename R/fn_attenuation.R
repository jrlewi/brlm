#' Computes attenuation factor
#'
#' Compute the attenuation factor of the Jacobian in the restricted likelihood MCMC step
#' @param ystst a data vector such that T(ystst)=T(y_obs)
#' @param X design matrix
#' @param proj the projection matrix onto the deviation space (aka the least squares residual space)
#' @param l1obs s1obs: observed statsitics
#' @param fn.psi: custom psi function returning psi(x) or psi'(x) defining M-estimator for location; default is fn.psi.huber
#' @param fn.chi: custom chi function returning chi(x) or chi'(x) defining M-estimator for scale; default is fn.chi.prop2
fn.attenuation <- function(ystst,X, proj,l1obs,s1obs,fn.psi=fn.psi.huber, fn.chi=fn.chi.prop2)
{
  p<-ncol(X)
  n<-nrow(X)
  #W<-t(Vt)[,(p+1):n]
  grads.ystst<-fn.grads(ystst,X, l1obs, s1obs,fn.psi, fn.chi)
  grad.y.B<-grads.ystst[1:p,]
  grad.y.s<-grads.ystst[(p+1),]
  unit.a<-grad.y.s/sqrt(sum(grad.y.s*grad.y.s)) #normalize
  zstst<-proj%*%ystst
  unit.b<- zstst/sqrt(sum(zstst*zstst))
  atten1 <-abs(sum(unit.a * unit.b))
  qr.A<-qr(cbind(grad.y.s,t(grad.y.B)))
  Qmat<-qr.Q(qr.A)
  tmpPrj<-((diag(n)-proj)%*%Qmat)[,-1]
  tmpPrj2<-t(tmpPrj)%*%tmpPrj
  atten2<-sqrt(abs(det(tmpPrj2)))
  return(atten1*atten2)
}

#' @rdname fn.attenuation
#' @param Qt  Q transpose where Q is the orthonormalized X
#' @return Attenuations of the transformation to ystst. These attenuations make up part of the Jacobian.
#' @details Two equivalent functions computing the attenuation piece. The second one uses the transpose of the orthonormalized \code{X} as an input. This makes it somewhat quicker.
fn.attenuation2 <- function(ystst,X,proj, Qt,l1obs,s1obs,fn.psi=fn.psi.huber, fn.chi=fn.chi.prop2)
{
  #Qt is just the t(Q) in Q%*%t(Q)+W%*%t(W)=I; C(X)=C(Q); Q orthonormalized X
  p<-ncol(X)
  n<-nrow(X)
  grads.ystst<-fn.grads(ystst,X, l1obs, s1obs,fn.psi, fn.chi)
  grad.y.B<-grads.ystst[1:p,]
  grad.y.s<-grads.ystst[(p+1),]
  unit.a<-grad.y.s/sqrt(sum(grad.y.s*grad.y.s)) #normalize
  zstst<-proj%*%ystst
  unit.b<- zstst/sqrt(sum(zstst*zstst))
  atten1 <-abs(sum(unit.a * unit.b))
  if(!is.na(atten1)){
    if(p==1){ #patch to allow this function to work when p=1
      qr.A<-qr(cbind(grad.y.s,grad.y.B))
    } else {
      qr.A<-qr(cbind(grad.y.s,t(grad.y.B)))
    }
    Qmat<-qr.Q(qr.A)
    tmpPrj<-Qt%*%Qmat
    sings<-svd(tmpPrj)$d
    atten2<-prod(sings)
    return(atten1*atten2)} else {
      return(Inf) #return inf-this will cause log.prop.den2=Inf so the proposed ystst will be rejected
    }

}

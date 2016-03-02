#' Some utility functions for brlm


#computing the estimators----
#used in rlDirectEval
estimators<-function(x){
  out<-rlm(x~1,psi=psi, scale.est=scale.est, k2=k2,maxit=maxit)
  return(c(as.numeric(out$coefficients), log(as.numeric(out$s))))
}

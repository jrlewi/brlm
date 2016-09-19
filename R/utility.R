#' Some utility functions for brlm

#'Computing the Estimators for Conditions----
#'Funtioned used to compute the location and scale estimators using the \link{rlm} function. This function is used within rlDirectEval and is not designed to be called directly (it uses several variables not in the arguments list)
#'@parm x vector of data
estimators<-function(x){
  out<-rlm(x~1,psi=psi, scale.est=scale.est, k2=k2,maxit=maxit)
  return(c(as.numeric(out$coefficients), log(as.numeric(out$s))))
}



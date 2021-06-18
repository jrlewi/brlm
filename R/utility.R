#' #' Some utility functions for brlm
#'
#' #'Computing the Estimators for Conditions----
#' #'Funtion used to compute the location and scale estimators using the \link{rlm} function. This function is used within rlDirectEval and is not designed to be called directly (it uses several variables not in the arguments list)
#' #'@param x vector of data
#' estimators<-function(x){
#'   out<-rlm(x~1,psi=psi, scale.est=scale.est, k2=k2,maxit=maxit)
#'   return(c(as.numeric(out$coefficients), log(as.numeric(out$s))))
#' }


define_psi_chi <- function(regEst, scaleEst){
  ############################
  #define the psi and chi functions
  ############################
  if (regEst == 'Huber') {
    psi <- get('psi.huber') #internal
    fn.psi <- get('fn.psi.huber')

  } else {
    if (regEst == 'Tukey') {
      psi <- get('psi.bisquare') #internal
      fn.psi <- get('fn.psi.bisquare')
    } else {
      stop("only set up for Huber or Tukey regression estimates")
    }
  }

  if (scaleEst != 'Huber') {
    stop('scale estimate only set up for Hubers Prop2 ')
  }
  fn.chi <- fn.chi.prop2
  list(fn.psi, psi, fn.chi)
}

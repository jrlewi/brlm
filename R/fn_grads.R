#' Computes the gradients in the linear regression MCMC
#' @param y numeric vector of responses
#' @param X design Matrix
#' @param l1obs,s1obs observed value of slope and scale statistics
#' @param fn.psi psi function for M-estimation
#' @param fn.chi chi function for M-estimation

fn.grads<-function(y,X,l1obs, s1obs, fn.psi=fn.psi.huber, fn.chi=fn.chi.prop2){

  fn.psi <- get('fn.psi', mode = "function")
  fn.chi<- get('fn.chi', mode = "function")
  p<-ncol(X)
  resids<-y-X%*%l1obs
  stResids<-resids/s1obs
  psiPrime<-fn.psi(stResids, deriv=1)
  chiPrime<-fn.chi(stResids, deriv=1)
  A<-matrix(NA,p+1,p+1)
  psiColX<-apply(X, 2, FUN=function(x) x*psiPrime)
  A[1:p,1:p]<-s1obs*t(X)%*%psiColX
  A[1:p,p+1]<-t(X)%*%(psiPrime*resids)
  A[p+1,1:p]<-s1obs*t(X)%*%chiPrime
  A[p+1,p+1]<-sum(chiPrime*resids)
  prime1<-psiColX #for the slopes
  prime2<-c(chiPrime) #for the scale
  rhs<-s1obs*t(cbind(prime1,prime2))
  tryCatch({return(solve(A,rhs))}
           , error=function(err){print(paste('The Error:', err))
             return(matrix(rep(0, (p+1)*length(y)), p+1,length(y)))
           }
  )
}


#
# @param ... arguments to pass to the psi and chi functions (i.e. to specify tuning parameters),...
# arguments<-list(...)
# if (length(arguments)) {
#   pm <- pmatch(names(arguments), names(formals(fn.psi)), nomatch = 0L)
#   pm_chi <- pmatch(names(arguments), names(formals(fn.chi)), nomatch = 0L)
#   if (any(pm == 0L) | any(pm_chi==0L)){
#     warning("some of ... do not match")
#   }
#   pm <- names(arguments)[pm > 0L]
#   pm_chi<-names(arguments[pm_chi>0L])
#   formals(fn.psi)[pm] <- unlist(arguments[pm])
#   formals(fn.chi)[pm_chi] <- unlist(arguments[pm_chi])
# }

#' Psi functions
#' @param u numeric vector of evaluation points
#' @param k,c,b tunning parameters for each objective function
#' @param deriv 0 or 1: returns psi(x) or psi'(x)
#' @name psi
#' @details see \link{psi.huber}/\link{psi.bisquare} for more details. Note the difference in the returned value for \code{deriv=0} for this package. Functions re-written out of convenience for the development of \link{brlm}.
NULL
#' @rdname psi
fn.psi.huber<-function(u,k=1.345, deriv=0){
  if (!deriv)
    return(pmax(-k, pmin(k,u)))
  abs(u) <= k
}


#' @rdname psi
fn.psi.bisquare<-function(u, c = 4.685, deriv = 0)
{
  if (!deriv)
    return(u*(1 - pmin(1, abs(u/c))^2)^2)
  t <- (u/c)^2
  ifelse(t < 1, (1 - t) * (1 - 5 * t), 0)
}

#'@rdname psi
fn.chi.prop2<-function(u, k=1.345, deriv=0, b=.71016){
  if (!deriv)
    return(fn.psi.huber(u,k, 0)^2-b)
  2*fn.psi.huber(u,k,1)*fn.psi.huber(u,k)
}


#'@rdname psi
psi.quantile <-function(u, tau = 0.99, deriv = 0){
  if(deriv == 0){
    #note - quantiles have 0 psi', so try expectiles.

    val <- 2*u*tau*(u > 0) + 2*u*(1 - tau)*(u < 0)
    val <- val/u
  }
  if(deriv == 1){
    val <- 2*u*tau*(u > 0) + 2*u*(1 - tau)*(u < 0)
  }
  val
}

#'@rdname psi
fn.psi.quantile <-function(u, tau = 0.99, deriv = 0){
  if(deriv == 0){
    val <- 2*u*tau*(u > 0) + 2*u*(1 - tau)*(u < 0)
  }
  if(deriv == 1){
    val <- 2*tau*(u > 0) + 2*(1 - tau)*(u < 0)
  }
  val
}






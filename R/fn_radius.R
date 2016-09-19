#' Computes radius
#'
#'Compute the radius in deviation space
#' @param ystst a data vector such that T(ystst)=T(y_obs)
#' @param X design matrix
fn.radius<- function(ystst,X)
{
  Bls<-solve(t(X)%*%X)%*%t(X)%*%ystst
  resid<-ystst-X%*%Bls
  rad <-sqrt(sum(resid^2))
  return(rad)
}

#' @rdname fn.radius
#' @param proj the projection matrix onto the deviation space (aka the least squares residual space)
#' @return The radius in deviation space
#' @details Two equivalent ways to calculate radius in deviation space. The second uses the  projection matrix onto the deviation space as an argument is a bit quicker.
fn.radius2<-function(ystst,proj){
  resid<-proj%*%ystst
  rad <-sqrt(sum(resid^2))
  return(rad)
}

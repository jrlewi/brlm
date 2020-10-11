#abc_helpers
#' @rdname hier_abc
#'
#'
proposal_group_abc <- function(zl,
                             zminus,
                             yl,
                             fitsl,
                             condsd,
                             condMeanMult,
                             stepzl,
                             a0,
                             b0,
                             Xl,
                             XtX.l,
                             Beta,
                             bstar,
                             Sigma0Inv
                             ){

  zl <- fn.sample.Zl(zl,zminus,yl, fitsl, condsd,condMeanMult, step=stepzl)
  sigma2l <- fn.compute.sigma2(zl,a0,b0)
  betal <- fn.sample.beta.l(yl, Xl,XtX.l,Beta,sigma2l,bstar, Sigma0Inv)
  yl <- rnorm(length(yl), Xl%*%betal, sqrt(sigma2l))

  ######################################
  statistic <- compute_statistics() ###
  #########################
  list(yl=yl, betal=betal, sigma2l = sigma2l, statistic = statistic)
}

update_group_abc <- function(propoposal, current){


  ###






}


mh_abc_group <- function(){
  props <- proposal_group_abc()
  update_group_abc(current, props)


}

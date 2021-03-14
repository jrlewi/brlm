#abc_helpers

#' @rdname hier_abc
#'
compute_statistics <- function(y, X, psi, scaleEst, maxit){
  robust <- MASS::rlm(X, y, psi = psi, scale.est = scaleEst, maxit = maxit)
  c(robust$coef, robust$s)
  }




#' @rdname hier_abc
#'
proposal_group_abc <- function(params,
                               psi, scaleEst, maxit){


    yl <- params$yl
    a0 <- params$a0
    b0 <- params$b0
    Xl <- params$Xl
    Beta <- params$Beta
    b <- params$bstar
    Sigma0 <- params$Sigma0

    sigma2l <- rinvgamma(1, a0, b0)
    betal <- mvrnorm(1, Beta, b*Sigma0)
    yl <- rnorm(length(yl), Xl%*%betal, sqrt(sigma2l))
    statistic <- compute_statistics(yl, Xl, psi = psi, scaleEst=scaleEst, maxit = maxit)

    list(yl = yl, betal = betal, sigma2l = sigma2l, statistic = statistic)
}


stat_distance <- function(x, y){
    p <- length(x)
    if(p != length(y)) {stop("x,y, must be same length")}
    #change sigma to log sigma
    x[p] <- log(x[p])
    y[p] <- log(y[p])
    sqrt(sum((x - y)^2))
}

abc_distance_kernel <- function(d, tol){
    dnorm(d, sd = tol)
}

update_group_abc <- function(proposal, current, stat_obs, tol){

    y_curr <- current$yl
    stat_curr <-  current$statistic
    beta_curr <- current$betal
    sigma2_curr <- current$sigma2l
    y_prop <- proposal$yl
    stat_prop <-  proposal$statistic
    beta_prop <- proposal$betal
    sigma2_prop <- proposal$sigma2l

    stat_dist_prop <- stat_distance(stat_obs, stat_prop)
    stat_dist_current <- stat_distance(stat_obs, stat_curr)

    log_rat <- log(abc_distance_kernel(stat_dist_prop, tol)) -
                log(abc_distance_kernel(stat_dist_current, tol))


    mh_rat <- min(1, exp(log_rat))
    if(runif(1) < mh_rat){
      list(sample = proposal, accept = 1)
    } else {
      list(sample = current, accept = 0)
    }
      }


mh_abc_group <- function(params, current, stat_obs, tol,
                         psi, scaleEst, maxit){
    proposal <- proposal_group_abc(params,
                                   psi = psi, scaleEst=scaleEst, maxit = maxit)
    update_group_abc(proposal, current, stat_obs, tol)
}

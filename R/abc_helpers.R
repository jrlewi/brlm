#abc_helpers

#' @rdname hier_abc
#'
compute_statistics <- function(y, X, psi, scaleEst, maxit){
  robust <- MASS::rlm(X, y, psi = psi, scale.est = scaleEst, maxit = maxit)
  c(robust$coef, robust$s)
}


#' @rdname hier_abc
#'
proposal_group_abc <- function(params){

    zl <- params$zl
    zminus  <- params$zminus
    yl <- params$yl
    fitsl <- params$fitsl
    condsd <- params$condsd
    condMeanMult <- params$condMeanMult
    stepzl <- params$stepzl
    a0 <- params$a0
    b0 <- params$b0
    Xl <- params$Xl
    XtX.l <- params$XtX.l
    Beta <- params$Beta
    bstar <- params$bstar
    Sigma0Inv <- params$Sigma0Inv

    zl <- fn.sample.Zl(zl,zminus,yl, fitsl, condsd,condMeanMult, step=stepzl)
    sigma2l <- fn.compute.sigma2(zl,a0,b0)
    betal <- fn.sample.beta.l(yl, Xl,XtX.l,Beta,sigma2l,bstar, Sigma0Inv)
    yl <- rnorm(length(yl), Xl%*%betal, sqrt(sigma2l))
    statistic <- compute_statistics(yl, Xl, psi = psi, scaleEst=ScaleEst, maxit = maxit) ###

    list(yl = yl, betal = betal, sigma2l = sigma2l, zl = zl, statistic = statistic)
}


stat_distance <- function(x, y){
    p <- length(x)
    if(p != length(y)) {stop("x,y, must be same length")}
    #change sigma to log sigma
    x[p_1] <- log(x[p])
    y[p_1] <- log(y[p])
    sqrt(sum((x - y)^2))
}

abc_distance_kernel <- function(d, tol){
    dnorm(d, sd = tol)
}

update_group_abc <- function(proposal, current, Xl, stat_obs, tol){

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

    log_rat1 <- sum(dnorm(y_curr, Xl%*%beta_curr, sqrt(sigma2_curr),
                          log = TRUE)) -
                sum(dnorm(y_curr, Xl%*%beta_prop, sqrt(sigma2_prop),
                          log = TRUE))

    log_rat2 <- log(abc_distance_kernel(stat_dist_prop, tol)) -
                log(abc_distance_kernel(stat_dist_current, tol))


    mh_rat <- min(1, exp(rat1 + rat2))

    if(runif(1) < mh_rat){
      list(sample = proposal, accept = 1)
    } else {
      list(sample = current, accept = 0)
    }
      }


mh_abc_group <- function(params, current, Xl, stat_obs, tol){
    proposal <- proposal_group_abc(params)
    update_group_abc(proposal, current, Xl, stat_obs, tol)
}

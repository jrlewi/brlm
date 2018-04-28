#' Sampling Functions for random effects model
#'
re_sample_mu <- function(sigma2_0, J, tau2, mu_0, theta){
  sum_precisions <- 1/sigma2_0 + J/tau2
  new_sd <- sqrt(1/sum_precisions)
  new_mean <- (mu_0/sigma2_0 + sum(theta)/tau2)/sum_precisions
  rnorm(1, new_mean, new_sd)
}

re_sample_tau2 <- function(alpha_t, J, beta_t, theta, mu){
  new_shape <- alpha_t + 0.5*J
  new_scale <- beta_t + 0.5*sum((theta - mu)^2)
  rinvgamma(1, new_shape, new_scale)
}
#define nj and yj
re_sample_thetaj <- function(nj,yj, tau2, sigma2, mu){
  sum_precisions <- 1/tau2 + nj/sigma2
  new_sd <- sqrt(1/sum_precisions)
  new_mean <- (mu/tau2 + sum(yj)/sigma2)/sum_precisions
  rnorm(1, new_mean, new_sd)
}


re_sample_sigma2 <- function(alpha_s, beta_s, n, Y, theta){
  new_shape <- alpha_s + 0.5*sum(n)
  ss <- sum(mapply(FUN = function(yj, thetaj){
    sum((yj-thetaj)^2)
    }, yj = Y, thetaj = theta)
  )
  new_scale <- beta_s + 0.5*ss
  rinvgamma(1, new_shape, new_scale)
}


fn_one_rep_re <- function(sigma2_0,
                          J,
                          tau2,
                          mu_0,
                          theta,
                          alpha_t,
                          beta_t,
                          mu,
                          n,
                          Y,
                          sigma2,
                          alpha_s,
                          beta_s){
  mu <- re_sample_mu(sigma2_0, J, tau2, mu_0, theta)

  tau2 <- re_sample_tau2(alpha_t,J, beta_t, theta, mu)

  theta <- mapply(FUN = re_sample_thetaj, nj = n, yj = Y,   MoreArgs = list(tau2 = tau2, sigma2 = sigma2, mu = mu))

  sigma2 <- re_sample_sigma2(alpha_s, beta_s, n, Y, theta)

  out <- list(mu = mu, tau2 = tau2, theta = theta, sigma2 = sigma2)
out
}


#' Fit the Normal Theory Random Effects Model
#'
#' Fits the random effects model using standard Gibbs Sampling
#' \eqn{\mu ~ N(\mu_0, \sigma^2_0)}; \eqn{\tau^2 ~ IG(\alpha_t, \beta_t)}
#' \eqn{\theta_j~N(\mu, \tau^2), j = 1,2,...,J}; \eqn{\sigma^2~IG(\alpha_s, \beta_s)},
#' \eqn{y_{ij}~N(\theta_j), \sigma^2), i = 1,2,...,n_j}
#'

#' @param Y list; jth element is a vector of the y_ij
#' @param mu_0,sigma2_0 prior mean and variance for \eqn{\mu}
#' @param alpha_t,beta_t prior shape and scale for \eqn{\tau^2}
#' @param alpha_s,beta_s prior shape and scale for \eqn{\sigma^2}
#' @param nkeep number of iterations to keep
#' @param nburn number of iterations to toss
#' @return list with mcmc sample and mean fitted values
#' @export
re_model <- function(Y, mu_0, sigma2_0, alpha_t, beta_t,alpha_s, beta_s, nkeep=1e3, nburn=1e3){
  ntot <- nkeep + nburn
  J <- length(Y)
  n <- sapply(Y, length)
  #saving the parameters
  theta_samps <- matrix(NA, nrow = ntot, ncol = J)
  mu_samps <- tau2_samps <- sigma2_samps <- numeric(ntot)

  # get some reasonable initial values
  ybar <-  sapply(Y, mean)
  var_y <- sapply(Y, var)
  theta <- ybar
  mu <- mean(theta)
  tau2 <- var(theta)
  sigma2 <- mean(var_y)

  for(i in 1:ntot){
  samp <- fn_one_rep_re(sigma2_0,
                        J,
                        tau2,
                        mu_0,
                        theta,
                        alpha_t,
                        beta_t,
                        mu,
                        n,
                        Y,
                        sigma2,
                        alpha_s,
                        beta_s)
  theta <- samp$theta
  mu <- samp$mu
  tau2 <- samp$tau2
  sigma2 <- samp$sigma2

  theta_samps[i, ] <- theta
  mu_samps[i] <- mu
  tau2_samps[i] <- tau2
  sigma2_samps[i] <- sigma2
  }

  out <- list(theta = theta_samps[-(1:nburn),], mu = mu_samps[-c(1:nburn)], tau2 = tau2_samps[-c(1:nburn)], sigma2 = sigma2_samps[-c(1:nburn)])
  out
}


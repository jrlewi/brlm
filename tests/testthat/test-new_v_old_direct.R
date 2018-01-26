library(brlm)
library(MASS)
library(mvtnorm)

context('Old v. new direct sampling functions')

# rl_direct vs. rlDirect -----

data(newcomb)
y <- newcomb

rlm_tukey <- function(y){
  fit <- rlm(y~1,psi = psi.bisquare,  scale.est = 'proposal 2', maxit=50)
  c(fit$coefficients, fit$s)
}

rlm_huber <- function(y){
  fit <- rlm(y~1,psi = psi.huber,  scale.est = 'proposal 2', maxit=50)
  c(fit$coefficients, fit$s)
}

eta <- 23.6; tau <- 2.04; alpha <- 5; beta <- 10;
mu_lims <- c(20, 40); length_mu <- 10;
sigma2_lims <- c(.01, 60); length_sigma2 <- 10;
N <- 100

set.seed(1)
fit_new <- brlm:::rl_direct(y, statistic = rlm_tukey, eta, tau, alpha, beta, mu_lims, sigma2_lims,length_mu, length_sigma2, smooth=1,N)

set.seed(1)
fit_old <- brlm:::rlDirectEval(y, psi = psi.bisquare, scale.est='Huber',k2=1.345, eta = eta, tau = tau, alpha = alpha, beta = beta, mu_lims = mu_lims, sigma2_lims= sigma2_lims, length_mu = length_mu, length_sigma2 = length_sigma2, smooth = 1,N = N, maxit=1000)

set.seed(1)
fit_new2 <- brlm:::rl_direct(y, statistic = rlm_huber, eta, tau, alpha, beta, mu_lims, sigma2_lims,length_mu, length_sigma2, smooth=1,N)

set.seed(1)
fit_old2 <- brlm:::rlDirectEval(y, psi = psi.huber, scale.est='Huber',k2=1.345, eta = eta, tau = tau, alpha = alpha, beta = beta, mu_lims = mu_lims, sigma2_lims= sigma2_lims, length_mu = length_mu, length_sigma2 = length_sigma2, smooth = 1,N = N, maxit=1000)


test_that('new version rl_direct matches old', {
  expect_equal(fit_new, fit_old)
  expect_equal(fit_new2, fit_old2)
})



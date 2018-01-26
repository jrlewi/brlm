# Compare rl_importance and rlImport

library(MASS)
library(mvtnorm)
library(MCMCpack)
library(brlm)

data(newcomb)
y <- newcomb

rlm_tukey <- function(y){
  fit <- rlm(y~1,psi = psi.bisquare,  scale.est = 'proposal 2', maxit=50)
  c(fit$coefficients, fit$s)
}


eta <- 23.6; tau <- 2.04; alpha <- 5; beta <- 10;
mu_lims <- c(20, 40); length_mu <- 300;
sigma2_lims <- c(2, 60); length_sigma2 <- 300;
N <- 5000

# direct ----
set.seed(1)
fit_direct <- brlm:::rl_direct(y, statistic = rlm_tukey, eta, tau, alpha, beta, mu_lims, sigma2_lims,length_mu, length_sigma2, smooth=1,N)


#importance sampling ----
fit_rlm <- rlm(y~1, psi = psi.bisquare,  scale.est = 'Huber', maxit=1000)

rlm_tukey <- function(y, X){
  fit <- rlm(X, y ,psi = psi.bisquare,  scale.est = 'Huber', maxit=1000)
  c(fit$coefficients, fit$s)
}

cov_b <- 5*vcov(fit_rlm)
scale <- 2*fit_rlm$s

N <- 3000
Nins <- 1e4
X <- as.matrix(rep(1, length(y)))


set.seed(1)
fit_importance <- brlm::rl_importance(y, X , statistic = rlm_tukey, mu0 = eta, Sigma0 = as.matrix(tau^2), alpha = alpha, beta = beta, cov_b = cov_b, scale = scale, N = N, Nins = Nins)

#plot compare direct and importance ----
plot(fit_importance$w, cex = .1)

plot(fit_direct$muPost, type = 'l')
lines(density(fit_importance$impSamps[,1], weights = fit_importance$w), type = 'l', col = '2')

plot(fit_direct$sigma2Post, type = 'l', xlim = c(2,50))
lines(density(fit_importance$impSamps[,2], weights = fit_importance$w), type = 'l', col = '2')



# rl_importance vs. mcmc function

library(MASS)
library(mvtnorm)
library(MCMCpack)
library(brlm)



#plot function ---
plot_kd <- function(variable = 1, ylim = NULL){
  plot(density(fit_importance$impSamps[,variable], weights =fit_importance$w), ylim = ylim, lwd = 2, main = paste0('variable ', variable))
  lines(density(fitmcmc$mcmc[,variable]), col = 4)
}

# For Newcomb data----

data("newcomb")
y <- (newcomb - mean(newcomb))/sd(newcomb)
X <-as.matrix(rep(1, length(y)))
rlm_tukey <- function(y, X){
  fit <- rlm(X, y ,psi = psi.bisquare,  scale.est = 'Huber', maxit=1000)
  c(fit$coefficients, fit$s)
}
alpha <- .01
beta <- .01
fit_rlm <- rlm(y~1, psi = psi.bisquare,  scale.est = 'Huber', maxit=1000)
cov_b <- 2*vcov(fit_rlm)
scale <- .25*fit_rlm$s
mu0 <- 0
Sigma0 <- as.matrix(5)
N <- 3000
Nins <- 3e4

set.seed(1)
fit_importance <- brlm::rl_importance(y, X , statistic = 'rlm_tukey', mu0 = mu0, Sigma0 = Sigma0, alpha = alpha, beta = beta, cov_b = cov_b, scale = scale, N = N, Nins = Nins)
range(fit_importance$impSamps[,2])
plot(fit_importance$w, cex = .1)

set.seed(12)
fitmcmc <- brlm::restrictedBayesLm(y, X, regEst = 'Tukey', scaleEst = 'Huber', mu0 = mu0, Sigma0 = Sigma0, a0 = alpha, b0 = beta,sigma2Int = fit_rlm$s^2, nkeep = 5000, nburn = 10000)
coda::traceplot(fitmcmc$mcmc)


plot_kd(1, c(0, 7))
plot_kd(2, c(0, 10))







# with the phones data ----

data(phones)
phones <- data.frame(phones)
phones <- as.data.frame(scale(phones))
y <- phones$calls
X <- cbind(rep(1, length(phones$year)), phones$year)


rlm_tukey <- function(y, X){
  fit <- rlm(X, y ,psi = psi.bisquare,  scale.est = 'Huber', maxit=1000)
  c(fit$coefficients, fit$s)
}

fit_rlm <- fit <- rlm(X, y ,psi = psi.bisquare,  scale.est = 'Huber', maxit=1000)


mu0 <- c(0,0)
Sigma0 <- diag(c(2,2))
alpha <- .01; beta <-.01
cov_b <- 2*vcov(fit_rlm)
scale <- .75*fit_rlm$s
N <- 3000
Nins <- 3e4
set.seed(1)
fit_importance <- brlm::rl_importance(y, X, statistic = 'rlm_tukey', mu0 = mu0, Sigma0 = Sigma0, alpha = alpha, beta = beta, cov_b = cov_b, scale = scale, N = N, Nins = Nins)
range(fit_importance$impSamps[,1])
range(fit_importance$impSamps[,2])
range(fit_importance$impSamps[,3])
plot(fit_importance$w, cex = .1)


fitmcmc <- brlm::restrictedBayesLm(y, X, regEst = 'Tukey', scaleEst = 'Huber', mu0 = mu0, Sigma0 = Sigma0, a0 = alpha, b0 = beta, sigma2Int = fit_rlm$s^2, nkeep = 1e4)

coda::traceplot(fitmcmc$mcmc)

plot_kd(1, ylim = c(0, 2.2))
plot_kd(2, c(0, 2.2))
plot_kd(3, c(0, 2))


apply(fit_importance$impSamps, 2, function(x) sum(x*fit_importance$w)/sum(fit_importance$w))
apply(fitmcmc$mcmc, 2, mean)




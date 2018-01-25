# Compare rl_importance and rlImport

library(MASS)
library(mvtnorm)
library(MCMCpack)
library(brlm)

#plot function ---

plot_kd <- function(variable = 1, ylim = NULL){
  plot(density(fit_new$impSamps[,variable], weights = fit_new$w), ylim = ylim)
  lines(density(fit_old$impSamps[,variable], weights = fit_old$w), col = 2)
  lines(density(fitmcmc$mcmc[,variable]), col = 4)
}




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
alpha <- 4; beta <-3
hist(rinvgamma(1000, alpha, beta))
cov_b <- 5*vcov(fit_rlm)
scale <- 2*fit_rlm$s
N <- 5000
Nins <-5000
set.seed(1)
fit_new <- brlm::rl_importance(y, X, statistic = 'rlm_tukey', mu0 = mu0, Sigma0 = Sigma0, alpha = alpha, beta = beta, cov_b = cov_b, scale = scale, N = N, Nins = Nins)


set.seed(1)
fit_old <- brlm:::rlImportSamp(X, y, psi = psi.bisquare, scale.est='Huber', mu0 = mu0, Sigma0 = Sigma0, alpha = alpha, beta = beta, N = N, Nins = Nins,  maxit = 1000)


fitmcmc <- brlm::restrictedBayesLm(y, X, regEst = 'Tukey', scaleEst = 'Huber', mu0 = mu0, Sigma0 = Sigma0, a0 = alpha, b0 = beta, sigma2Int = fit_rlm$s^2, nkeep = 1e4)

coda::traceplot(fitmcmc$mcmc)

plot_kd(1, ylim = c(0, 2.2))
plot_kd(2, c(0, 2.2))
plot_kd(3, c(0, 2))


apply(fit_old$impSamps, 2, function(x) sum(x*fit_old$w)/sum(fit_old$w))
apply(fit_new$impSamps, 2, function(x) sum(x*fit_new$w)/sum(fit_new$w))

apply(fitmcmc$mcmc, 2, mean)


# plot(log(fit_new$w +.01), cex = .1)
# plot(log(fit_old$w +.01), cex = .1)



# Newcomb ----

data("newcomb")
set.seed(1)
y <- (newcomb - mean(newcomb))
X <-as.matrix(rep(1, length(y)))
rlm_tukey <- function(y, X){
  fit <- rlm(X, y ,psi = psi.bisquare,  scale.est = 'Huber', maxit=1000)
  c(fit$coefficients, fit$s)
}
alpha <- .01
beta <- .01
fit_rlm <- fit <- rlm(y~1, psi = psi.bisquare,  scale.est = 'Huber', maxit=1000)
cov_b <- 5*vcov(fit_rlm)
scale <- 2*fit_rlm$s
mu0 <- 0
Sigma0 <- as.matrix(5)
N <- 3000
Nins <-3000

fit_new <- brlm::rl_importance(y, X , statistic = 'rlm_tukey', mu0 = mu0, Sigma0 = Sigma0, alpha = alpha, beta = beta, cov_b = cov_b, scale = scale, N = N, Nins = Nins)

fit_old <-  brlm:::rlImportSamp(X, y, psi = psi.bisquare, scale.est='Huber', mu0 = mu0, Sigma0 = Sigma0, alpha = alpha, beta = beta, N = N, Nins = Nins,  maxit = 1000)

fitmcmc <- brlm::restrictedBayesLm(y, X, regEst = 'Tukey', scaleEst = 'Huber', mu0 = mu0, Sigma0 = Sigma0, a0 = alpha, b0 = beta,sigma2Int = fit_rlm$s^2, nkeep = 5000, nburn = 10000)
coda::traceplot(fitmcmc$mcmc)


plot_kd(1, c(0, 1))
plot_kd(2, c(0,.1))

range(fit_new$impSamps)
range(fit_old$impSamps)

apply(fit_old$impSamps, 2, function(x) sum(x*fit_old$w)/sum(fit_old$w))
apply(fit_new$impSamps, 2, function(x) sum(x*fit_new$w)/sum(fit_new$w))

apply(fitmcmc$mcmc, 2, mean)


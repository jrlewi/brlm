#rl_import vs. rlImport ----
library(brlm)
library(MASS)
library(mvtnorm)

context('Old v. new importance sampling functions')
data(phones)
phones <- data.frame(phones)
phones <- as.data.frame(scale(phones))
y <- phones$calls
X <- cbind(rep(1, length(phones$year)), phones$year)


rlm_tukey <- function(y, X){
  fit <- rlm(X, y ,psi = psi.bisquare,  scale.est = 'Huber', maxit=1000)
  c(fit$coefficients, fit$s)
}

fit_rlm <- rlm(X, y ,psi = psi.bisquare,  scale.est = 'Huber', maxit=1000)


mu0 <- c(0,0)
Sigma0 <- diag(c(2,2))
alpha <- .01; beta <-.01
cov_b <- 5*vcov(fit_rlm)
scale <- 2*fit_rlm$s
N <- 50
Nins <-50
set.seed(1)
fit_new <- brlm::rl_importance(y, X, statistic = rlm_tukey, mu0 = mu0, Sigma0 = Sigma0, alpha = alpha, beta = beta, cov_b = cov_b, scale = scale, N = N, Nins = Nins)


set.seed(1)
fit_old <- brlm:::rlImportSamp(X, y, psi = psi.bisquare, scale.est='Huber', mu0 = mu0, Sigma0 = Sigma0, alpha = alpha, beta = beta, N = N, Nins = Nins,  maxit = 1000)

test_that('new version rl_import matches old', {
  expect_equal(fit_new$impSamps, fit_old$impSamps, check.attributes = FALSE)
  expect_equal(fit_new$w, fit_old$w,  check.attributes = FALSE)
})


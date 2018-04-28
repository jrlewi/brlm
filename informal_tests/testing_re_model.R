library(rstan)
library(MCMCpack)
library(brlm)


gen_random_effects <- function(n, mu, tau2, sigma2){
  sim_data <- sapply(n, function(nn){
    theta_j <- rnorm(1, mu, sqrt(tau2))
    yj <- rnorm(nn, theta_j, sqrt(sigma2))
    c(theta = theta_j, y = yj)
  })
  theta <- sapply(sim_data, function(x) x["theta"])
  Y <- lapply(sim_data, function(x) x[-1])
  list(theta = theta, Y = Y)
  }
Ns <- c(rep(10, 10), rep(50, 20), rep(100, 20))
sim_data <- gen_random_effects(n =  Ns, mu = 0, tau2 = 2, sigma2 = 1)


# my code

post_samps <- re_model(sim_data$Y, mu_0 = 0, sigma2_0 = 3, alpha_t = 2, beta_t = 1,alpha_s = 2, beta_s = 1, nkeep=1e3, nburn=1e3)


# stan code

Y <- unlist(sim_data$Y)
ll <- NA
for(i in 1:length(Ns)){
  ll <- c(ll, rep(i, Ns[i]))
}
ll <-ll[-1]

stan_data <- list(Y = Y, N = length(Y), mu_0 = 0,  sigma_0 = sqrt(3), alpha_t = 2, beta_t = 1,alpha_s = 2, beta_s = 1, J = length(Ns), ll = ll )


# fit re model ----
  fit_re <- stan(file='informal_tests/re_stan.stan', data=stan_data, chains=1, iter = 2e3, refresh = -1, verbose = FALSE)

theta_stan <- extract(fit_re, 'theta')$theta

jj <- 50
plot(density(post_samps$theta[,jj]))
lines(density(theta_stan[,jj]), col = 2)
abline(v = sim_data$theta[jj])

mu_stan <- extract(fit_re, 'mu')$mu
plot(density(post_samps$mu), type = 'l')
lines(density(mu_stan), col = 2)

tau_stan <- extract(fit_re, 'tau')$tau
plot(density(sqrt(post_samps$tau2)), type = 'l')
lines(density(tau_stan), col = 2)

sigma_stan <- extract(fit_re, 'sigma')$sigma
plot(density(sqrt(post_samps$sigma2)), type = 'l')
lines(density(sigma_stan), col = 2)







mns <- apply(post_samps$theta,2, mean)
sds <- sapply(Y, sd)

matplot(post_samps$theta, type = 'l')
plot(sim_data$theta, apply(post_samps$theta, 2, mean), pch = 19)
abline(0,1, col = 2, lwd = 2)

plot(post_samps$mu, type = 'l')
abline(h = mean(sapply(Y, mean)), col = 2)
abline(h = 0, col = 2)
plot(sqrt(post_samps$tau2), type = 'l')
abline(h = sd(mns), col = 2)
abline(h = sqrt(2), col = 2)


plot(sqrt(post_samps$sigma2), type = 'l')
abline(h = 1, col = 2)
abline(h = mean(sds), col = 2)

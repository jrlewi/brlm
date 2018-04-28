data {
 int<lower = 0> N; //total length
 vector[N] Y;
 real mu_0;
 real sigma_0;
 real alpha_t;
 real beta_t;
 real alpha_s;
 real beta_s;
 int J; // number of groups
 int<lower = 1, upper = J> ll[N]; //group label
}

parameters {
  real mu;
  real<lower = 0> tau2;
  vector[J] theta;
  real<lower = 0> sigma2;
}

transformed parameters{
  real<lower=0> sigma;
  real<lower=0> tau;
    sigma =  pow(sigma2, 0.5);
    tau =  pow(tau2, 0.5);
}

model {
 mu ~ normal(mu_0, sigma_0);
 tau2 ~ inv_gamma(alpha_t, beta_t);
 sigma2 ~ inv_gamma(alpha_s, beta_s);
 theta ~ normal(mu, tau);
 for(k in 1:N){
   Y[k] ~ normal(theta[ll[k]], sigma);
 }
}


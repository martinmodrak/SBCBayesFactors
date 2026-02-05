data {
  int<lower=0> N;
  array[N] int<lower=0> y;
  real<lower=0> intercept_prior_mean;
  real<lower=0> intercept_prior_sd;
}

parameters {
  real intercept;
}


model {
  target += normal_lpdf(intercept | intercept_prior_mean, intercept_prior_sd);
  target += poisson_log_lpmf(y | intercept);
}

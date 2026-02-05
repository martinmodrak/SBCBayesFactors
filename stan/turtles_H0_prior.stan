data {
  int<lower = 1> N;
  array[N] int<lower = 0, upper = 1> y;
  vector[N] x;
  real<lower=0> alpha0_prior_width;
  real<lower=0> alpha1_prior_width;
}

parameters {
  real alpha0;
  real alpha1;
}

model {
  // priors
  // priors
  if(alpha0_prior_width > 0) {
    target += normal_lpdf(alpha0 | 0, alpha0_prior_width);
  }
  if(alpha1_prior_width > 0) {
    target += normal_lpdf(alpha1 | 0, alpha1_prior_width);
  }
  // likelihood
  for (i in 1:N) {
    target += bernoulli_lpmf(y[i] | Phi(alpha0 + alpha1 * x[i]));
  }
}

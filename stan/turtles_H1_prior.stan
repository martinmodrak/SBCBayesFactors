data {
  int<lower = 1> N;
  array[N] int<lower = 0, upper = 1> y;
  vector[N] x;
  int<lower = 1> C;
  array[N] int<lower = 1, upper = C> clutch;
  real<lower=0> alpha0_prior_width;
  real<lower=0> alpha1_prior_width;
}

parameters {
  real alpha0;
  real alpha1;
  vector[C] b_raw;
  real<lower = 0> sigma2;
}

transformed parameters {
  vector[C] b;
  real<lower = 0> sigma = sqrt(sigma2);
  b = sigma * b_raw;
}

model {
  // priors
  if(alpha0_prior_width > 0) {
    target += normal_lpdf(alpha0 | 0, alpha0_prior_width);
  }
  if(alpha1_prior_width > 0) {
    target += normal_lpdf(alpha1 | 0, alpha1_prior_width);
  }
  target += - 2 * log1p(sigma2); // p(sigma2) = 1 / (1 + sigma2) ^ 2

  // random effects
  target += normal_lpdf(b_raw | 0, 1);

  // likelihood
  for (i in 1:N) {
    target += bernoulli_lpmf(y[i] | Phi(alpha0 + alpha1 * x[i] + b[clutch[i]]));
  }
}

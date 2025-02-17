data {
  int<lower = 1> N;
  array[N] int<lower = 0, upper = 1> y;
  vector[N] x;
  int<lower = 1> C;
  array[N] int<lower = 1, upper = C> clutch;
  real<lower=0> prior_width;
  int<lower=0,upper=1> link;
}

parameters {
  real alpha0_raw;
  real alpha1_raw;
  vector[C] b_raw;
  real<lower = 0> sigma;
}

transformed parameters {
  vector[C] b;
  real alpha0 = prior_width * alpha0_raw;
  real alpha1 = prior_width * alpha1_raw;
  b = sigma * b_raw;
}

model {
  // priors
  target += normal_lpdf(sigma | 0, 1);

  target += normal_lpdf(alpha0_raw | 0, 1);
  target += normal_lpdf(alpha1_raw | 0, 1);

  // random effects
  target += normal_lpdf(b_raw | 0, 1);

  // likelihood
  for (i in 1:N) {
    if(link == 0) {
      target += bernoulli_lpmf(y[i] | Phi(alpha0 + alpha1 * x[i] + b[clutch[i]]));
    } else {
      target += bernoulli_logit_lpmf(y[i] | alpha0 + alpha1 * x[i] + b[clutch[i]]);
    }
  }
}

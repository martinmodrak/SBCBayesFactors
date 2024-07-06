data {
  int<lower=1> N;
  vector[N] y;
  real r;
}

parameters {
  real delta;
  real log_sigma;
}

transformed parameters {
  real sigma = exp(log_sigma);
  real mu = delta * sigma;
}

model {
  target += normal_lpdf(y | mu, sigma);
  target += cauchy_lpdf(delta | 0, r);
  target += -2 * log_sigma;
}

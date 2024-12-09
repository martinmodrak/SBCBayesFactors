data {
  int<lower=1> N;
  vector[N] y;
  real r;
}

parameters {
  real delta;
  real<lower=0> sigma2;
}

transformed parameters {
  real sigma = sqrt(sigma2);
  real mu = delta * sigma;
}

model {
  target += normal_lpdf(y | mu, sigma);
  target += cauchy_lpdf(delta | 0, r);
  target += -log(sigma2);
}

data {
  int<lower=1> N;
  vector[N] y;
}

parameters {
  real<lower=0> sigma;
}

model {
  target += normal_lpdf(y | 0, sigma);
  target += -2 * log(sigma);
}

generated quantities {
  real mu = 0;
}

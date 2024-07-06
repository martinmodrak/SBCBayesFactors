data {
  int<lower=1> N;
  vector[N] y;
}

parameters {
  real log_sigma;
}

transformed parameters {
  real sigma = exp(log_sigma);
}

model {
  target += normal_lpdf(y | 0, sigma);
  target += -2 * log_sigma;
}

generated quantities {
  real mu = 0;
}

data {
  int<lower=1> N;
  vector[N] y;
}

parameters {
  real<lower=0> sigma2;
}

transformed parameters {
  real sigma = sqrt(sigma2);
}

model {
  target += normal_lpdf(y | 0, sigma);
  target += -log(sigma2);
}

generated quantities {
  real mu = 0;
}

data {
  int<lower=0> N;
  vector[N] y;
  real<lower=0> intercept_prior_sd;
  real<lower=0> sigma_prior;
}

parameters {
  real intercept;
  real<lower=0> sigma;
}


model {
  if(intercept_prior_sd != 0) {
    target +=  normal_lpdf(intercept | 0, intercept_prior_sd);
  }
  if(sigma_prior != 0) {
    target +=  normal_lpdf(sigma | 0, sigma_prior) - normal_lccdf(0 | 0, 1); //norm. constant for the sigma prior;;
  }
  target += normal_lpdf(y | intercept, sigma);
}

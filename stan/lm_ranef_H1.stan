data {
  int<lower=0> N;
  vector[N] y;
  int<lower=1> N_groups;
  array[N] int<lower=1,upper=N_groups> groups;
  real<lower=0> intercept_prior_sd;
  real<lower=0> sigma_prior;
}

parameters {
  real intercept;
  real<lower=0> sigma_raw;
  vector[N_groups] group_z;
  real<lower=0> group_sd;
}

transformed parameters {
  vector[N_groups] group_r = group_z * group_sd;
  real sigma;
  if(sigma_prior == 0) {
    sigma = sqrt(sigma_raw);
  } else {
    sigma = sigma_raw;
  }
}

model {
  if(intercept_prior_sd != 0) {
    target +=  normal_lpdf(intercept | 0, intercept_prior_sd);
  }
  if(sigma_prior != 0) {
    target += normal_lpdf(sigma | 0, sigma_prior) -normal_lccdf(0 | 0, 1); //norm. constant for the sigma prior;;
  } else {
    target += -1 * log(sigma_raw);
  }
  target += std_normal_lpdf(group_z);
  target += normal_lpdf(group_sd | 0, 1)
            -normal_lccdf(0 | 0, 1); //norm. constant for the sigma prior;
  target += normal_lpdf(y | intercept + group_r[groups], sigma);
}

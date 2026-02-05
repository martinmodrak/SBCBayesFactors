data {
  int<lower=0> N;
  int<lower=1> N_groups;
  array[N] int<lower=1,upper=N_groups> groups;
  array[N] int<lower=0> y;
  real<lower=0> intercept_prior_mean;
  real<lower=0> intercept_prior_sd;
}

parameters {
  real intercept;
  vector[N_groups] group_z;
  real<lower=0> group_sd;
}

transformed parameters {
  vector[N_groups] group_r = group_z * group_sd;
}

model {
  target += normal_lpdf(intercept | intercept_prior_mean, intercept_prior_sd);
  target += std_normal_lpdf(group_z);
  target += normal_lpdf(group_sd | 0, 1)
            -normal_lccdf(0 | 0, 1); //norm. constant for the sigma prior;
  target += poisson_log_lpmf(y | intercept + group_r[groups]);
}

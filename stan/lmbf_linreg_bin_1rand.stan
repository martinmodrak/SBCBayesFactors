data {
  int<lower=0> N;   // number of data items
  int<lower=0> K;   // number of binary predictors
  array[N, K] int<lower=0, upper=1> X;   //predictor matrix

  int<lower=0> G;  // N of Random groups
  array[N] int<lower=1,upper=G> group;

  vector[N] y;      // outcome vector
  real<lower=0> r_fixed;
  real<lower=0> r_random;
}

transformed data {
  matrix[N, K] X_prime = -(to_matrix(X) - 0.5) * sqrt(2);
}

parameters {
  real alpha;           // intercept
  vector[K] beta_std;       // coefficients for predictors
  real<lower=0> sigma2;  // error scale
  vector[G] ranef_std;
  real<lower=0> g_ranef;
}

transformed parameters {
  real<lower=0> sigma = sqrt(sigma2);
  vector[K] beta = beta_std * sigma;       // coefficients for predictors
  vector[G] ranef = ranef_std * sigma * r_random * sqrt(g_ranef);
}

model {
  vector[N] Xbeta;
  target += -log(sigma2);
  if(K > 0) {
    target += cauchy_lpdf(beta_std | 0, r_fixed);
    Xbeta = X_prime * beta;
  } else {
    Xbeta = rep_vector(0, N);
  }
  target += inv_gamma_lpdf(g_ranef | 0.5, 0.5);
  target += std_normal_lpdf(ranef_std);
  target += normal_lpdf(y | Xbeta + alpha + ranef[group], sigma);
}

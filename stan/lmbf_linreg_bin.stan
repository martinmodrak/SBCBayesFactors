data {
  int<lower=0> N;   // number of data items
  int<lower=0> K;   // number of binary predictors
  array[N, K] int<lower=0, upper=1> X;   //predictor matrix

  vector[N] y;      // outcome vector
  real<lower=0> r_fixed;
}

transformed data {
  matrix[N, K] X_prime = -(to_matrix(X) - 0.5) * sqrt(2);
}

parameters {
  real alpha;           // intercept
  vector[K] beta_std;       // coefficients for predictors
  real<lower=0> sigma2;  // error scale
}

transformed parameters {
  real<lower=0> sigma = sqrt(sigma2);
  vector[K] beta = beta_std * sigma;       // coefficients for predictors
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
  target += normal_lpdf(y | Xbeta + alpha, sigma);
}

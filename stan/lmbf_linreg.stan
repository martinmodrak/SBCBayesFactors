data {
  int<lower=0> N;   // number of data items
  int<lower=0> K;   // number of predictors
  matrix[N, K] X;   // predictor matrix
  vector[N] y;      // outcome vector
  real<lower=0> r;
  int prior_type;
}

transformed data {
  matrix[K, K] inv_xt_x = inverse(transpose(X) * X);
  vector[K] beta0 = rep_vector(0, K);
}

parameters {
  real alpha;           // intercept
  vector[K] beta_std;       // coefficients for predictors
  real<lower=0> sigma2;  // error scale
  real<lower=0> g;
}

transformed parameters {
  real<lower=0> sigma = sqrt(sigma2);
  vector[K] beta = beta_std * sigma;       // coefficients for predictors
}

model {
  vector[N] Xbeta;
  target += -log(sigma2);
  target += inv_gamma_lpdf(g | 0.5, 0.5 * r);
  if(K > 0) {
    if(prior_type == 0) {
      target += multi_normal_lpdf(beta_std | beta0, N * g * sigma2 * inv_xt_x);
    } else {
      target += cauchy_lpdf(beta_std | 0, r);
    }
    Xbeta = X * beta;
  } else {
    Xbeta = rep_vector(0, N);
  }
  target += normal_lpdf(y | Xbeta + alpha, sigma);
}

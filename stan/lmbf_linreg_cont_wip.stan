data {
  int<lower=0> N;   // number of data items
  int<lower=0> K_cont;   // number of continous predictors
  matrix[N, K_cont] X_cont;   // continous predictor matrix
  // int<lower=0> K_disc;   // number of discrete/group predictors
  // array[K_disc] int<lower=2> G_disc; // number of categories
  // array[N, K_disc] int<lower=1,upper= max(K_disc)> X_dics; //group membership index

  vector[N] y;      // outcome vector
  real<lower=0> r_cont;
  int prior_type;
}

transformed data {
  matrix[N, K_cont] X_cont_centered;
  matrix[K_cont, K_cont] inv_xt_x;
  vector[K_cont] beta0 = rep_vector(0, K_cont);

  for(k in 1:K_cont) {
    X_cont_centered[, k] = X_cont[, k] - mean(X_cont[,k]);
  }
  inv_xt_x = inverse(transpose(X_cont_centered) * X_cont_centered);
}

parameters {
  real alpha;           // intercept
  vector[K_cont] beta_std;       // coefficients for predictors
  real<lower=0> sigma2;  // error scale
  real<lower=0> g;
}

transformed parameters {
  real<lower=0> sigma = sqrt(sigma2);
  vector[K_cont] beta = beta_std * sigma;       // coefficients for predictors
}

model {
  vector[N] Xbeta;
  target += -log(sigma2);
  target += inv_gamma_lpdf(g | 0.5, 0.5 * r_cont);
  if(K_cont > 0) {
    if(prior_type == 0) {
      target += multi_normal_lpdf(beta_std | beta0, N * g * sigma2 * inv_xt_x);
    } else {
      target += cauchy_lpdf(beta_std | 0, r_cont);
    }
    Xbeta = X_cont_centered * beta;
  } else {
    Xbeta = rep_vector(0, N);
  }
  target += normal_lpdf(y | Xbeta + alpha, sigma);
}

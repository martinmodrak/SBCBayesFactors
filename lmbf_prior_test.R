model <- stan_model("stan/lmbf_linreg.stan")

r_scale <- sqrt(2) / 4

# Null dataset
set.seed(455223)
N <- 40
df <- data.frame(x1 = rnorm(N), x2 = exp(rnorm(N))) %>%
  mutate(y = rnorm(N))


bf_H1 <- lmBF(y ~ x1 + x2, data = df, rscaleCont =  r_scale)
bf_H1
#BF: 0.1562293 ±0.01%

# Fit with rstan
data_stan_base <- list(N = N, r_cont = r_scale, y = df$y, prior_type = 0)

data_stan_H0 <- data_stan_base
data_stan_H0$X_cont <- array(dim = c(N, 0))
data_stan_H0$K_cont <- 0

data_stan_H1 <- data_stan_base
mm_H1 <- model.matrix(y ~ x1 + x2, df)[, -1, drop = FALSE]
data_stan_H1$X_cont <- mm_H1
data_stan_H1$K_cont <- ncol(mm_H1)

res_H0 <- sampling(model, data = data_stan_H0, cores = 1, iter = 16000, warmup = 1000, thin = 10)
res_H1 <- sampling(model, data = data_stan_H1, cores = 1, iter = 16000, warmup = 1000, thin = 10)


# Bridgesampling
bs_H0 <- bridgesampling::bridge_sampler(res_H0)
bs_H1 <- bridgesampling::bridge_sampler(res_H1)
bf_bs <- bridgesampling::bayes_factor(bs_H1, bs_H0)

bf_bs
bf_H1


library(bridgesampling)
library(rstan)

m_H0 <- stan_model("stan/lm_ranef_H0.stan")
m_H1 <- stan_model("stan/lm_ranef_H1.stan")

# The seed gives some variability, and not all seeds exhibit the problem
# but the problem is quite reliable
seed <- 14611187
set.seed(seed)

N <- 40
N_groups <- 20

groups <- rep(1:N_groups, length.out = N)

intercept <- 1
group_sd <- 0.1
group_r <- rnorm(N_groups, sd = group_sd)

# Reducing sigma seems to result in stronger changes in BF, but not heavily tested
sigma <- 0.1

mu <- intercept + group_r[groups]
y <- rnorm(N, mu, sigma)

stan_data <- list(N = N,
                  N_groups = N_groups,
                  groups = groups,
                  y = y,
                  # zeroes make the prior improper - flat for intercept and 1/sigma^2 for residual sigma
                  intercept_prior_sd = 0,
                  sigma_prior = 0)

fit_H0 <- sampling(m_H0, data = stan_data,
                   iter = 15500, warmup = 500, chains = 4, seed = 1, refresh = 0)

fit_H1 <- sampling(m_H1, data = stan_data,
                   iter = 15500, warmup = 500, chains = 4, seed = 1, refresh = 0, control = list(adapt_delta = 0.95))

bridge_H0 <- bridge_sampler(fit_H0)
bridge_H1 <- bridge_sampler(fit_H1)

stan_data_prior <- stan_data
stan_data_prior$intercept_prior_sd <- 1
stan_data_prior$sigma_prior <- 1
fit_H0_prior <- sampling(m_H0, data = stan_data_prior,
                         iter = 15500, warmup = 500, chains = 4, seed = 1, refresh = 0)

fit_H1_prior <- sampling(m_H1, data = stan_data_prior,
                         iter = 15500, warmup = 500, chains = 4, seed = 1, refresh = 0, control = list(adapt_delta = 0.95))

bridge_H0_prior <- bridge_sampler(fit_H0_prior)
bridge_H1_prior <- bridge_sampler(fit_H1_prior)

error_measures(bridge_H0)$percentage
error_measures(bridge_H1)$percentage


error_measures(bridge_H0_prior)$percentage
error_measures(bridge_H1_prior)$percentage



bf(bridge_H1, bridge_H0)
bf(bridge_H1_prior, bridge_H0_prior)

#set.seed(168455)
N <- 100

# Reliably worse SBC
# probs_true <- seq(0.01,0.99, length.out = N)
# probs_computed <- pmin(1, probs_true + 0.2)


probs_true <- runif(N, 0.45,0.55)#rep(0.5, length.out = N)#seq(0.01,0.99, length.out = N)
probs_computed <-  probs_true + if_else(runif(N) < 0.5, -1, 1) * 0.2

#K <- 7; probs_true <- c(rep(0.8, K), probs_computed[K:N])
#probs_true <-  plogis(qlogis(probs_computed) + 2 * if_else(probs_computed < 0.5, 1, - 1))

#probs_true <- c(seq(0.01,0.5, length.out = N/2),  rep(0.5, length.out = N/2))
# probs_computed <- probs_true + c(rep(0.1, N/2), rep(0, N/2))


true_model <- rbinom(N, size = 1, prob = probs_true)

N_post <- 99
post_draws <- rbinom(N, size = N_post, prob = probs_computed)
ranks <- numeric(N)
ranks[true_model == 0] <- SBC:::rdunif(sum(true_model == 0), a = 0, N_post - post_draws[true_model == 0])
ranks[true_model == 1] <- SBC:::rdunif(sum(true_model == 1), a = N_post - post_draws[true_model == 1], b = N_post)

ranks_df_example <- data.frame(sim_id = 1:N, variable = "model", rank = ranks, max_rank = N_post)

lgamma_thres <- log(SBC:::adjust_gamma(N = N, L = 1, K = N_post + 1))
lgamma_stat <- log_gamma_stat(ranks_df_example$rank, max_rank = N_post)

plot_rank_hist(ranks_df_example) | plot_ecdf_diff(ranks_df_example) + ggtitle(0.05 * exp(-lgamma_thres + lgamma_stat)) + theme(legend.position = "bottom") |
  my_reliability_diag(data.frame(prob = probs_computed, simulated_value = true_model)) + ggtitle(miscalibration_resampling_p(probs_computed, true_model, B = 10000)
  )

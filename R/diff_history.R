compute_diff_history_norm <- function(ranks, max_rank, p = 2) {
  rank_t <- rep(0, max_rank + 1)
  history_norm <- rep(NA_real_, length(ranks))
  uniform_cdf <- seq(0, 1, length.out = max_rank + 2)[2:(max_rank + 2)]
  for(i in 1:length(ranks)) {
    rank_t[ranks[i] + 1] <- rank_t[ranks[i] + 1] + 1
    diff <- uniform_cdf - cumsum(rank_t) / i
    history_norm[i] <- sum((abs(diff)^p)/(max_rank + 1))^(1/p)
  }
  history_norm
}

get_precomputed_gamma_thresholds <- function(K, min_sims = 1, max_sims = 1000) {
  if(!dir.exists("./_SBC_cache")) {
    dir.create("./_SBC_cache")
  }

  gamma_file <- paste0("./_SBC_cache/gamma_thresholds_",K,"_",min_sims,"_", max_sims, ".rds")
  if(!file.exists(gamma_file)) {
    message("Precomputing gamma thresholds")
    if(max_sims > 1000) {
      sims_above_1000 <- round(exp(seq(log(1001), log(max_sims), by = 0.02)))
      #Ensure max_sims is present
      sims_above_1000 <- unique(c(sims_above_1000, max_sims))
      if(min_sims <= 1000) {
        sims_to_compute <- c(min_sims:1000, sims_above_1000)
      } else {
        sims_to_compute <- sims_above_1000
      }
      interpolate <- TRUE
    } else {
      sims_to_compute <- min_sims:max_sims
      interpolate <- FALSE
    }
    # Computing time differs by N, randomly permute to get optimal division of work
    sims_to_compute <- sample(sims_to_compute)
    thresholds <- future.apply::future_sapply(sims_to_compute, FUN = function(n_sims) {
      SBC:::adjust_gamma(N = n_sims, L = 1, K = K)
    })


    if(interpolate) {
      thres_approx <- approx(log(sims_to_compute), log(thresholds), xout = log(min_sims:max_sims))
      thresholds_df <- data.frame(N_sims = min_sims:max_sims, log_gamma_threshold = thres_approx$y)
    } else {
      thresholds_df <-
        dplyr::arrange(data.frame(N_sims = sims_to_compute, log_gamma_threshold = log(thresholds)), N_sims)
    }

    saveRDS(thresholds_df, gamma_file)
  } else {
    thresholds_df <- readRDS(gamma_file)
  }
  thresholds_df
}

rank_counts_from_ranks <- function(ranks, max_rank) {
  stopifnot(all(ranks >= 0 & ranks <= max_rank))
  ranks_inflated_by_1 <- table(c(0:max_rank, ranks))
  return(as.numeric(ranks_inflated_by_1 - 1))
}


log_gamma_stat <- function(rank_counts, ranks_to_check = NULL) {
  max_rank <- length(rank_counts) - 1
  if(is.null(ranks_to_check)) {
    rank_ids_to_check = 1:max_rank
  } else {
    stopifnot(all(ranks_to_check >= 0 & ranks_to_check <= max_rank))
    rank_ids_to_check = setdiff(unique(ranks_to_check), max_rank) + 1
  }
  # Note that the ECDF follows the notation in Säilynoja, Bürkner & Vehtari 2022
  # i.e. it is the ECDF of a uniform distribution over j/(max_rank + 1) where j = 1:max_rank + 1
  scaled_ecdf <- cumsum(rank_counts)[rank_ids_to_check]
  z = rank_ids_to_check / (max_rank + 1)

  N_trials <- sum(rank_counts)
  return(log(2) + min(
    pbinom(scaled_ecdf, N_trials, z, log = TRUE),
    pbinom(scaled_ecdf - 1, N_trials, z, lower.tail = FALSE, log = TRUE)
  )
  )

}

compute_log_gamma_history <- function(ranks, max_rank) {
  rank_t <- rep(0, max_rank + 1)
  log_gamma <- rep(NA_real_, length(ranks))
  dummy <- rep(NA_real_, length(ranks))

  for(i in 1:length(ranks)) {
    rank_t[ranks[i] + 1] <- rank_t[ranks[i] + 1] + 1
    log_gamma[i] <- log_gamma_stat(rank_t)
    #dummy[i] <- log(2) + pbinom(scaled_ecdf[max_rank + 1] - 1, i, z[max_rank + 1], lower.tail = FALSE, log = TRUE)
  }
  #print(dummy)
  log_gamma
}

compute_ks_test_history <- function(ranks, max_rank) {
  ranks_cont <- (ranks + runif(length(ranks))) / (max_rank + 1)
  ks_p <- rep(NA_real_, length(ranks))
  for(i in 1:length(ranks)) {
    ks_p[i] <- ks.test(ranks_cont[1:i], "punif")$p.value
  }
  ks_p
}


compute_log_gamma_q_history <- function(pm1, true_model, K = 101) {
  log_gamma <- rep(NA_real_, length(pm1))

  z = unique(pm1)

  gamma_thresholds_df <- get_precomputed_gamma_thresholds(K = K,
                                                          min_sims = 1, max_sims = max(length(pm1), 1000))

  for(i in 1:length(pm1)) {
    stopifnot(gamma_thresholds_df$N_sims[i] == i)
    q <- compute_q(data.frame(pm1 = pm1[1:i], true_model = true_model[1:i]), include_01 = FALSE, K = K,
                   gamma = exp(gamma_thresholds_df$log_gamma_threshold[i]))


    log_gamma[i] <- log(2) + min(
      pbinom(q$q_sum, i, q$q_x, log = TRUE),
      pbinom(q$q_sum - 1, i, q$q_x, lower.tail = FALSE, log = TRUE)
    )
  }
  log_gamma
}

plot_log_gamma_history <- function(stats, min_sim_id = 0, max_sim_id = Inf, wrap_cols = 4, variables_regex = NULL, ylim = NULL) {
  unique_max_rank <- unique(stats$max_rank)
  if(length(unique_max_rank) > 1) {
    stop("Requires all max_rank to be equal")
  }

  max_sim_id_to_show <- min(max_sim_id, max(stats$sim_id))

  gamma_thresholds_df <- get_precomputed_gamma_thresholds(K = unique_max_rank + 1,
                                                          min_sims = 1, max_sims = max(max_sim_id_to_show, 1000))

  stats <- stats
  if(!is.null(variables_regex)) {
     stats <- stats %>% filter(grepl(variables_regex, variable))
  }

  stats %>%
    filter(sim_id <= max_sim_id) %>%
    group_by(variable) %>%
    mutate(log_gamma = compute_log_gamma_history(rank, unique_max_rank)) %>%
    filter(sim_id >= min_sim_id) %>%
    inner_join(gamma_thresholds_df, by = c("sim_id" = "N_sims")) %>%
    ggplot(aes(x = sim_id, y = log_gamma - log_gamma_threshold)) +
    geom_hline(yintercept = 0, color = "lightblue") +
    geom_line() +
    scale_y_continuous("Log Gamma - Threshold", limits = ylim) +
    scale_x_continuous("Number of simulations") +
    facet_wrap(~variable, ncol = wrap_cols)
}

plot_ks_test_history <- function(stats, min_sim_id = 0, max_sim_id = Inf, wrap_cols = 4, ylim = NULL) {
  stats %>%
    filter(sim_id <= max_sim_id) %>%
    group_by(variable) %>%
    mutate(ks_p = compute_ks_test_history(rank, unique(max_rank))) %>%
    filter(sim_id >= min_sim_id) %>%
    ggplot(aes(x = sim_id, y = ks_p)) +
    geom_hline(yintercept = 0.05,  color = "lightblue") +
    geom_line() +
    scale_y_log10("P - value (KS test)", limits = ylim) +
    scale_x_continuous("Number of simulations") +
    facet_wrap(~variable, ncol = wrap_cols)
}


plot_log_gamma_q_history <- function(stats, min_sim_id = 0, max_sim_id = Inf, wrap_cols = 4, variables_regex = NULL, ylim = NULL, K = 101) {
  max_sim_id_to_show <- min(max_sim_id, max(stats$sim_id))

  # TODO: this is only approximate
  gamma_thresholds_df <- get_precomputed_gamma_thresholds(K = K,
                                                          min_sims = 1, max_sims = max(max_sim_id_to_show, 1000))


  if(!is.null(variables_regex)) {
    stats <- stats %>% filter(grepl(variables_regex, variable))
  }

  stats %>%
    filter(sim_id <= max_sim_id) %>%
    group_by(variable) %>%
    mutate(log_gamma_q = compute_log_gamma_q_history(pm1, true_model)) %>%
    filter(sim_id >= min_sim_id) %>%
    inner_join(gamma_thresholds_df, by = c("sim_id" = "N_sims")) %>%
    ggplot(aes(x = sim_id, y = log_gamma_q - log_gamma_threshold)) +
    geom_hline(yintercept = 0, color = "lightblue") +
    geom_line() +
    scale_y_continuous("Log Gamma - Threshold", limits = ylim) +
    scale_x_continuous("Number of simulations") +
    facet_wrap(~variable, ncol = wrap_cols)
}

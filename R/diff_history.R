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

include_step <- function(sim_id, step) {
  sim_id %% step == 0
}

step_id <- function(sim_id, step) {
  sim_id[include_step(sim_id, step)] / step
}

compute_log_gamma_history_single <- function(ranks, max_rank, step = 1) {
  rank_t <- rep(0, max_rank + 1)
  log_gamma <- rep(NA_real_, floor(length(ranks) / step))

  for(i in 1:length(ranks)) {
    rank_t[ranks[i] + 1] <- rank_t[ranks[i] + 1] + 1
    if(include_step(i, step)) {
      log_gamma[step_id(i, step)] <- log_gamma_stat(rank_t)
    }
  }
  log_gamma
}

compute_log_gamma_history <- function(stats, step = step, min_sim_id = 1) {
  unique_max_rank <- unique(stats$max_rank)
  if(length(unique_max_rank) > 1) {
    stop("Requires all max_rank to be equal")
  }

  max_n_sims <- stats %>% group_by(variable) %>% summarise(n = n()) %>% pull(n) %>% max()
  gamma_thresholds_df <- get_precomputed_gamma_thresholds(K = unique_max_rank + 1,
                                                          min_sims = 1, max_sims = max(max_n_sims, 1000))
  stats %>%
    group_by(variable) %>%
    reframe(
      sim_id = sim_id[include_step(sim_id, step)],
      log_gamma = compute_log_gamma_history_single(rank, unique_max_rank, step = step)
      ) %>%
    filter(sim_id >= min_sim_id) %>%
    inner_join(gamma_thresholds_df, by = c("sim_id" = "N_sims"))
}

compute_schad_history_single <- function(probs, step = 1) {
  steps_to_include <- include_step(1:length(probs), step)
  sums <- cumsum(probs)[steps_to_include]
  ns <- (1:length(probs))[steps_to_include]

  sds <- numeric(floor(length(probs)/step))
  for(i in which(steps_to_include)) {
    sds[step_id(i,step)] <- sd(probs[1:i])
  }
  t_stat <- ((sums/ns) - 0.5) / (sds / sqrt(ns))

  dfs <- ns - 1
  log_ps <- log(2) + pt(-abs(t_stat), dfs, log.p = TRUE)

  return(log_ps)
}

compute_schad_history <- function(stats, step = 1, min_sim_id = 1) {
  if(!("prob" %in% names(stats))) {
    stop("Stats must contain prob - maybe you forgot to call `binary_probabilities_from_stats`?")
  }
  stats %>%
    group_by(variable) %>%
    reframe(
      sim_id = sim_id[include_step(sim_id, step)],
      log_p = compute_schad_history_single(prob, step = step)
    ) %>%
    filter(sim_id >= min_sim_id)
}

compute_giviti_history_single <- function(probs, true_model, step = 1) {
  steps_to_include <- which(include_step(1:length(probs), step))

  ps <- numeric(length(steps_to_include))
  for(i in 1:length(steps_to_include)) {
    probs_to_test <- probs[1:steps_to_include[i]]
    true_to_test <- true_model[1:steps_to_include[i]]
    if(length(true_to_test) < 5 || all(true_to_test == 0) || all(true_to_test == 1)){
      ps[i] <- 1
    } else {
      print(true_to_test)
      givi <- givitiR::givitiCalibrationTest(true_to_test, probs_to_test, devel = "external")
      ps[i] <- givi$p.value
    }
  }
  return(log(ps))
}

compute_giviti_history <- function(stats, step = 1, min_sim_id = 1) {
  if(!("prob" %in% names(stats))) {
    stop("Stats must contain prob - maybe you forgot to call `binary_probabilities_from_stats`?")
  }
  stats %>%
    group_by(variable) %>%
    reframe(
      sim_id = sim_id[include_step(sim_id, step)],
      log_p = compute_giviti_history_single(prob, simulated_value, step = step)
    ) %>%
    filter(sim_id >= min_sim_id)
}

compute_bootstrapped_histories <- function(stats, history_length, n_histories, history_fun, step = 1, min_sim_id = 1) {
  res_df <- list()
  for(h in 1:n_histories) {
    stats_boot <- stats %>%
      group_by(variable) %>%
      sample_n(history_length) %>%
      mutate(sim_id = 1:n()) %>%
      ungroup()
    res_df[[h]] <- history_fun(stats_boot, step = step, min_sim_id = min_sim_id)
    res_df[[h]]$history_id = h
  }
  do.call(rbind, res_df)
}


plot_log_gamma_histories <- function(histories_df, min_sim_id = 0, wrap_cols = 4, variables_regex = NULL, ylim = NULL) {
  alpha <- sqrt(1/length(unique(histories_df$history_id)))

  histories_df %>%
    filter(sim_id >= min_sim_id) %>%
    ggplot(aes(x = sim_id, y = log_gamma - log_gamma_threshold, group = history_id)) +
    geom_line(alpha = alpha) +
    geom_hline(yintercept = 0, color = "lightblue", linewidth = 1) +
    scale_y_continuous("Log Gamma - Threshold", limits = ylim) +
    scale_x_continuous("Number of simulations") +
    facet_wrap(~variable, ncol = wrap_cols)
}

plot_log_p_histories <- function(histories_df, title, min_sim_id = 0, wrap_cols = 4, variables_regex = NULL, ylim = NULL) {
  alpha <- sqrt(1/length(unique(histories_df$history_id)))

  histories_df %>%
    filter(sim_id >= min_sim_id) %>%
    ggplot(aes(x = sim_id, y = log_p, group = history_id)) +
    geom_line(alpha = alpha) +
    geom_hline(yintercept = log(0.05), color = "lightblue", linewidth = 1) +
    scale_y_continuous(paste0("Log p (", title,")"), limits = ylim) +
    scale_x_continuous("Number of simulations") +
    facet_wrap(~variable, ncol = wrap_cols)
}

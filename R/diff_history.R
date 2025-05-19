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

#' @export
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

#' @export
log_gamma_stat <- function(ranks, max_rank, ranks_to_check = NULL) {
  rank_counts <- rep(0, max_rank + 1)
  for(i in 0:max_rank) {
    rank_counts[i + 1] <- sum(ranks == i)
  }
  stopifnot(sum(rank_counts) == length(ranks))
  log_gamma_stat_counts(rank_counts, ranks_to_check)
}

#' @export
log_gamma_stat_counts <- function(rank_counts, ranks_to_check = NULL) {
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
      log_gamma[step_id(i, step)] <- log_gamma_stat_counts(rank_t)
    }
  }
  log_gamma
}

#' @export
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

compute_schad_history_single <- function(probs, true_model, step = 1, expected = 0.5) {
  steps_to_include <- include_step(1:length(probs), step)
  sums <- cumsum(probs)[steps_to_include]
  ns <- (1:length(probs))[steps_to_include]

  sds <- numeric(floor(length(probs)/step))
  for(i in which(steps_to_include)) {
    sds[step_id(i,step)] <- sd(probs[1:i])
  }
  means <- sums/ns
  t_stat <- (means - expected) / (sds / sqrt(ns))

  dfs <- ns - 1
  log_ps <- case_when(
    sds == 0 & (abs(means - expected) < 1e-8) ~ 1,
    sds == 0 ~ NA_real_,
    TRUE ~ log(2) + pt(-abs(t_stat), dfs, log.p = TRUE)
  )

  return(log_ps)
}

#' @export
compute_schad_history <- function(...) {
  compute_calibration_history(compute_schad_history_single, ...)
}

compute_schad_bf_history_single <- function(probs, true_model, step = 1, expected = 0.5) {
  steps_to_include <- include_step(1:length(probs), step)

  log_ps <- numeric(floor(length(probs)/step))
  uniform_sd <- sqrt(1/12) # Assume alternative matching the variance of the uniform distribution
  for(i in which(steps_to_include)) {
    # Handling some edge cases typically with few sims where the Bayes factors
    # don't play nicely
    probs_to_test <- probs[1:i]
    out_i <- step_id(i,step)
    if(i <= 5) {
      log_ps[out_i] <- 0
    } else if(sd(probs_to_test) == 0) {
      if(abs(mean(probs_to_test) - expected) < 1e-8) {
        log_ps[out_i] <- 0
      } else {
        log_ps[out_i] <- NA_real_
      }
    } else if(sd(probs_to_test) < 0.0001) {
      if(abs(mean(probs_to_test) - expected) < 0.001) {
        log_ps[out_i] <- 0
      } else {
        log_ps[out_i] <- pnorm(-abs(mean(probs_to_test) - expected), 0, 0.001, log.p = TRUE) + log(2)
      }
    } else {
      bf <- tryCatch(BayesFactor::ttestBF(probs_to_test, mu = expected, rscale = uniform_sd),
                     error = \(e) { message("ttestBF failed for c(",paste0(probs_to_test, collapse = ", "), ")");  NA_real_ })
      log_ps[out_i] <-  plogis(-bf@bayesFactor[1, "bf"], log.p = TRUE)
    }
  }

  return(log_ps)
}

#' @export
compute_schad_bf_history <- function(...) {
  compute_calibration_history(compute_schad_bf_history_single, ...)
}

compute_ttest_history_single <- function(probs, true_model, step = 1, expected = 0.5) {
  steps_to_include <- include_step(1:length(probs), step)

  stopifnot(is.numeric(expected) || expected == "avg_true")
  log_ps <- numeric(floor(length(probs)/step))
  for(i in which(steps_to_include)) {
    # Handling some edge cases typically with few sims where t-test
    # don't play nicely
    probs_to_test <- probs[1:i]
    out_i <- step_id(i,step)
    if(i < 2) {
      log_ps[out_i] <- 0
    } else if(sd(probs_to_test) == 0 && expected != "avg_true") {
      if(abs(mean(probs_to_test) - expected) < 1e-8) {
        log_ps[out_i] <- 0
      } else {
        log_ps[out_i] <- NA_real_
      }
    } else if(sd(probs_to_test) == 0 && expected == "avg_true" && sd(true_model[1:i]) == 0) {
      if(abs(mean(probs_to_test) - mean(true_model[1:i])) < 1e-8) {
        log_ps[out_i] <- 0
      } else {
        log_ps[out_i] <- NA_real_
      }

    } else {
      if(expected == "avg_true") {
        t <- t.test(probs_to_test, true_model[1:i])
      } else {
        t <- t.test(probs_to_test, mu = expected)
      }
      log_ps[out_i] <- log(t$p.value)
    }
  }

  return(log_ps)
}

#' @export
compute_ttest_history <- function(...) {
  compute_calibration_history(compute_ttest_history_single, ...)
}


compute_giviti_history_single <- function(probs, true_model, step = 1) {
  steps_to_include <- which(include_step(1:length(probs), step))

  ps <- numeric(length(steps_to_include))
  for(i in 1:length(steps_to_include)) {
    probs_to_test <- probs[1:steps_to_include[i]]
    true_to_test <- true_model[1:steps_to_include[i]]
    if(length(true_to_test) < 10 || all(true_to_test == 0) || all(true_to_test == 1)){
      ps[i] <- 1
    } else {
      ps[i] <- tryCatch({
        givi <- givitiR::givitiCalibrationTest(true_to_test, probs_to_test, devel = "external")
        givi$p.value
      }, error = \(e) {
        message("giviti failed for c(",paste0(probs[1:i], collapse = ", "), "), c(",paste0(true_model[1:i], collapse = ", "), ")");
        #message(e)
        NA_real_
      })
    }
  }
  return(log(ps))
}


#' @export
compute_giviti_history <- function(...) {
  compute_calibration_history(compute_giviti_history_single, ...)
}

compute_dimitriadis_history_single <- function(probs, true_model, step = 1, saddle = TRUE, adjust.method = "holm") {
  steps_to_include <- which(include_step(1:length(probs), step))

  ps <- numeric(length(steps_to_include))
  for(i in 1:length(steps_to_include)) {
    probs_to_test <- probs[1:steps_to_include[i]]
    true_to_test <- true_model[1:steps_to_include[i]]
    if(saddle) {
      ps[i] <- calibration_p_saddle(probs_to_test, true_to_test, adjust.method = adjust.method)
    } else {
      ps[i] <- calibration_p(probs_to_test, true_to_test, adjust.method = adjust.method)
    }
  }
  return(log(ps))
}

#' @export
compute_dimitriadis_history <- function(...) {
  compute_calibration_history(compute_dimitriadis_history_single, ...)
}


#' @export
compute_miscalibration_history <- function(...) {
  single_func <- function(...) {
    compute_indep_history_single(miscalibration_resampling_p, ...)
  }
  compute_calibration_history(single_func, ...)
}

#' @export
compute_hosmer_lemeshow_history <- function(...) {
  single_func <- function(...) {
    compute_indep_history_single(hosmer_lemeshow_p, ...)
  }
  compute_calibration_history(single_func, ...)
}

compute_brier_history_single <- function(probs, true_model, step = 1, B = 2000) {
  hist_all <- brier_resampling_history(probs, true_model, B = B)
  hist <- hist_all[include_step(1:length(probs), step)]
  hist[hist == 0] <- 0.5/B
  return(log(hist))
}

#' @export
compute_brier_history <- function(...) {
  compute_calibration_history(compute_brier_history_single, ...)
}


#' @export
compute_indep_history_single <- function(func, probs, true_model, step = 1, ...) {
  steps_to_include <- which(include_step(1:length(probs), step))

  ps <- numeric(length(steps_to_include))
  for(i in 1:length(steps_to_include)) {
    probs_to_test <- probs[1:steps_to_include[i]]
    true_to_test <- true_model[1:steps_to_include[i]]
    ps[i] <- func(probs_to_test, true_to_test, ...)
  }
  return(log(ps))
}

#' @export
compute_calibration_history <- function(history_single_func, stats, step = 1, min_sim_id = 1, ...) {
  if(!("prob" %in% names(stats))) {
    stop("Stats must contain prob - maybe you forgot to call `binary_probabilities_from_stats`?")
  }
  stats |>
    dplyr::group_by(variable) |>
    dplyr::reframe(
      sim_id = sim_id[include_step(sim_id, step)],
      log_p = history_single_func(prob, simulated_value, step = step, ...)
    ) |>
    dplyr::filter(sim_id >= min_sim_id)
}




#' @export
compute_bootstrapped_histories <- function(stats, history_length, n_histories, history_fun, step = 1, min_sim_id = 1, cache = TRUE, ...) {
  if(cache) {
    cache_dir <- here::here("cache", "hist_single")
    if(!dir.exists(cache_dir)) {
      dir.create(cache_dir)
    }
    dotlist <- list(...)
    if(length(dotlist) == 0) {
      dotlist_hash <- ""
    } else {
      dotlist_hash <- paste0("_", rlang::hash(dotlist))
    }
    cache_basename <- paste0(rlang::as_name(enquo(history_fun)), "_len", history_length, "_nhist",
                             n_histories, "_step", step, "_min", min_sim_id,
                             rlang::hash(stats), dotlist_hash, ".rds")
    cache_file <- file.path(cache_dir, cache_basename)
    if(file.exists(cache_file)) {
      message("Read from cache '", cache_file, "'")
      return(readRDS(cache_file))
    }
  }
  res_list <- future.apply::future_lapply(1:n_histories, future.seed = TRUE,
                                          FUN = function(history_id) {
      stats_boot <- dplyr::group_by(stats, variable) |>
        dplyr::sample_n(history_length) |>
        dplyr::mutate(sim_id = 1:dplyr::n()) |>
        dplyr::ungroup()
      res_df <- history_fun(stats_boot, step = step, min_sim_id = min_sim_id, ...)
      res_df$history_id = history_id
      return(res_df)
  })
  res <- do.call(rbind, res_list)

  if(cache) {
    saveRDS(res, cache_file)
  }

  return(res)
  # res_df <- list()
  # for(h in 1:n_histories) {
  #   stats_boot <- stats %>%
  #     group_by(variable) %>%
  #     sample_n(history_length) %>%
  #     mutate(sim_id = 1:n()) %>%
  #     ungroup()
  #   res_df[[h]] <- history_fun(stats_boot, step = step, min_sim_id = min_sim_id, ...)
  #   res_df[[h]]$history_id = h
  # }
  # do.call(rbind, res_df)
}


plot_log_gamma_histories <- function(histories_df, min_sim_id = 0, wrap_cols = 4, variables_regex = NULL, ylim = NULL) {
  alpha <- sqrt(1/length(unique(histories_df$history_id)))

  variable_order <- as.integer(factor(histories_df$variable))
  variable_order[histories_df$variable == "model"] <- -1
  histories_df$variable <- forcats::fct_reorder(histories_df$variable, variable_order)

  power_df <- histories_df %>% group_by(sim_id, variable) %>%
    summarise(power = mean(log_gamma < log_gamma_threshold), .groups = "drop") %>%
    filter(power >= 0.8) %>%
    group_by(variable) %>%
    summarise(first_power_sim_id = min(c(Inf,sim_id)))


  histories_df %>%
    filter(sim_id >= min_sim_id) %>%
    ggplot(aes(x = sim_id, y = log_gamma - log_gamma_threshold, group = history_id)) +
    geom_vline(aes(xintercept = first_power_sim_id), data = power_df, color = "orangered") +
    geom_line(alpha = alpha) +
    geom_hline(yintercept = 0, color = "lightblue", linewidth = 1) +
    scale_y_continuous("Log Gamma - Threshold", limits = ylim) +
    scale_x_continuous("Number of simulations") +
    facet_wrap(~variable, ncol = wrap_cols)
}

log_p_breaks <- function(lim) {
  min_log10 <- min(lim) / log(10)
  if(min_log10 > -2) {
    return(log(c(0.5, 0.05,0.01)))
  }
  n_extra_breaks <- 4
  step <- floor(min_log10 / n_extra_breaks)
  extra_breaks_log10 <- -1 + step * (1:n_extra_breaks)
  return(c(log(0.05), extra_breaks_log10 * log(10)))
}

log_p_labels <- function(breaks) {
  signif(exp(breaks), digits = 3)
}

plot_log_p_histories <- function(histories_df, title = NULL, min_sim_id = 0, wrap_cols = 4, variables_regex = NULL, ylim = NULL) {
  alpha <- sqrt(1/length(unique(histories_df$history_id)))

  if(is.null(title)) {
    y_label <- "p-value"
  } else {
    y_label <- paste0("p - ", title)
  }


  if(length(unique(histories_df$variable)) == 1) {
    histories_df$variable <- " "
  } else {
    variable_order <- as.integer(factor(histories_df$variable))
    variable_order[histories_df$variable == "model"] <- -1
    histories_df$variable <- forcats::fct_reorder(histories_df$variable, variable_order)
  }

  power_df <- histories_df %>% group_by(sim_id, variable) %>%
    summarise(power = mean(log_p < log(0.05)), .groups = "drop") %>%
    filter(power >= 0.8) %>%
    group_by(variable) %>%
    summarise(first_power_sim_id = min(c(Inf,sim_id)))


  histories_df %>%
    filter(sim_id >= min_sim_id) %>%
    ggplot(aes(x = sim_id, y = log_p, group = history_id)) +
    geom_vline(aes(xintercept = first_power_sim_id), data = power_df, color = "orangered") +
    geom_line(alpha = alpha) +
    geom_hline(yintercept = log(0.05), color = "lightblue", linewidth = 1) +
    scale_y_continuous(y_label, limits = ylim, breaks = log_p_breaks, labels = log_p_labels) +
    scale_x_continuous("Number of simulations") +
    facet_wrap(~variable, ncol = wrap_cols)
}

save_histories <- function(name, ...) {
  saveRDS(list(...), here::here("cache", paste0("hist_", name, ".rds")))
}

load_precomputed_file <- function(filename, producer_script) {
  if(!file.exists(filename)) {
    stop("File: `", filename, "` does not exist. Run ", producer_script, " to generate it.")
  }
  readRDS(filename)
}

load_histories <- function(name, producer_script) {
  filename <- here::here("cache", paste0("hist_", name, ".rds"))
  load_precomputed_file(filename, producer_script)
}

